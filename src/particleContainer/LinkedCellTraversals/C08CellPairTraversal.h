/*
 * C08CellPairTraversal.h
 *
 *  Created on: 15 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08CELLPAIRTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08CELLPAIRTRAVERSAL_H_

#include "particleContainer/LinkedCellTraversals/C08BasedTraversals.h"
#include "utils/threeDimensionalMapping.h"
#include "utils/mardyn_assert.h"


template <class CellTemplate>
class C08CellPairTraversal: public C08BasedTraversals<CellTemplate> {
public:
	C08CellPairTraversal(
			std::vector<CellTemplate>& cells,
			const std::array<unsigned long, 3>& dims) :
			C08BasedTraversals<CellTemplate>(cells, dims) {
	}
	~C08CellPairTraversal() {
	}

	void traverseCellPairs(CellProcessor& cellProcessor) const;
	void traverseCellPairsOuter(CellProcessor& cellProcessor) const;
	void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) const;

protected:
	void traverseCellPairsBackend(CellProcessor& cellProcessor,
			const std::array<unsigned long, 3> & start,
			const std::array<unsigned long, 3> & end,
			const std::array<unsigned long, 3> & stride) const;
};


template<class CellTemplate>
void C08CellPairTraversal<CellTemplate>::traverseCellPairs(
		CellProcessor& cellProcessor) const {

	using std::array;
	const array<unsigned long, 3> strides = { 2, 2, 2 };
	array<unsigned long, 3> end;
	for (int d = 0; d < 3; ++d) {
		end[d] = this->_dims[d] - 1;
	}

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		for (unsigned long col = 0; col < 8; ++col) {
			std::array<unsigned long, 3> begin = threeDimensionalMapping::oneToThreeD(col, strides);

			traverseCellPairsBackend(cellProcessor, begin, end, strides);
			#if defined(_OPENMP)
			#pragma omp barrier
			#endif
			// this barrier is needed, since we have a nowait in the backend
		}
	}
}

template<class CellTemplate>
void C08CellPairTraversal<CellTemplate>::traverseCellPairsOuter(
		CellProcessor& cellProcessor) const {
	using std::array;
	const array<unsigned long, 3> strides2 = { 2, 2, 2 };

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		for (unsigned long col = 0; col < 8; ++col) {
			// halo & boundaries in z direction
			const array< unsigned long, 3> begin = threeDimensionalMapping::oneToThreeD(col, strides2);

			// values, which are modified by 2 are actually modified by the respective stride, which is always 2
			array<unsigned long, 3> startZ = { begin[0], begin[1], begin[2] }; 									// default
			array<unsigned long, 3> endZ = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 1 };		// default
			array<unsigned long, 3> stridesZ = {strides2[0], strides2[1], this->_dims[2] - 3};					// mod z
			traverseCellPairsBackend(cellProcessor, startZ, endZ, stridesZ);

			// halo & boundaries in y direction
			// boundaries in z direction are excluded!
			array<unsigned long, 3> startY = { begin[0], begin[1], begin[2] + 2 };								// mod z
			array<unsigned long, 3> endY = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 3 };		// mod z
			array<unsigned long, 3> stridesY = {strides2[0], this->_dims[1] - 3, strides2[2]};					// mod y
			traverseCellPairsBackend(cellProcessor, startY, endY, stridesY);

			// halo & boundaries in x direction
			// boundaries in z and y direction are excluded!
			array<unsigned long, 3> startX = { begin[0], begin[1] + 2, begin[2] + 2 };							// mod yz
			array<unsigned long, 3> endX = { this->_dims[0] - 1, this->_dims[1] - 3, this->_dims[2] - 3 };		// mod yz
			array<unsigned long, 3> stridesX = {this->_dims[0] - 3, strides2[1], strides2[2]};					// mod x
			traverseCellPairsBackend(cellProcessor, startX, endX, stridesX);

			#if defined(_OPENMP)
				// this barrier is needed, since we have a nowait in the backend
				// except at the last run - there we have the barrier at the end of the parallel region
				if (col < 7) {
					#pragma omp barrier
				}
			#endif
		}
	} // end pragma omp parallel
}

template<class CellTemplate>
void C08CellPairTraversal<CellTemplate>::traverseCellPairsInner(
		CellProcessor& cellProcessor, unsigned stage,
		unsigned stageCount) const {
	using std::array;

	unsigned long splitdim = 0;
	unsigned long maxcellsize = this->_dims[0];
	for (unsigned long i = 1; i < 3; i++){
		if (this->_dims[i] > maxcellsize) {
			splitdim = i;
			maxcellsize = this->_dims[i];
		}
	}
	unsigned long splitsize = maxcellsize - 5;
	unsigned long minsize = min(this->_dims[0], min(this->_dims[1], this->_dims[2]));

	mardyn_assert(minsize >= 4);  // there should be at least 4 cells in each dimension, otherwise we did something stupid!

	if (minsize <= 5) {
		return;  // we can not iterate over any inner cells, that do not depend on boundary or halo cells
	}

	array<unsigned long, 3> lower;
	array<unsigned long, 3> upper;
	for (unsigned long i = 0; i < 3; i++) {
		lower[i] = 2;
		upper[i] = this->_dims[i] - 3;
	}
	lower[splitdim] = 2 + splitsize * stage / stageCount;  // at least 2
	upper[splitdim] = 2 + splitsize * (stage+1) / stageCount;  // at most _cellsPerDimension[i] - 3


	#if defined(_OPENMP)
		#pragma omp parallel
	#endif
	{
		array<unsigned long, 3> strides = {2, 2, 2};

		for (unsigned long col = 0; col < 8; ++col) {
			array<unsigned long, 3> startIndices = threeDimensionalMapping::oneToThreeD(col, strides);
			for (int i = 0; i < 3; i++) {
				startIndices[i] = startIndices[i] + lower[i];
			}

			traverseCellPairsBackend(cellProcessor, startIndices, upper, strides);
			#if defined(_OPENMP)
				// this barrier is needed, since we have a nowait in the backend
				// except at the last run - there we have the barrier at the end of the parallel region
				if (col < 7) {
					#pragma omp barrier
				}
			#endif
		}
	} // end pragma omp parallel
}


template<class CellTemplate>
void C08CellPairTraversal<CellTemplate>::traverseCellPairsBackend(
		CellProcessor& cellProcessor,
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end,
		const std::array<unsigned long, 3>& stride) const {

	// note parallel region is open outside

	#if defined(_OPENMP)
	#pragma omp for schedule(dynamic, 1) collapse(3) nowait
	#endif
	for (unsigned long z = start[2]; z < end[2]; z += stride[2]) {
		for (unsigned long y = start[1]; y < end[1]; y += stride[1]) {
			for (unsigned long x = start[0]; x < end[0]; x += stride[0]) {
				unsigned long baseIndex = threeDimensionalMapping::threeToOneD( x, y, z, this->_dims);
				this->processBaseCell(cellProcessor, baseIndex);
			}
		}
	}
}


#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08CELLPAIRTRAVERSAL_H_ */