#include "particleContainer/LinkedCells.h"

#include <cmath>

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "ParticleCell.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include <array>
using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################

LinkedCells::LinkedCells(double bBoxMin[3], double bBoxMax[3],
		double cutoffRadius, double LJCutoffRadius, double cellsInCutoffRadius) :
		ParticleContainer(bBoxMin, bBoxMax) {
	int numberOfCells = 1;
	_cutoffRadius = cutoffRadius;
	_LJCutoffRadius = LJCutoffRadius;

	global_log->debug() << "cutoff: " << cutoffRadius << endl;
	global_log->debug() << "LJ cutoff:" << LJCutoffRadius << endl;
	global_log->debug() << "# cells in cutoff: " << cellsInCutoffRadius << endl;

	_cellsInCutoff = ceil(cellsInCutoffRadius);

	if (_cellsInCutoff != 1) {
		global_log->error()
				<< "With the recent release only 1 cell per cutoff radius is supported,"
				<< " but the input file prescribes " << _cellsInCutoff
				<< " cells per cutoff radius." << endl
				<< "\tThe support has been dropped, since no speedup can be expected using"
				<< " multiple cells per cutoff radius." << endl
				<< "\tIf you can provide a case, where this is not true, please contact us."
				<< endl;
		global_simulation->exit(-1);
	}

	for (int d = 0; d < 3; d++) {
		/* first calculate the cell length for this dimension */
		_boxWidthInNumCells[d] = floor(
				(_boundingBoxMax[d] - _boundingBoxMin[d]) / cutoffRadius
						* cellsInCutoffRadius);
		// in each dimension at least one layer of (inner+boundary) cells is necessary
		if (_boxWidthInNumCells[d] == 0) {
			_boxWidthInNumCells[d] = 1;
		}
		_cellLength[d] = (_boundingBoxMax[d] - _boundingBoxMin[d])
				/ _boxWidthInNumCells[d];
		_haloWidthInNumCells[d] = ceil(cellsInCutoffRadius);
		_haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloLength[d];

		_cellsPerDimension[d] = _boxWidthInNumCells[d]
				+ 2 * _haloWidthInNumCells[d];

		numberOfCells *= _cellsPerDimension[d];
		assert(numberOfCells > 0);
	}
	global_log->debug() << "Cell size (" << _cellLength[0] << ", "
			<< _cellLength[1] << ", " << _cellLength[2] << ")" << endl;
	global_log->debug() << "Cells per dimension (incl. halo): "
			<< _cellsPerDimension[0] << " x " << _cellsPerDimension[1] << " x "
			<< _cellsPerDimension[2] << endl;

	_cells.resize(numberOfCells);

	// If the width of the inner region is less than the width of the halo
	// region a parallelization is not possible (with the used algorithms).
	// If a particle leaves this box, it would need to be communicated to the two next neighbors.
	if (_boxWidthInNumCells[0] < 2 * _haloWidthInNumCells[0]
			|| _boxWidthInNumCells[1] < 2 * _haloWidthInNumCells[1]
			|| _boxWidthInNumCells[2] < 2 * _haloWidthInNumCells[2]) {
		global_log->error_always_output()
				<< "LinkedCells (constructor): bounding box too small for calculated cell length"
				<< endl;
		global_log->error_always_output() << "_cellsPerDimension: " << _cellsPerDimension[0]
				<< " / " << _cellsPerDimension[1] << " / "
				<< _cellsPerDimension[2] << endl;
		global_log->error_always_output() << "_haloWidthInNumCells: "
				<< _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1]
				<< " / " << _haloWidthInNumCells[2] << endl;
		global_log->error_always_output() << "_boxWidthInNumCells: " << _boxWidthInNumCells[0]
				<< " / " << _boxWidthInNumCells[1] << " / "
				<< _boxWidthInNumCells[2] << endl;
		global_simulation->exit(5);
	}

	initializeCells();
	calculateNeighbourIndices();
	calculateCellPairOffsets();

	_cellsValid = false;
}

LinkedCells::~LinkedCells() {
	std::vector<ParticleCell>::iterator it;
	for (it = _cells.begin(); it != _cells.end(); ++it) {
		it->deallocateAllParticles();
	}
}

void LinkedCells::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("cellsInCutoffRadius", _cellsInCutoff);
	if (_cellsInCutoff != 1) {
		global_log->error()
				<< "With the recent release only 1 cell per cutoff radius is supported,"
				<< " but the input file prescribes " << _cellsInCutoff
				<< " cells per cutoff radius." << endl
				<< "\tThe support has been dropped, since no speedup can be expected using"
				<< " multiple cells per cutoff radius." << endl
				<< "\tIf you can provide a case, where this is not true, please contact us."
				<< endl;
		global_simulation->exit(-1);
	}
	global_log->info() << "Cells in cut-off radius: " << _cellsInCutoff << endl;
}

void LinkedCells::rebuild(double bBoxMin[3], double bBoxMax[3]) {
	for (int i = 0; i < 3; i++) {
		this->_boundingBoxMin[i] = bBoxMin[i];
		this->_boundingBoxMax[i] = bBoxMax[i];
		_haloWidthInNumCells[i] = ::ceil(_cellsInCutoff); /* TODO: Single value?! */
	}
	int numberOfCells = 1;

	for (int dim = 0; dim < 3; dim++) {
		_boxWidthInNumCells[dim] = floor(
				(_boundingBoxMax[dim] - _boundingBoxMin[dim]) / _cutoffRadius
						* _haloWidthInNumCells[dim]);
		// in each dimension at least one layer of (inner+boundary) cells is necessary
		if (_boxWidthInNumCells[dim] == 0) {
			_boxWidthInNumCells[dim] = 1;
		}

		_cellsPerDimension[dim] = (int) floor(
				(_boundingBoxMax[dim] - _boundingBoxMin[dim])
						/ (_cutoffRadius / _haloWidthInNumCells[dim]))
				+ 2 * _haloWidthInNumCells[dim];
		// in each dimension at least one layer of (inner+boundary) cells necessary
		if (_cellsPerDimension[dim] == 2 * _haloWidthInNumCells[dim]) {
			global_log->error_always_output() << "LinkedCells::rebuild: region to small"
					<< endl;
			global_simulation->exit(1);
		}
		numberOfCells *= _cellsPerDimension[dim];
		_cellLength[dim] = (_boundingBoxMax[dim] - _boundingBoxMin[dim])
				/ (_cellsPerDimension[dim] - 2 * _haloWidthInNumCells[dim]);
		_haloBoundingBoxMin[dim] = this->_boundingBoxMin[dim]
				- _haloWidthInNumCells[dim] * _cellLength[dim];
		_haloBoundingBoxMax[dim] = this->_boundingBoxMax[dim]
				+ _haloWidthInNumCells[dim] * _cellLength[dim];
		_haloLength[dim] = _haloWidthInNumCells[dim] * _cellLength[dim];
	}

	_cells.resize(numberOfCells);

	// If the with of the inner region is less than the width of the halo region
	// a parallelisation isn't possible (with the used algorithms).
	// In this case, print an error message
	// _cellsPerDimension is 2 times the halo width + the inner width
	// so it has to be at least 3 times the halo width
	if (_boxWidthInNumCells[0] < 2 * _haloWidthInNumCells[0]
			|| _boxWidthInNumCells[1] < 2 * _haloWidthInNumCells[1]
			|| _boxWidthInNumCells[2] < 2 * _haloWidthInNumCells[2]) {
		global_log->error_always_output()
				<< "LinkedCells (rebuild): bounding box too small for calculated cell Length"
				<< endl;
		global_log->error_always_output() << "cellsPerDimension " << _cellsPerDimension[0]
				<< " / " << _cellsPerDimension[1] << " / "
				<< _cellsPerDimension[2] << endl;
		global_log->error_always_output() << "_haloWidthInNumCells" << _haloWidthInNumCells[0]
				<< " / " << _haloWidthInNumCells[1] << " / "
				<< _haloWidthInNumCells[2] << endl;
		global_simulation->exit(5);
	}

	initializeCells();
	calculateNeighbourIndices();
	calculateCellPairOffsets();

	// TODO: We loose particles here as they are not communicated to the new owner
	// delete all Particles which are outside of the halo region
	deleteParticlesOutsideBox(_haloBoundingBoxMin, _haloBoundingBoxMax);

	_cellsValid = false;
}

void LinkedCells::update() {
	const vector<ParticleCell>::size_type numCells = _cells.size();

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (vector<ParticleCell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++) {
			_cells[cellIndex].preUpdateLeavingMolecules();
		}

		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (vector<ParticleCell>::size_type cellIndex = 0; cellIndex < numCells; cellIndex++) {
			ParticleCell& cell = _cells[cellIndex];

			for (unsigned long j = 0; j < _backwardNeighbourOffsets.size(); j++) {
				const unsigned long neighbourIndex = cellIndex - _backwardNeighbourOffsets[j];
				if (neighbourIndex >= _cells.size()) {
					// handles cell_index < 0 (indices are unsigned!)
					assert(cell.isHaloCell());
					continue;
				}
				cell.updateLeavingMolecules(_cells[neighbourIndex]);
			}

			for (unsigned long j = 0; j < _forwardNeighbourOffsets.size(); j++) {
				const unsigned long neighbourIndex = cellIndex + _forwardNeighbourOffsets[j];
				if (neighbourIndex >= numCells) {
					assert(cell.isHaloCell());
					continue;
				}
				cell.updateLeavingMolecules(_cells[neighbourIndex]);
			}
		}

		#if defined(_OPENMP)
		#pragma omp for schedule(static)
		#endif
		for (vector<ParticleCell>::size_type cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
			_cells[cellIndex].postUpdateLeavingMolecules();
		}
	} // end pragma omp parallel

	/*
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (int i=0; i < numCells; i++) {
		std::vector<Molecule> & molsToSort = _cells[i].filterLeavingMolecules();

		for (auto it = molsToSort.begin(); it != molsToSort.end(); ++it) {
			addParticle(*it);
		}
		molsToSort.clear();
	}
	*/
	_cellsValid = true;
}

bool LinkedCells::addParticle(Molecule& particle, bool inBoxCheckedAlready, bool checkWhetherDuplicate, const bool& rebuildCaches) {
	const bool inBox = inBoxCheckedAlready or particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax);

	bool wasInserted = false;

	if (inBox) {
		int cellIndex = getCellIndexOfMolecule(&particle);
		wasInserted = _cells[cellIndex].addParticle(particle, checkWhetherDuplicate);
		if (rebuildCaches) {
			_cells[cellIndex].buildSoACaches();
		}
	}

	return wasInserted;
}

int LinkedCells::addParticles(vector<Molecule>& particles, bool checkWhetherDuplicate) {
	int oldNumberOfParticles = getNumberOfParticles();

	typedef vector<Molecule>::size_type index_t;
	static vector< vector<ParticleCell>::size_type > index_vector;
	index_vector.resize(particles.size());

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int thread_id = mardyn_get_thread_num();
		const int num_threads = mardyn_get_num_threads();

		const index_t start = (thread_id    ) * particles.size() / num_threads;
		const index_t end   = (thread_id + 1) * particles.size() / num_threads;
		for (index_t i = start; i < end; i++) {
			Molecule particle = particles[i];
			//there could be particles meant for another process that get here...
			/*if(!particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax)){
				index_vector[i] = -1; //this particle is not meant for this process...
				continue;
			}*/
			#ifndef NDEBUG
				if(!particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax)){
					global_log->error()<<"At particle with ID "<<particle.id()<<" assertion failed..."<<endl;
				}
				assert(particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax));
			#endif

			const unsigned long cellIndex = getCellIndexOfMolecule(&particle);
			assert(cellIndex < _cells.size());
			index_vector[i] = cellIndex;
		}

		#if defined(_OPENMP)
		#pragma omp barrier
		#endif

		typedef vector<ParticleCell>::size_type cell_index_t;
		const cell_index_t cells_start = (thread_id    ) * _cells.size() / num_threads;
		const cell_index_t cells_end   = (thread_id + 1) * _cells.size() / num_threads;

		for (index_t i = 0; i < particles.size(); i++) {
			index_t index = index_vector[i];
			if (cells_start <= index and index < cells_end) {
				_cells[index].addParticle(particles[i], checkWhetherDuplicate);
				//delete particles[i]; //replace at end with particles.clear()?
			}
		}
	} // end pragma omp parallel

	index_vector.clear();

	int numberOfAddedParticles = getNumberOfParticles() - oldNumberOfParticles;
	global_log->debug()<<"In LinkedCells::addParticles :"<<endl;
	global_log->debug()<<"\t#Particles to be added = "<<particles.size()<<endl;
	global_log->debug()<<"\t#Particles actually added = "<<numberOfAddedParticles<<endl;

	return numberOfAddedParticles;
}

/**
 * @todo replace this by a call to component->getNumMolecules() !?
 */
unsigned LinkedCells::countParticles(unsigned int cid) {
	unsigned N = 0;
	for (unsigned i = 0; i < _cells.size(); i++) {
		ParticleCell& currentCell = _cells[i];
		if (!currentCell.isHaloCell()) {
			const int numMols = currentCell.getMoleculeCount();
			for (int i = 0; i < numMols; ++i) {
				Molecule& m = currentCell.moleculesAt(i);
				if(m.componentid() == cid)
					++N;
			}
		}
	}
	return N;
}

// @todo: couldn't this use getRegion?
unsigned LinkedCells::countParticles(unsigned int cid, double* cbottom, double* ctop) {
	int minIndex[3];
	int maxIndex[3];
	for (int d = 0; d < 3; d++) {
		if (cbottom[d] < this->_haloBoundingBoxMin[d])
			minIndex[d] = 0;
		else
			minIndex[d] = (int) floor((cbottom[d] - this->_haloBoundingBoxMin[d]) / _cellLength[d]);

		if (ctop[d] > this->_haloBoundingBoxMax[d])
			maxIndex[d] = (int) floor((this->_haloBoundingBoxMax[d] - _haloBoundingBoxMin[d]) / this->_cellLength[d]);
		else
			maxIndex[d] = (int) floor((ctop[d] - this->_haloBoundingBoxMin[d]) / _cellLength[d]);

		if (minIndex[d] < 0)
			minIndex[d] = 0;
		if (maxIndex[d] >= _cellsPerDimension[d])
			maxIndex[d] = _cellsPerDimension[d] - 1;
	}

	unsigned N = 0;
	int cix[3];
	bool individualCheck;
	int cellid;

	for (cix[0] = minIndex[0]; maxIndex[0] >= cix[0]; (cix[0])++) {
		for (cix[1] = minIndex[1]; maxIndex[1] >= cix[1]; (cix[1])++) {
			for (cix[2] = minIndex[2]; maxIndex[2] >= cix[2]; (cix[2])++) {
				individualCheck = (cix[0] == minIndex[0])
						|| (cix[0] == minIndex[0] + 1)
						|| (cix[0] == maxIndex[0])
						|| (cix[0] == maxIndex[0] - 1)
						|| (cix[1] == minIndex[1])
						|| (cix[1] == minIndex[1] + 1)
						|| (cix[1] == maxIndex[1])
						|| (cix[1] == maxIndex[1] - 1)
						|| (cix[2] == minIndex[2])
						|| (cix[2] == minIndex[2] + 1)
						|| (cix[2] == maxIndex[2])
						|| (cix[2] == maxIndex[2] - 1);
				cellid = this->cellIndexOf3DIndex(cix[0], cix[1], cix[2]);
				ParticleCell& currentCell = _cells[cellid];
				if (currentCell.isHaloCell())
					continue;

				const int numMols = currentCell.getMoleculeCount();
				if (individualCheck) {
					for (int i = 0; i < numMols; ++i) {
						Molecule& m = currentCell.moleculesAt(i);
						if(m.inBox(cbottom, ctop) and m.componentid() == cid)
							++N;
					}
				} else {
					for (int i = 0; i < numMols; ++i) {
						Molecule& m = currentCell.moleculesAt(i);
						if(m.componentid() == cid)
							++N;
					}
				}
			}
		}
	}

	return N;
}

void LinkedCells::traverseNonInnermostCells(CellProcessor& cellProcessor) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (traversePairs) invalid, call update first" << endl;
		global_simulation->exit(1);
	}
	// loop over all inner cells and calculate forces to forward neighbours

	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		ParticleCell& currentCell = _cells[cellIndex];
		if (!currentCell.isInnerMostCell()) {
			traverseCell(cellIndex, cellProcessor);
		}
	} // loop over all cells
}

void LinkedCells::traversePartialInnermostCells(CellProcessor& cellProcessor, unsigned int stage, int stageCount) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (traversePairs) invalid, call update first" << endl;
		global_simulation->exit(1);
	}

	// loop over parts of innermost cells and calculate forces to forward neighbours
	// _innerMostCellIndices
	const long int lower =  _innerMostCellIndices.size() * stage / stageCount;
	const long int upper =  _innerMostCellIndices.size() * (stage+1) / stageCount;

#ifndef NDEBUG
	global_log->debug() << "LinkedCells::traversePartialInnermostCells: Processing cells." << endl;
	global_log->debug() << "lower=" << lower << "; upper="
			<< upper << endl;
#endif

	for (long int cellIndex = lower; cellIndex < upper; cellIndex++) {
		traverseCell(_innerMostCellIndices[cellIndex], cellProcessor);
	} // loop over all cells
}

void LinkedCells::traverseCell(const long int cellIndex, CellProcessor& cellProcessor) {

	ParticleCell& currentCell = _cells[cellIndex];
	if (currentCell.isInnerCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	} else if (currentCell.isHaloCell()) {
		cellProcessor.processCell(currentCell);
		for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			long int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
			if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) ((_cells.size()))))
				continue;

			ParticleCell& neighbourCell = _cells[neighbourCellIndex];
			if (!neighbourCell.isHaloCell())
				continue;

			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	} else
	// loop over all boundary cells and calculate forces to forward and backward neighbours
	if (currentCell.isBoundaryCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
		// loop over all backward neighbours. calculate only forces
		// to neighbour cells in the halo region, all others already have been calculated
		for (auto neighbourOffsetsIter = _backwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			ParticleCell& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];
			if (neighbourCell.isHaloCell()) {
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}
	} // if ( isBoundaryCell() )
}

void LinkedCells::traverseCells(CellProcessor& cellProcessor) {
	if (_cellsValid == false) {
		global_log->error()
				<< "Cell structure in LinkedCells (traversePairs) invalid, call update first"
				<< endl;
		global_simulation->exit(1);
	}

	#if defined(_OPENMP)
		traverseCellsC08(cellProcessor);
	#else
		traverseCellsOrig(cellProcessor);
	#endif
}

void LinkedCells::traverseCellsOrig(CellProcessor& cellProcessor) {
	vector<long int>::iterator neighbourOffsetsIter;

#ifndef NDEBUG
	global_log->debug()
			<< "LinkedCells::traverseCells: Processing pairs."
			<< endl;
	global_log->debug() << "_minNeighbourOffset=" << _minNeighbourOffset
			<< "; _maxNeighbourOffset=" << _maxNeighbourOffset << endl;
#endif

	cellProcessor.initTraversal();

	// loop over all inner cells and calculate forces to forward neighbours
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		traverseCell(cellIndex, cellProcessor);
	} // loop over all cells

	cellProcessor.endTraversal();
}

void LinkedCells::traverseCellsC08(CellProcessor& cellProcessor) {
	cellProcessor.initTraversal();

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int strides[3] = {2, 2, 2};
		for (int col = 0; col < 8; ++col) {
			int startIndices[3];
			threeDIndexOfCellIndex(col, startIndices, strides);

			#if defined (_OPENMP)
			#pragma omp for schedule(dynamic, 1) collapse(3)
			#endif
			for (int z = startIndices[2]; z < _cellsPerDimension[2]-1 ; z+= strides[2]) {
				for (int y = startIndices[1]; y < _cellsPerDimension[1]-1; y += strides[1]) {
					for (int x = startIndices[0]; x < _cellsPerDimension[0]-1; x += strides[0]) {
						long int baseIndex = cellIndexOf3DIndex(x, y, z);

						const int num_pairs = _cellPairOffsets.size();
						for(int j = 0; j < num_pairs; ++j) {
							pair<long int, long int> current_pair = _cellPairOffsets[j];

							long int offset1 = current_pair.first;
							long int cellIndex1 = baseIndex + offset1;
							if ((cellIndex1 < 0) || (cellIndex1 >= (int) (_cells.size())))
								continue;

							long int offset2 = current_pair.second;
							long int cellIndex2 = baseIndex + offset2;
							if ((cellIndex2 < 0) || (cellIndex2 >= (int) (_cells.size())))
								continue;

							ParticleCell& cell1 = _cells[cellIndex1];
							ParticleCell& cell2 = _cells[cellIndex2];

							if(cell1.isHaloCell() and cell2.isHaloCell()) {
								continue;
							}

							if(cellIndex1 == cellIndex2) {
								cellProcessor.processCell(cell1);
							}
							else {
								if(!cell1.isHaloCell()) {
									cellProcessor.processCellPair(cell1, cell2);
								}
								else {
									cellProcessor.processCellPair(cell2, cell1);
								}
							}
						}
					}
				}
			}
		}
	} // end pragma omp parallel

	cellProcessor.endTraversal();
}

unsigned long LinkedCells::getNumberOfParticles() {
	unsigned long N = 0;
	std::vector<ParticleCell>::iterator it;
	for (it = _cells.begin(); it != _cells.end(); ++it) {
		N += it->getMoleculeCount();
	}
	return N;
}

void LinkedCells::clear() {
	vector<ParticleCell>::iterator cellIter;
	for (cellIter = _cells.begin(); cellIter != _cells.end(); cellIter++) {
		cellIter->deallocateAllParticles();
	}
}

void LinkedCells::deleteParticlesOutsideBox(double boxMin[3], double boxMax[3]) {
	// This should be unimportant

	ParticleIterator it;
	for (it = iteratorBegin(); it != iteratorEnd();) {
		bool keepMolecule = it->inBox(boxMin, boxMax);
		if (keepMolecule) {
			++it;
		} else {
			it.deleteCurrentParticle();
		}
	}
}

void LinkedCells::deleteOuterParticles() {
	if (_cellsValid == false) {
		global_log->error()
				<< "Cell structure in LinkedCells (deleteOuterParticles) invalid, call update first"
				<< endl;
		global_simulation->exit(1);
	}

	const int numHaloCells = _haloCellIndices.size();

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (int i = 0; i < numHaloCells; i++) {
		ParticleCell& currentCell = _cells[_haloCellIndices[i]];
		currentCell.deallocateAllParticles();
	}
}

double LinkedCells::get_halo_L(int index) const {
	return _haloLength[index];
}

void LinkedCells::getHaloParticles(list<Molecule*> &haloParticlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (getHaloParticles) invalid, call update first" << endl;
		global_simulation->exit(1);
	}

	vector<unsigned long>::iterator cellIndexIter;

	// loop over all halo cells
	for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter != _haloCellIndices.end(); cellIndexIter++) {
		ParticleCell& currentCell = _cells[*cellIndexIter];
		// loop over all molecules in the cell
		const int numMols = currentCell.getMoleculeCount();
		for (int i = 0; i < numMols; ++i) {
			haloParticlePtrs.push_back(&(currentCell.moleculesAt(i)));
		}
	}
}

void LinkedCells::getHaloParticlesDirection(int direction, vector<Molecule>& v, bool removeFromContainer) {
	assert(direction != 0);

	int startIndex[3] = { 0, 0, 0 };
	int stopIndex[3] = { _cellsPerDimension[0] - 1, _cellsPerDimension[1] - 1, _cellsPerDimension[2] - 1 };

	// get dimension in 0, 1, 2 format from direction in +-1, +-2, +-3 format
	unsigned dim = abs(direction) - 1;
	if (direction < 0) {
		stopIndex[dim] = startIndex[dim] + (_haloWidthInNumCells[dim] - 1); // -1 needed for equality below
	}
	else {
		startIndex[dim] = stopIndex[dim] - (_haloWidthInNumCells[dim] - 1); // -1 needed for equality below
	}

	#if 1
		unsigned int startCellIndex = cellIndexOf3DIndex(startIndex[0], startIndex[1], startIndex[2]);
		unsigned int endCellIndex = cellIndexOf3DIndex(stopIndex[0], stopIndex[1], stopIndex[2]);
		vector<vector<Molecule>> threadData;
		vector<int> prefixArray;

		#if defined (_OPENMP)
		#pragma omp parallel shared(v, threadData, startCellIndex, endCellIndex)
		#endif
		{
			const int numThreads = mardyn_get_num_threads();
			const int threadNum = mardyn_get_thread_num();
			RegionParticleIterator begin = iterateRegionBegin(startCellIndex, endCellIndex, IterateType::ALL, removeFromContainer);
			RegionParticleIterator end = iterateRegionEnd();

			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				threadData.resize(numThreads);
				prefixArray.resize(numThreads + 1);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			for(RegionParticleIterator i = begin; i != end; i++){
				//traverse and gather all halo particles in the cells
				//i is a pointer to a Molecule; (*i) is the Molecule
				threadData[threadNum].push_back(*i);
			}

			prefixArray[threadNum + 1] = threadData[threadNum].size();

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//build the prefix array and resize the send buffer
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				int totalNumMols = 0;
				//build the prefix array
				prefixArray[0] = 0;
				for(int i = 1; i <= numThreads; i++){
					prefixArray[i] += prefixArray[i - 1];
					totalNumMols += threadData[i - 1].size();
				}

				//resize the send buffer
				v.resize(totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in v
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				v[prefixArray[threadNum] + i] = threadData[threadNum][i];
			}
		}
	#else
		vector<vector<Molecule>> threadData;
		vector<int> prefixArray;

		#if defined (_OPENMP)
		#pragma omp parallel shared(v, threadData)
		#endif
		{
			const int numThreads = mardyn_get_num_threads();
			const int threadNum = mardyn_get_thread_num();
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				threadData.resize(numThreads);
				prefixArray.resize(numThreads + 1);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			#if defined (_OPENMP)
			#pragma omp for schedule(static) collapse(3)
			#endif
			for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
				for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
					for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
						const int cellIndex = cellIndexOf3DIndex(ix, iy, iz);
						ParticleCell & cell = _cells[cellIndex];

						const int numMols = cell.getMoleculeCount();
						for (int i = 0; i < numMols; ++i) {
							threadData[threadNum].push_back(cell.moleculesAt(i));
						}

						if (removeFromContainer == true) {
							cell.deallocateAllParticles();
						}

						prefixArray[threadNum + 1] += numMols;
					}
				}
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//build the prefix array and allocate memory for v
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				int totalNumMols = 0;
				//build the prefix array
				prefixArray[0] = 0;
				for(int i = 1; i <= numThreads; i++){
					prefixArray[i] += prefixArray[i - 1];
					totalNumMols += threadData[i - 1].size();
				}

				//allocate memory for v
				v.resize(totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in v
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				v[prefixArray[threadNum] + i] = threadData[threadNum].back();
				threadData[threadNum].pop_back();
			}
		}
		/*
		for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
			for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
				for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
					const int cellIndex = cellIndexOf3DIndex(ix, iy, iz);
					ParticleCell & cell = _cells[cellIndex];

					const int numMols = cell.getMoleculeCount();
					for (int i = 0; i < numMols; ++i) {
						v.push_back(cell.moleculesAt(i));
					}

					if (removeFromContainer == true) {
						cell.deallocateAllParticles();
					}
				}
			}
		}
		*/
	#endif
}

void LinkedCells::getBoundaryParticlesDirection(int direction, vector<Molecule>& v) {
	assert(direction != 0);

	int startIndex[3] = { 0, 0, 0 };
	int stopIndex[3] = { _cellsPerDimension[0] - 1, _cellsPerDimension[1] - 1, _cellsPerDimension[2] - 1 };

	// get dimension in 0, 1, 2 format from direction in +-1, +-2, +-3 format
	unsigned dim = abs(direction) - 1;
	if (direction < 0) {
		startIndex[dim] = _haloWidthInNumCells[dim];
		stopIndex[dim] = startIndex[dim] + (_haloWidthInNumCells[dim] - 1); // -1 needed for equality below
	}
	else {
		stopIndex[dim] = _boxWidthInNumCells[dim];
		startIndex[dim] = stopIndex[dim] - (_haloWidthInNumCells[dim] - 1); // -1 needed for equality below
	}

	#if 1
		unsigned int startCellIndex = cellIndexOf3DIndex(startIndex[0], startIndex[1], startIndex[2]);
		unsigned int endCellIndex = cellIndexOf3DIndex(stopIndex[0], stopIndex[1], stopIndex[2]);
		vector<vector<Molecule>> threadData;
		vector<int> prefixArray;

		#if defined (_OPENMP)
		#pragma omp parallel shared(v, threadData, startCellIndex, endCellIndex)
		#endif
		{
			const int numThreads = mardyn_get_num_threads();
			const int threadNum = mardyn_get_thread_num();
			RegionParticleIterator begin = iterateRegionBegin(startCellIndex, endCellIndex);
			RegionParticleIterator end = iterateRegionEnd();

			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				threadData.resize(numThreads);
				prefixArray.resize(numThreads + 1);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			for(RegionParticleIterator i = begin; i != end; i++){
				//traverse and gather all halo particles in the cells
				//i is a pointer to a Molecule; (*i) is the Molecule
				threadData[threadNum].push_back(*i);
			}

			prefixArray[threadNum + 1] = threadData[threadNum].size();

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//build the prefix array and resize the send buffer
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				int totalNumMols = 0;
				//build the prefix array
				prefixArray[0] = 0;
				for(int i = 1; i <= numThreads; i++){
					prefixArray[i] += prefixArray[i - 1];
					totalNumMols += threadData[i - 1].size();
				}

				//resize the send buffer
				v.resize(totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in v
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				v[prefixArray[threadNum] + i] = threadData[threadNum][i];
			}
		}
	#else
		vector<vector<Molecule>> threadData;
		vector<int> prefixArray;

		#if defined (_OPENMP)
		#pragma omp parallel shared(v, threadData)
		#endif
		{
			const int numThreads = mardyn_get_num_threads();
			const int threadNum = mardyn_get_thread_num();
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				threadData.resize(numThreads);
				prefixArray.resize(numThreads + 1);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			#if defined (_OPENMP)
			#pragma omp for schedule(static) collapse(3)
			#endif
			for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
				for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
					for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
						const int cellIndex = cellIndexOf3DIndex(ix, iy, iz);
						ParticleCell & cell = _cells[cellIndex];

						const int numMols = cell.getMoleculeCount();
						for (int i = 0; i < numMols; ++i) {
							threadData[threadNum].push_back(cell.moleculesAt(i));
						}

						prefixArray[threadNum + 1] += numMols;
					}
				}
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//build the prefix array and allocate memory for v
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				int totalNumMols = 0;
				//build the prefix array
				prefixArray[0] = 0;
				for(int i = 1; i <= numThreads; i++){
					prefixArray[i] += prefixArray[i - 1];
					totalNumMols += threadData[i - 1].size();
				}

				//allocate memory for v
				v.resize(totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in v
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				v[prefixArray[threadNum] + i] = threadData[threadNum].back();
				threadData[threadNum].pop_back();
			}
		}
		/*
		for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
			for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
				for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
					const int cellIndex = cellIndexOf3DIndex(ix, iy, iz);
					ParticleCell & cell = _cells[cellIndex];

					const int numMols = cell.getMoleculeCount();
					for (int i = 0; i < numMols; ++i) {
						v.push_back(cell.moleculesAt(i));
					}
				}
			}
		}
		*/
	#endif
}

void LinkedCells::getRegionSimple(double lowCorner[3], double highCorner[3], vector<Molecule*> &particlePtrs, bool removeFromContainer) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (getRegionSimple) invalid, call update first" << endl;
		global_simulation->exit(1);
	}

	int startIndex[3];
	int stopIndex[3];

	for (int dim = 0; dim < 3; dim++) {
		if (lowCorner[dim] <= this->_haloBoundingBoxMax[dim] && highCorner[dim] >= this->_haloBoundingBoxMin[dim]) {
			startIndex[dim] = (int) floor( (lowCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]);
			stopIndex[dim] = (int) floor( (highCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]);
			if (startIndex[dim] < 0)
				startIndex[dim] = 0;
			if (stopIndex[dim] > _cellsPerDimension[dim] - 1)
				stopIndex[dim] = _cellsPerDimension[dim] - 1;
		}
		else {
			// No Part of the given region is owned by this process
			return;
		}
	}

	#if 1
		const int prevNumMols = particlePtrs.size();
		vector<vector<Molecule*>> threadData;
		vector<int> prefixArray;

		#if defined (_OPENMP)
		#pragma omp parallel shared(particlePtrs, threadData)
		#endif
		{
			const int numThreads = mardyn_get_num_threads();
			const int threadNum = mardyn_get_thread_num();
			RegionParticleIterator begin = iterateRegionBegin(lowCorner, highCorner);
			RegionParticleIterator end = iterateRegionEnd();

			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				threadData.resize(numThreads);
				prefixArray.resize(numThreads + 1);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			for(RegionParticleIterator i = begin; i != end; ){
				//traverse and gather all halo particles in the cells
				//i is a pointer to a Molecule; (*i) is the Molecule
				if ((*i).inBox(lowCorner, highCorner)) {
					Molecule* m = &(*i);
					if (not removeFromContainer) {
						threadData[threadNum].push_back(m);
					}
					else {
						threadData[threadNum].push_back(new Molecule(*m));
						i.removeCurrentMoleculeFromContainer();
						// i is already at next molecule, so continue without incrementing
						continue;
					}
				}
				i++;
			}

			prefixArray[threadNum + 1] = threadData[threadNum].size();

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//build the prefix array and resize the send buffer
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				int totalNumMols = 0;
				//build the prefix array
				prefixArray[0] = 0;
				for(int i = 1; i <= numThreads; i++){
					prefixArray[i] += prefixArray[i - 1];
					totalNumMols += threadData[i - 1].size();
				}

				//resize the send buffer
				particlePtrs.resize(prevNumMols + totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in particlePtrs
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				particlePtrs[prevNumMols + prefixArray[threadNum] + i] = threadData[threadNum][i];
			}
		}
	#else
		const int prevNumMols = particlePtrs.size();
		vector<vector<Molecule*>> threadData;
		vector<int> prefixArray;

		#if defined (_OPENMP)
		#pragma omp parallel shared(particlePtrs, threadData)
		#endif
		{
			const int numThreads = mardyn_get_num_threads();
			const int threadNum = mardyn_get_thread_num();
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				threadData.resize(numThreads);
				prefixArray.resize(numThreads + 1);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			#if defined (_OPENMP)
			#pragma omp for schedule(static) collapse(3)
			#endif
			for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
				for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
					for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
						const int cellIndex = cellIndexOf3DIndex(ix, iy, iz);
						ParticleCell & currentCell = _cells[cellIndex];

						// traverse all molecules in the current cell, save to-go-molecules in threadData[threadNum]
						currentCell.getRegion(lowCorner, highCorner, threadData[threadNum], removeFromContainer);
					}
				}
			}

			prefixArray[threadNum + 1] = threadData[threadNum].size();

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//build the prefix array and allocate memory for particlePtrs
			#if defined (_OPENMP)
			#pragma omp master
			#endif
			{
				int totalNumMols = 0;
				//build the prefix array
				prefixArray[0] = 0;
				for(int i = 1; i <= numThreads; i++){
					prefixArray[i] += prefixArray[i - 1];
					totalNumMols += threadData[i - 1].size();
				}

				//allocate / extend memory for particlePtrs
				particlePtrs.resize(prevNumMols + totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in particlePtrs
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				particlePtrs[prevNumMols + prefixArray[threadNum] + i] = threadData[threadNum].back();
				threadData[threadNum].pop_back();
			}
		}
		/*
		int globalCellIndex;
		for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
			for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
				for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
					// globalCellIndex is the cellIndex of the molecule on the coarse Cell level.
					globalCellIndex = cellIndexOf3DIndex(ix, iy, iz);
					// loop over all subcells (either 1 or 8)
					// traverse all molecules in the current cell
					ParticleCell & currentCell = _cells[globalCellIndex];
					currentCell.getRegion(lowCorner, highCorner, particlePtrs, removeFromContainer);
				}
			}
		}
		*/
	#endif
}

void LinkedCells::getRegion(double lowCorner[3], double highCorner[3], vector<Molecule*> &particlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (getRegion) invalid, call update first" << endl;
		global_simulation->exit(1);
	}

	int startIndex[3];
	int stopIndex[3];
	int globalCellIndex;

	for (int dim = 0; dim < 3; dim++) {
		if (lowCorner[dim] < this->_boundingBoxMax[dim] && highCorner[dim] > this->_boundingBoxMin[dim]) {
			startIndex[dim] = (int) floor( (lowCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) - 1;
			stopIndex[dim] = (int) floor( (highCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) + 1;
			if (startIndex[dim] < 0)
				startIndex[dim] = 0;
			if (stopIndex[dim] > _cellsPerDimension[dim] - 1)
				stopIndex[dim] = _cellsPerDimension[dim] - 1;
		} else {
			// No Part of the given region is owned by this process
			return;
		}
	}

	for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
		for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
			for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
				// globalCellIndex is the cellIndex of the molecule on the coarse Cell level.
				globalCellIndex = cellIndexOf3DIndex(ix, iy, iz);
				// loop over all subcells (either 1 or 8)
				// traverse all molecules in the current cell
				ParticleCell & currentCell = _cells[globalCellIndex];
				currentCell.getRegion(lowCorner, highCorner, particlePtrs);
			}
		}
	}
}

int LinkedCells::countNeighbours(ParticlePairsHandler* /*particlePairsHandler*/, Molecule* m1, CellProcessor& cellProcessor, double RR) {
        int m1neigh = 0;
        assert(_cellsValid);
        unsigned long cellIndex = getCellIndexOfMolecule(m1);
        ParticleCell& currentCell = _cells[cellIndex];

        cellProcessor.initTraversal();

        // extend the window of cells with cache activated
        for (unsigned int windowCellIndex = cellIndex - _minNeighbourOffset; windowCellIndex < cellIndex + _maxNeighbourOffset+1 ; windowCellIndex++) {
                cellProcessor.preprocessCell(_cells[windowCellIndex]);
        }

        m1neigh += cellProcessor.countNeighbours(m1, currentCell, RR);

        // forward neighbours
        for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++)
        {
                ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
                m1neigh += cellProcessor.countNeighbours(m1, neighbourCell, RR);
        }
        // backward neighbours
        for (auto neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++)
        {
                ParticleCell& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];  // minus oder plus?
                m1neigh += cellProcessor.countNeighbours(m1, neighbourCell, RR);
        }

        // close the window of cells activated
        for (unsigned int windowCellIndex = cellIndex - _minNeighbourOffset; windowCellIndex < cellIndex + _maxNeighbourOffset+1; windowCellIndex++) {
                cellProcessor.postprocessCell(_cells[windowCellIndex]);
        }

        cellProcessor.endTraversal();
        return m1neigh;
}

unsigned long LinkedCells::numCavities(CavityEnsemble* ce, DomainDecompBase* comm) {
   return ce->communicateNumCavities(comm);
}

void LinkedCells::cavityStep(CavityEnsemble* ce, double /*T*/, Domain* domain, CellProcessor& cellProcessor) {
   ParticlePairs2PotForceAdapter particlePairsHandler(*domain);
   map<unsigned long, Molecule*>* pc = ce->particleContainer();
   double RR = ce->getRR();
   
   for(map<unsigned long, Molecule*>::iterator pcit = pc->begin(); pcit != pc->end(); pcit++) {
      assert(pcit->second != NULL);
      Molecule* m1 = pcit->second;
      unsigned neigh = this->countNeighbours(&particlePairsHandler, m1, cellProcessor, RR);
      unsigned long m1id = pcit->first;
      assert(m1id == m1->id());
      ce->decideActivity(neigh, m1id);
   }
}

RegionParticleIterator LinkedCells::iterateRegionBegin(const unsigned int startCellIndex, const unsigned int endCellIndex, IterateType type, bool removeFromContainer) {
	// parameters "type" and "removeFromContainer" not yet used
	// add functionality in a future version...

	int start3DIndices[3], end3DIndices[3];
	int regionDimensions[3];
	threeDIndexOfCellIndex(startCellIndex, start3DIndices, _cellsPerDimension);
	threeDIndexOfCellIndex(endCellIndex, end3DIndices, _cellsPerDimension);
	for(int d = 0; d < 3; d++){
		regionDimensions[d] = end3DIndices[d] - start3DIndices[d] + 1;
	}
	ParticleIterator::CellIndex_T offset = mardyn_get_thread_num(); // starting position
	ParticleIterator::CellIndex_T stride = mardyn_get_num_threads(); // stride

	return RegionParticleIterator(&_cells, offset, stride, startCellIndex, regionDimensions, _cellsPerDimension);
}

RegionParticleIterator LinkedCells::iterateRegionBegin(const double startCorner[3], const double endCorner[3], IterateType type, bool removeFromContainer) {
	unsigned int startRegionCellIndex;
	unsigned int endRegionCellIndex;

	getCellIndicesOfRegion(startCorner, endCorner, startRegionCellIndex, endRegionCellIndex);

	return iterateRegionBegin(startRegionCellIndex, endRegionCellIndex, type, removeFromContainer);
}

RegionParticleIterator LinkedCells::iterateRegionEnd() {
	return RegionParticleIterator::invalid();
}

//################################################
//############ PRIVATE METHODS ###################
//################################################

void LinkedCells::initializeCells() {
	_innerMostCellIndices.clear();
	_innerCellIndices.clear();
	_boundaryCellIndices.clear();
	_haloCellIndices.clear();

	long int cellIndex;
	double cellBoxMin[3], cellBoxMax[3];

	for (int iz = 0; iz < _cellsPerDimension[2]; ++iz) {
		cellBoxMin[2] = iz * _cellLength[2] + _haloBoundingBoxMin[2];
		cellBoxMax[2] = (iz + 1) * _cellLength[2] + _haloBoundingBoxMin[2];

		for (int iy = 0; iy < _cellsPerDimension[1]; ++iy) {
			cellBoxMin[1] = iy * _cellLength[1] + _haloBoundingBoxMin[1];
			cellBoxMax[1] = (iy + 1) * _cellLength[1] + _haloBoundingBoxMin[1];

			for (int ix = 0; ix < _cellsPerDimension[0]; ++ix) {
				cellBoxMin[0] = ix * _cellLength[0] + _haloBoundingBoxMin[0];
				cellBoxMax[0] = (ix + 1) * _cellLength[0] + _haloBoundingBoxMin[0];

				cellIndex = cellIndexOf3DIndex(ix, iy, iz);
				ParticleCell & cell = _cells[cellIndex];

				cell.skipCellFromHaloRegion();
				cell.skipCellFromBoundaryRegion();
				cell.skipCellFromInnerRegion();
				cell.skipCellFromInnerMostRegion();

				cell.setBoxMin(cellBoxMin);
				cell.setBoxMax(cellBoxMax);
				_cells[cellIndex].setCellIndex(cellIndex); //set the index of the cell to the index of it...
				if (ix < _haloWidthInNumCells[0] ||
					iy < _haloWidthInNumCells[1] ||
					iz < _haloWidthInNumCells[2] ||
					ix >= _cellsPerDimension[0] - _haloWidthInNumCells[0] ||
					iy >= _cellsPerDimension[1] - _haloWidthInNumCells[1] ||
					iz >= _cellsPerDimension[2] - _haloWidthInNumCells[2]) {

					cell.assignCellToHaloRegion();
					_haloCellIndices.push_back(cellIndex);
				}
				else{
					if (ix < 2 * _haloWidthInNumCells[0] ||
						iy < 2 * _haloWidthInNumCells[1] ||
						iz < 2 * _haloWidthInNumCells[2] ||
						ix >= _cellsPerDimension[0] - 2 * _haloWidthInNumCells[0] ||
						iy >= _cellsPerDimension[1] - 2 * _haloWidthInNumCells[1] ||
						iz >= _cellsPerDimension[2] - 2 * _haloWidthInNumCells[2]) {

						cell.assignCellToBoundaryRegion();
						_boundaryCellIndices.push_back(cellIndex);
					}
					else{
						if (ix < 3 * _haloWidthInNumCells[0] ||
							iy < 3 * _haloWidthInNumCells[1] ||
							iz < 3 * _haloWidthInNumCells[2] ||
							ix >= _cellsPerDimension[0] - 3 * _haloWidthInNumCells[0] ||
							iy >= _cellsPerDimension[1] - 3 * _haloWidthInNumCells[1] ||
							iz >= _cellsPerDimension[2] - 3 * _haloWidthInNumCells[2]) {

							cell.assignCellToInnerRegion();
							_innerCellIndices.push_back(cellIndex);
						}
						else {
							cell.assignCellToInnerMostAndInnerRegion();
							_innerMostCellIndices.push_back(cellIndex);
							_innerCellIndices.push_back(cellIndex);
						}
					}
				}
			}
		}
	}

	/*********************  Compute Border indices *********************/
	for (int dim = 0; dim < 3; ++dim) {
		for (int dir = 0; dir < 2; ++dir) {
			for (int typ = 0; typ < 2; ++typ) {
				_borderCellIndices[dim][dir][typ].clear();
				int low[3] = { 0, 0, 0 };
				int high[3] = { _cellsPerDimension[0] - 1, _cellsPerDimension[1] - 1, _cellsPerDimension[2] - 1 };

				if (typ == 1) {
					if (dir == 0){
						low[dim]++;
					}
					else {
						high[dim]--;
					}
				}

				if (dir == 0){
					high[dim] = low[dim];
				}
				else {
					low[dim] = high[dim];
				}

				for (int iz = low[2]; iz <= high[2]; ++iz) {
					for (int iy = low[1]; iy <= high[1]; ++iy) {
						for (int ix = low[0]; ix <= high[0]; ++ix) {
							cellIndex = cellIndexOf3DIndex(ix, iy, iz);
#ifndef NDEBUG
							ParticleCell & cell = _cells[cellIndex];

							assert(not cell.isInnerCell());

							if (typ == 0){
								assert(cell.isHaloCell());
							} else {
								 /* assert(cell.isBoundaryCell()) is not always true, as we have some halo cells in there */
							}
#endif

							_borderCellIndices[dim][dir][typ].push_back(cellIndex);
						}
					}
				}
			}
		}
	}
}

void LinkedCells::calculateNeighbourIndices() {
	global_log->debug() << "Setting up cell neighbour indice lists." << endl;
	std::fill(_forwardNeighbourOffsets.begin(), _forwardNeighbourOffsets.end(), 0);
	std::fill(_backwardNeighbourOffsets.begin(), _backwardNeighbourOffsets.end(), 0);
	int forwardNeighbourIndex = 0, backwardNeighbourIndex = 0;

	_maxNeighbourOffset = 0;
	_minNeighbourOffset = 0;
	double xDistanceSquare;
	double yDistanceSquare;
	double zDistanceSquare;
	double cutoffRadiusSquare = pow(_cutoffRadius, 2);
	for (int zIndex = -_haloWidthInNumCells[2];
			zIndex <= _haloWidthInNumCells[2]; zIndex++) {
		// The distance in one dimension is the width of a cell multiplied with the number
		// of cells between the two cells (this is received by substracting one of the
		// absolute difference of the cells, if this difference is not zero)
		if (zIndex != 0) {
			zDistanceSquare = pow((abs(zIndex) - 1) * _cellLength[2], 2);
		} else {
			zDistanceSquare = 0;
		}
		for (int yIndex = -_haloWidthInNumCells[1];
				yIndex <= _haloWidthInNumCells[1]; yIndex++) {
			if (yIndex != 0) {
				yDistanceSquare = pow((abs(yIndex) - 1) * _cellLength[1], 2);
			} else {
				yDistanceSquare = 0;
			}
			for (int xIndex = -_haloWidthInNumCells[0];
					xIndex <= _haloWidthInNumCells[0]; xIndex++) {
				if (xIndex != 0) {
					xDistanceSquare = pow((abs(xIndex) - 1) * _cellLength[0],
							2);
				} else {
					xDistanceSquare = 0;
				}
				if (xDistanceSquare + yDistanceSquare + zDistanceSquare
						<= cutoffRadiusSquare) {
					long int offset = cellIndexOf3DIndex(xIndex, yIndex,
							zIndex);
					if (offset > 0) {
						_forwardNeighbourOffsets[forwardNeighbourIndex] = offset;
						++forwardNeighbourIndex;
						if (offset > _maxNeighbourOffset) {
							_maxNeighbourOffset = offset;
						}
					}
					if (offset < 0) {
						_backwardNeighbourOffsets[backwardNeighbourIndex] = abs(offset);
						++backwardNeighbourIndex;
						if (abs(offset) > _minNeighbourOffset) {
							_minNeighbourOffset = abs(offset);
						}
					}
				}
			}
		}
	}

	assert(forwardNeighbourIndex == 13);
	assert(backwardNeighbourIndex == 13);

	global_log->info() << "Neighbour offsets are bounded by "
			<< _minNeighbourOffset << ", " << _maxNeighbourOffset << endl;
}

void LinkedCells::calculateCellPairOffsets() {
	long int o   = cellIndexOf3DIndex(0,0,0); // origin
	long int x   = cellIndexOf3DIndex(1,0,0); // displacement to the right
	long int y   = cellIndexOf3DIndex(0,1,0); // displacement ...
	long int z   = cellIndexOf3DIndex(0,0,1);
	long int xy  = cellIndexOf3DIndex(1,1,0);
	long int yz  = cellIndexOf3DIndex(0,1,1);
	long int xz  = cellIndexOf3DIndex(1,0,1);
	long int xyz = cellIndexOf3DIndex(1,1,1);

	// minimize number of cells simultaneously in memory:

	_cellPairOffsets[ 0] = make_pair(o, xyz);
	// evict xyz

	_cellPairOffsets[ 1] = make_pair(o, yz );
	_cellPairOffsets[ 2] = make_pair(x, yz );
	// evict yz

	_cellPairOffsets[ 3] = make_pair(o, x  );

	_cellPairOffsets[ 4] = make_pair(o, xy );
	_cellPairOffsets[ 5] = make_pair(xy, z );
	// evict xy

	_cellPairOffsets[ 6] = make_pair(o, z  );
	_cellPairOffsets[ 7] = make_pair(x, z  );
	_cellPairOffsets[ 8] = make_pair(y, z  );
	// evict z

	_cellPairOffsets[ 9] = make_pair(o, y  );
	_cellPairOffsets[10] = make_pair(x, y  );
	// evict x

	_cellPairOffsets[11] = make_pair(o, xz );
	_cellPairOffsets[12] = make_pair(y, xz );
	// evict xz

	_cellPairOffsets[13] = make_pair(o, o  );
}

unsigned long int LinkedCells::getCellIndexOfMolecule(Molecule* molecule) const {
	int cellIndex[3]; // 3D Cell index

	for (int dim = 0; dim < 3; dim++) {
		#ifndef NDEBUG
		if (molecule->r(dim) < _haloBoundingBoxMin[dim] || molecule->r(dim) >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "Molecule is outside of bounding box" << endl;
			global_log->error() << "Molecule:\n" << *molecule << endl;
			global_log->error() << "_haloBoundingBoxMin = (" << _haloBoundingBoxMin[0] << ", " << _haloBoundingBoxMin[1] << ", " << _haloBoundingBoxMin[2] << ")" << endl;
			global_log->error() << "_haloBoundingBoxMax = (" << _haloBoundingBoxMax[0] << ", " << _haloBoundingBoxMax[1] << ", " << _haloBoundingBoxMax[2] << ")" << endl;
			global_simulation->exit(1);
		}
		#endif
		//this version is sensitive to roundoffs, if we have molecules (initialized) precisely at position 0.0:
		//cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);
		cellIndex[dim] = ((int) floor((molecule->r(dim) - _boundingBoxMin[dim]) / _cellLength[dim])) + _haloWidthInNumCells[dim];

	}
	int cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
	// in very rare cases rounding is stupid, thus we need a check...
	//TODO: check if this can in any way be done better...
	if (_cells[cellIndex1d].testInBox(*molecule)) {
		return cellIndex1d;
	} else {
		for (int dim = 0; dim < 3; dim++) {
			if(molecule->r(dim)<_cells[cellIndex1d].getBoxMin(dim)){
				cellIndex[dim]--;
			}else if(molecule->r(dim)>=_cells[cellIndex1d].getBoxMax(dim)){
				cellIndex[dim]++;
			}
		}
		cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
		assert(_cells[cellIndex1d].testInBox(*molecule));
		return cellIndex1d;
	}
}

unsigned long int LinkedCells::getCellIndexOfPoint(const double point[3]) const {
	int cellIndex[3]; // 3D Cell index
	double localPoint[] = {point[0], point[1], point[2]};

	for (int dim = 0; dim < 3; dim++) {
		// different than getCellIndexOfMolecule!!!

		// ignore a bit of rounding, if the point is outside of the box.
		if (localPoint[dim] < _haloBoundingBoxMin[dim]){
			localPoint[dim] += _cellLength[dim]/2;
		} 
		else if(localPoint[dim] >= _haloBoundingBoxMax[dim]){
			localPoint[dim] -= _cellLength[dim]/2;
		}

		#ifndef NDEBUG
		//this should never ever happen!
		if (localPoint[dim] < _haloBoundingBoxMin[dim] || localPoint[dim] >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "Point is outside of halo bounding box" << endl;
			global_log->error() << "Point p = (" << localPoint[0] << ", " << localPoint[1] << ", " << localPoint[2] << ")" << endl;
			global_log->error() << "_haloBoundingBoxMin = (" << _haloBoundingBoxMin[0] << ", " << _haloBoundingBoxMin[1] << ", " << _haloBoundingBoxMin[2] << ")" << endl;
			global_log->error() << "_haloBoundingBoxMax = (" << _haloBoundingBoxMax[0] << ", " << _haloBoundingBoxMax[1] << ", " << _haloBoundingBoxMax[2] << ")" << endl;
			global_simulation->exit(1);
		}
		#endif
		
		//this version is sensitive to roundoffs, if we have molecules (initialized) precisely at position 0.0:
		//cellIndex[dim] = (int) floor(point[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]);
		cellIndex[dim] = ((int) floor((localPoint[dim] - _boundingBoxMin[dim]) / _cellLength[dim])) + _haloWidthInNumCells[dim];
	}
	
	int cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
	// in very rare cases rounding is stupid, thus we need a check...
	//TODO: check if this can in any way be done better...
	if (_cells[cellIndex1d].testPointInCell(localPoint)) {
		return cellIndex1d;
	}
	else {
		for (int dim = 0; dim < 3; dim++) {
			if(localPoint[dim]<_cells[cellIndex1d].getBoxMin(dim)){
				cellIndex[dim]--;
			}
			else{
				if(localPoint[dim]>=_cells[cellIndex1d].getBoxMax(dim)){
					cellIndex[dim]++;
				}
			}
		}
		cellIndex1d = this->cellIndexOf3DIndex(cellIndex[0], cellIndex[1], cellIndex[2]);
		assert(_cells[cellIndex1d].testPointInCell(localPoint));
		return cellIndex1d;
	}
}

long int LinkedCells::cellIndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const {
	return (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
}

void LinkedCells::threeDIndexOfCellIndex(int ind, int r[3], const int dim[3]) const {
	r[2] = ind / (dim[0] * dim[1]);
	r[1] = (ind - r[2] * dim[0] * dim[1]) / dim[0];
	r[0] = ind - dim[0] * (r[1] + dim[1] * r[2]);
}

void LinkedCells::getCellIndicesOfRegion(const double startRegion[3], const double endRegion[3], unsigned int &startIndex, unsigned int &endIndex) {
	//find the cellIndices of the cells which contain the start and end corners of the region
	startIndex = getCellIndexOfPoint(startRegion);
	endIndex = getCellIndexOfPoint(endRegion);
}

void LinkedCells::deleteMolecule(unsigned long molid, double x, double y, double z, const bool& rebuildCaches) {
	int ix = (int) floor( (x - this->_haloBoundingBoxMin[0]) / this->_cellLength[0]);
	int iy = (int) floor( (y - this->_haloBoundingBoxMin[1]) / this->_cellLength[1]);
	int iz = (int) floor( (z - this->_haloBoundingBoxMin[2]) / this->_cellLength[2]);

	unsigned long hash = this->cellIndexOf3DIndex(ix, iy, iz);
	if (hash >= _cells.size()) {
		global_log->error_always_output()
				<< "coordinates for atom deletion lie outside bounding box."
				<< endl;
		global_simulation->exit(1);
	}

	bool found = this->_cells[hash].deleteMoleculeByID(molid);

	if (!found) {
		global_log->error_always_output() << "could not delete molecule " << molid << "."
				<< endl;
		global_simulation->exit(1);
	}
	else if (rebuildCaches) {
		_cells[hash].buildSoACaches();
	}
}

double LinkedCells::getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessorI) {
	CellProcessor* cellProcessor;
	if (dynamic_cast<LegacyCellProcessor*>(&cellProcessorI)) {
		cellProcessor = &cellProcessorI;
	} else {
		cellProcessor = new LegacyCellProcessor(cellProcessorI.getCutoffRadius(), cellProcessorI.getLJCutoffRadius(),
				particlePairsHandler);
	}

	double u = 0.0;

	static ParticleCell dummyCell;
	{
		// (potentially re-) initialize dummyCell
		dummyCell.assignCellToInnerRegion();
		double l[3] = {m1->r(0)-10., m1->r(1)-10., m1->r(2)-10.}, u[3] = {m1->r(0)+10., m1->r(1)+10., m1->r(2)+10.};
		dummyCell.setBoxMin(l);
		dummyCell.setBoxMax(u);
		dummyCell.addParticle(*m1);

		dummyCell.buildSoACaches();
		dummyCell.setCellIndex(_cells.back().getCellIndex() * 10);
	}

	unsigned long cellIndex = getCellIndexOfMolecule(m1);

	ParticleCell& currentCell = _cells[cellIndex];

	assert(not currentCell.isHaloCell());

	cellProcessor->initTraversal();

	Molecule * molWithSoA = &(dummyCell.moleculesAt(0));

	u += cellProcessor->processSingleMolecule(molWithSoA, currentCell);

	// forward neighbours
	for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
			neighbourOffsetsIter != _forwardNeighbourOffsets.end();
			neighbourOffsetsIter++) {
		ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
		u += cellProcessor->processSingleMolecule(molWithSoA, neighbourCell);
	}
	// backward neighbours
	for (auto neighbourOffsetsIter = _backwardNeighbourOffsets.begin();
			neighbourOffsetsIter != _backwardNeighbourOffsets.end();
			neighbourOffsetsIter++) {
		ParticleCell& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];
		u += cellProcessor->processSingleMolecule(molWithSoA, neighbourCell);
	}

	cellProcessor->endTraversal();

	if (!dynamic_cast<LegacyCellProcessor*>(&cellProcessorI)) {
		delete cellProcessor;
	}

	dummyCell.deallocateAllParticles();

    assert(not std::isnan(u)); // catches NaN

	return u;
}

void LinkedCells::updateInnerMoleculeCaches() {
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		if(_cells[cellIndex].isInnerCell()){
			_cells[cellIndex].buildSoACaches();
		}
	}
}

void LinkedCells::updateBoundaryAndHaloMoleculeCaches() {
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		if (_cells[cellIndex].isHaloCell() or _cells[cellIndex].isBoundaryCell()) {
			_cells[cellIndex].buildSoACaches();
		}
	}
}

void LinkedCells::updateMoleculeCaches() {
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++) {
		_cells[cellIndex].buildSoACaches();
	}
}
