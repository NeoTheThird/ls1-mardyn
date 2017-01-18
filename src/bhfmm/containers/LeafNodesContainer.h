/*
 * LeafNodesContainer.h
 *
 *  Created on: Feb 6, 2015
 *      Author: tchipev
 *
 *  A bare-bones secondary Linked-Cells structure to allow
 *  carrying out P2P, P2M, L2P operations when the tree is deeper
 *  than what would be allowed by the Lennard-Jones cutoff radius.
 *
 *  Assumption: there is an LJ Linked-Cells structure with cells bigger or equal than these.
 */

#ifndef LEAFNODESCONTAINER_H_
#define LEAFNODESCONTAINER_H_

#include "bhfmm/containers/ParticleCellPointers.h"
#include <vector>

namespace bhfmm {

class SimpleCellProcessor;
class VectorizedChargeP2PCellProcessor;

class LeafNodesContainer {
public:
	LeafNodesContainer(double bBoxMin[3], double bBoxMax[3],
		double LJCellLength[3], unsigned subdivisionFactor,
		bool periodic = true);

	~LeafNodesContainer();

	double* getCellLength() {
		return _cellLength;
	}

	void addParticle(Molecule& particle);
	void clearParticles();
	void traverseCells(SimpleCellProcessor& cellProcessor);
	void traverseCellPairs(VectorizedChargeP2PCellProcessor& cellProcessor);
	void traverseCellPairsOrig(VectorizedChargeP2PCellProcessor& cellProcessor);
	void traverseCellPairsC08(VectorizedChargeP2PCellProcessor& cellProcessor);

private:
	void initializeCells();
	long int cellIndexOf3DIndex(int xIndex, int yIndex, int zIndex) const;
	void calculateNeighbourIndices();
	void calculateCellPairOffsets();
	unsigned long int getCellIndexOfMolecule(Molecule* molecule) const;
	//! @brief given the index in the cell vector, return the index of the cell.
	void threeDIndexOfCellIndex(int ind, int r[3], int dim[3]) const;

	bool _periodicBC;

	double _boundingBoxMin[3];
	double _boundingBoxMax[3];
	double _haloBoundingBoxMin[3];
	double _haloBoundingBoxMax[3];
	double _cellLength[3];
	int _numInnerCellsPerDimension[3];
	int _numCellsPerDimension[3];
	std::vector<ParticleCellPointers> _cells;
	std::vector<unsigned long> _forwardNeighbourOffsets;
	std::vector<unsigned long> _backwardNeighbourOffsets;
	unsigned _maxNeighbourOffset;
	unsigned _minNeighbourOffset;

	// addition for compact SimpleMD-style traversal
	std::vector<std::pair<unsigned long, unsigned long> > _cellPairOffsets;

	unsigned int _numActiveColours;
	std::vector<std::vector<long int> > _cellIndicesPerColour;

};

} /* namespace bhfmm */

#endif /* LEAFNODESCONTAINER_H_ */
