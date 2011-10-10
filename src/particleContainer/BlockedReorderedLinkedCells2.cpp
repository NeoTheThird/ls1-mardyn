/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "particleContainer/BlockedReorderedLinkedCells2.h"

#include <cmath>

#include "molecules/potforce.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "BlockedCell.h"
#include "molecules/MoleculeTypes.h"
#include "parallel/DomainDecompBase.h"
#include "ensemble/GrandCanonical.h"
#include "Domain.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


BlockedReorderedLinkedCells2::BlockedReorderedLinkedCells2(
		double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double LJCutoffRadius,
		double tersoffCutoffRadius, double cellsInCutoffRadius
)
		: ParticleContainer(bBoxMin, bBoxMax),
			_blockTraverse(this, _cells, _innerCellIndices, _boundaryCellIndices, _haloCellIndices )
{
	int numberOfCells = 1;
	_cutoffRadius = cutoffRadius;
	_LJCutoffRadius = LJCutoffRadius;
	_tersoffCutoffRadius = tersoffCutoffRadius;

	for (int d = 0; d < 3; d++) {
		/* first calculate the cell length for this dimension */
		_boxWidthInNumCells[d] = floor((_boundingBoxMax[d] - _boundingBoxMin[d]) / cutoffRadius * cellsInCutoffRadius);
		// in each dimension at least one layer of (inner+boundary) cells is necessary
		if( _boxWidthInNumCells[d] == 0 ) {
			_boxWidthInNumCells[d] = 1;
		}
		_cellLength[d] = (_boundingBoxMax[d] - _boundingBoxMin[d]) / _boxWidthInNumCells[d];

		_haloWidthInNumCells[d] = ceil(cellsInCutoffRadius);
		_haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloLength[d];

		_cellsPerDimension[d] = _boxWidthInNumCells[d] + 2 * _haloWidthInNumCells[d];

		numberOfCells *= _cellsPerDimension[d];
	}
	global_log->debug() << "Cell size (" << _cellLength[0] << ", " << _cellLength[1] << ", " << _cellLength[2] << ")" << endl;

	_cells.resize(numberOfCells);

	// If the width of the inner region is less than the width of the halo region
	// a parallelisation isn't possible (with the used algorithms).
	// In this case, print an error message
	// _cellsPerDimension is 2 times the halo width + the inner width
	// so it has to be at least 3 times the halo width
	if (_boxWidthInNumCells[0] < _haloWidthInNumCells[0] ||
	    _boxWidthInNumCells[1] < _haloWidthInNumCells[1] ||
	    _boxWidthInNumCells[2] < _haloWidthInNumCells[2]) {
		global_log->error() << "BlockedReorderedLinkedCells (constructor): bounding box too small for calculated cell length" << endl;
		global_log->error() << "_cellsPerDimension" << _cellsPerDimension[0] << " / " << _cellsPerDimension[1] << " / " << _cellsPerDimension[2] << endl;
		global_log->error() << "_haloWidthInNumCells" << _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1] << " / " << _haloWidthInNumCells[2] << endl;
		exit(5);
	}
	this->_localInsertionsMinusDeletions = 0;

	initializeCells();
	calculateNeighbourIndices();
	_cellsValid = false;
}


BlockedReorderedLinkedCells2::~BlockedReorderedLinkedCells2() {
}

void BlockedReorderedLinkedCells2::rebuild(double bBoxMin[3], double bBoxMax[3]) {
	for (int i = 0; i < 3; i++) {
		this->_boundingBoxMin[i] = bBoxMin[i];
		this->_boundingBoxMax[i] = bBoxMax[i];
	}

	int numberOfCells = 1;

	for (int dim = 0; dim < 3; dim++) {
		_cellsPerDimension[dim] = (int) floor((this->_boundingBoxMax[dim] - this->_boundingBoxMin[dim]) / (_cutoffRadius / _haloWidthInNumCells[dim]))
		    + 2 * _haloWidthInNumCells[dim];
		// in each dimension at least one layer of (inner+boundary) cells necessary
		if (_cellsPerDimension[dim] == 2 * _haloWidthInNumCells[dim]) {
			global_log->error() << "BlockedReorderedLinkedCells::rebuild: region to small" << endl;
			exit(1);
		}
		numberOfCells *= _cellsPerDimension[dim];
		_cellLength[dim] = (this->_boundingBoxMax[dim] - this->_boundingBoxMin[dim]) / (_cellsPerDimension[dim] - 2 * _haloWidthInNumCells[dim]);
		_haloBoundingBoxMin[dim] = this->_boundingBoxMin[dim] - _haloWidthInNumCells[dim] * _cellLength[dim];
		_haloBoundingBoxMax[dim] = this->_boundingBoxMax[dim] + _haloWidthInNumCells[dim] * _cellLength[dim];
		_haloLength[dim] = _haloWidthInNumCells[dim] * _cellLength[dim];
	}

	_cells.resize(numberOfCells);

	// If the with of the inner region is less than the width of the halo region
	// a parallelisation isn't possible (with the used algorithms).
	// In this case, print an error message
	// _cellsPerDimension is 2 times the halo width + the inner width
	// so it has to be at least 3 times the halo width
	if (_cellsPerDimension[0] < 3*_haloWidthInNumCells[0] ||
	    _cellsPerDimension[1] < 3*_haloWidthInNumCells[1] ||
	    _cellsPerDimension[2] < 3*_haloWidthInNumCells[2]) {
		global_log->error() << "BlockedReorderedLinkedCells (rebuild): bounding box too small for calculated cell Length" << endl;
		global_log->error() << "cellsPerDimension" << _cellsPerDimension[0] << " / " << _cellsPerDimension[1] << " / " << _cellsPerDimension[2] << endl;
		global_log->error() << "_haloWidthInNumCells" << _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1] << " / " << _haloWidthInNumCells[2] << endl;
		exit(5);
	}

	initializeCells();
	calculateNeighbourIndices();

	// delete all Particles which are outside of the halo region
	bool erase_mol;
	for (unsigned int i = 0; i < _cells.size(); i++) {
		MoleculeArray::iterator molIter =
				_cells[i].getParticles().begin();
		while (molIter != _cells[i].getParticles().end()) {
			erase_mol = false;
			for (unsigned short d = 0; d < 3; ++d) {
				const double& rd = molIter->r(d);
				// The molecules has to be within the domain of the process
				// If it is outside in at least one dimension, it has to be
				// erased /
				if (rd < this->_haloBoundingBoxMin[d] || rd >= this->_haloBoundingBoxMax[d])
					erase_mol = true;
			}
			if (erase_mol) {
				molIter = _cells[i].deleteMolecule(molIter);
			} else {
				++molIter;
			}
		}
	}

	_cellsValid = false;
}

void BlockedReorderedLinkedCells2::update() {

	unsigned long index;
	int counter = 0;

	for (unsigned i = 0; i < _cells.size(); i++) {
			BlockedCell& currentCell = _cells[i];
			MoleculeArray& molecules = currentCell.getParticles();
			MoleculeArray::iterator moleculeIterator = molecules.begin();

			while (moleculeIterator != molecules.end()) {
				index = getCellIndexOfMolecule(moleculeIterator.operator ->());
				if (index != i) {
					counter++;
					global_log->debug() << "moving particle " << moleculeIterator->id() << " from cell " << i << " to " << index << endl;
					_cells[index].addParticle(*moleculeIterator);
					moleculeIterator = molecules.erase(moleculeIterator);
				} else {
					++moleculeIterator;
				}
			}
	}

	global_log->debug() << "BlockedReorderedLinkedCell update moved " << counter << " particle(s)." << endl;
	_cellsValid = true;
}

void BlockedReorderedLinkedCells2::addParticle(Molecule& particle) {

	double x = particle.r(0);
	double y = particle.r(1);
	double z = particle.r(2);

	if ( ( x >= _haloBoundingBoxMin[0]) && (x < _haloBoundingBoxMax[0]) &&
			( y >= _haloBoundingBoxMin[1]) && (y < _haloBoundingBoxMax[1]) &&
			( z >= _haloBoundingBoxMin[2]) && (z < _haloBoundingBoxMax[2]) ) {

		int cellIndex = getCellIndexOfMolecule(&particle);
		if (cellIndex < (int) 0 || cellIndex >= (int) _cells.size()) {
			global_log->error() << "BlockedReorderedLinkedCells::addParticle(): INDEX ERROR " << endl;
			exit(1);
		}
		_cells[cellIndex].addParticle(particle);
	}
}

void BlockedReorderedLinkedCells2::traversePairs(ParticlePairsHandler* particlePairsHandler) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in BlockedReorderedLinkedCells (traversePairs) invalid, call update first" << endl;
		exit(1);
	}
	_blockTraverse.traversePairs(particlePairsHandler);
}

unsigned long BlockedReorderedLinkedCells2::getNumberOfParticles() {
	unsigned long size = 0;
	for (unsigned int i = 0; i < _cells.size(); i++) {
		size += _cells[i].getMoleculeCount();
	}

	return size;
}

Molecule* BlockedReorderedLinkedCells2::begin() {
	_cellIteratorIndex = 0;
	_moleculeIteratorIndex = 0;

	while (_cellIteratorIndex < _cells.size()) {
		if (_cells[_cellIteratorIndex].getMoleculeCount() > 0) {
			return &(_cells[_cellIteratorIndex].getParticles()[0]);
		}
		_cellIteratorIndex++;
	}
	return NULL;
}

Molecule* BlockedReorderedLinkedCells2::next() {
	++_moleculeIteratorIndex;
	if (_moleculeIteratorIndex < (unsigned int)_cells[_cellIteratorIndex].getMoleculeCount()) {
		return &(_cells[_cellIteratorIndex].getParticles()[_moleculeIteratorIndex]);
	}
	else {
		++_cellIteratorIndex;
		_moleculeIteratorIndex = 0;
		while (_cellIteratorIndex < _cells.size()) {
			if (_cells[_cellIteratorIndex].getMoleculeCount() > 0) {
				return &(_cells[_cellIteratorIndex].getParticles()[0]);
			}
			++_cellIteratorIndex;
		}
	}
	return NULL;
}

Molecule* BlockedReorderedLinkedCells2::end() {
	return NULL;
}



void BlockedReorderedLinkedCells2::deleteOuterParticles() {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in BlockedReorderedLinkedCells (deleteOuterParticles) invalid, call update first" << endl;
		exit(1);
	}

	vector<unsigned long>::iterator cellIndexIter;
	for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter != _haloCellIndices.end(); cellIndexIter++) {
		BlockedCell& currentCell = _cells[*cellIndexIter];
		currentCell.removeAllParticles();
	}
}


double BlockedReorderedLinkedCells2::get_halo_L(int index) const {
	return _haloLength[index];
}


void BlockedReorderedLinkedCells2::getHaloParticles(list<Molecule*> &haloParticlePtrs) {
#ifndef NDEBUG
	global_log->error() << "TODO: method getHaloParticles is not implemented / refactored right!" << endl;
	exit(-1);
#endif

/*	if (_cellsValid == false) {
		global_log->error() << "Cell structure in BlockedReorderedLinkedCells (getHaloParticles) invalid, call update first" << endl;
		exit(1);
	}

	std::vector<Molecule*>::iterator particleIter;
	vector<unsigned long>::iterator cellIndexIter;

	// loop over all halo cells
	for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter != _haloCellIndices.end(); cellIndexIter++) {
		BlockedCell& currentCell = _cells[*cellIndexIter];
		// loop over all molecules in the cell
		for (particleIter = currentCell.getParticlePointers().begin(); particleIter != currentCell.getParticlePointers().end(); particleIter++) {
			haloParticlePtrs.push_back(*particleIter);
		}
	}
*/
}

void BlockedReorderedLinkedCells2::getRegion(double lowCorner[3], double highCorner[3], list<Molecule*> &particlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in BlockedReorderedLinkedCells (getRegion) invalid, call update first" << endl;
		exit(1);
	}

	int startIndex[3];
	int stopIndex[3];
	int globalCellIndex;

	for (int dim = 0; dim < 3; dim++) {
		if (lowCorner[dim] < this->_boundingBoxMax[dim] && highCorner[dim] > this->_boundingBoxMin[dim]) {
			startIndex[dim] = (int) floor((lowCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) - 1;
			stopIndex[dim] = (int) floor((highCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) + 1;
			if (startIndex[dim] < 0)
				startIndex[dim] = 0;
			if (stopIndex[dim] > _cellsPerDimension[dim] - 1)
				stopIndex[dim] = _cellsPerDimension[dim] - 1;
		}
		else {
			// No Part of the given region is owned by this process
			// --> chose some startIndex which is higher than the stopIndex
			startIndex[dim] = 1;
			stopIndex[dim] = 0;
		}
	}

	for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
		for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
			for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
				// globalCellIndex is the cellIndex of the molecule on the coarse Cell level.
				globalCellIndex = (iz * _cellsPerDimension[1] + iy) * _cellsPerDimension[0] + ix;
				// loop over all subcells (either 1 or 8)
				// traverse all molecules in the current cell
				MoleculeArray& particles = _cells[globalCellIndex].getParticles();
				for (size_t i = 0; i < particles.size(); i++) {
					if (particles[i].r(0) >= lowCorner[0] && particles[i].r(0) < highCorner[0] &&
						particles[i].r(1) >= lowCorner[1] && particles[i].r(1) < highCorner[1] &&
						particles[i].r(2) >= lowCorner[2] && particles[i].r(2) < highCorner[2]) {
						particlePtrs.push_back(&(particles[i]));
					}
				}
			}
		}
	}
}


void BlockedReorderedLinkedCells2::getRegion(double lowCorner[3], double highCorner[3], utils::DynamicArray<BasicMolecule, true, false>& particleArray) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in BlockedReorderedLinkedCells (getRegion) invalid, call update first" << endl;
		exit(1);
	}

	if (!IsSame<Molecule, BasicMolecule>::Result::value) {
		global_log->error() << "BlockedReorderedLinkedCells::getRegion(DynamicArray<BasicMolecule>&) must only be called for BasicMolecules!" << endl;
		return;
	}

	// cast to make this code compile for all combinations of Molecule/HandlerTypeMolecule
	MoleculeArray& particleArrayToReturn = reinterpret_cast<MoleculeArray&>(particleArray);

	int startIndex[3];
	int stopIndex[3];
	int globalCellIndex;

	for (int dim = 0; dim < 3; dim++) {
		if (lowCorner[dim] < this->_boundingBoxMax[dim] && highCorner[dim] > this->_boundingBoxMin[dim]) {
			startIndex[dim] = (int) floor((lowCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) - 1;
			stopIndex[dim] = (int) floor((highCorner[dim] - _haloBoundingBoxMin[dim]) / _cellLength[dim]) + 1;
			if (startIndex[dim] < 0)
				startIndex[dim] = 0;
			if (stopIndex[dim] > _cellsPerDimension[dim] - 1)
				stopIndex[dim] = _cellsPerDimension[dim] - 1;
		}
		else {
			// No Part of the given region is owned by this process
			// --> chose some startIndex which is higher than the stopIndex
			startIndex[dim] = 1;
			stopIndex[dim] = 0;
		}
	}

	// TODO Optimisation: for those cells which are in the inner of the requested area,
	//                    we can just append the whole particle array instead of pushing them back individually.
	for (int iz = startIndex[2]; iz <= stopIndex[2]; iz++) {
		for (int iy = startIndex[1]; iy <= stopIndex[1]; iy++) {
			for (int ix = startIndex[0]; ix <= stopIndex[0]; ix++) {
				// globalCellIndex is the cellIndex of the molecule on the coarse Cell level.
				globalCellIndex = (iz * _cellsPerDimension[1] + iy) * _cellsPerDimension[0] + ix;
				// loop over all subcells (either 1 or 8)
				// traverse all molecules in the current cell
				MoleculeArray& particles = _cells[globalCellIndex].getParticles();
				for (size_t i = 0; i < particles.size(); i++) {
					if (particles[i].r(0) >= lowCorner[0] && particles[i].r(0) < highCorner[0] &&
						particles[i].r(1) >= lowCorner[1] && particles[i].r(1) < highCorner[1] &&
						particles[i].r(2) >= lowCorner[2] && particles[i].r(2) < highCorner[2]) {
						particleArrayToReturn.push_back((particles[i]));
					}
				}
			}
		}
	}
}

void BlockedReorderedLinkedCells2::linearize(MoleculeArray& molecules, std::vector<int>& cellStartIndices) {
	int cellStartIndex = 0;
	for (size_t i = 0; i < _cells.size(); i++) {
		//std::cout << "Cell[" << i<< "].size()=" << _cells[i].getMoleculeCount() << endl;
		MoleculeArray& cellMolecules = _cells[i].getParticles();
		molecules.insert(molecules.end(), cellMolecules.begin(), cellMolecules.end());
		cellStartIndices.push_back(cellStartIndex);
		cellStartIndex += cellMolecules.size();
	}
}


void BlockedReorderedLinkedCells2::delinearize(MoleculeArray& molecules) {
	// TODO maybe replace the following with copying subarrays, instead of assignment?
	assert(molecules.size() == getNumberOfParticles());

	int moleculeIndex = 0;
	for (size_t i = 0; i < _cells.size(); i++) {
		MoleculeArray& cellMolecules = _cells[i].getParticles();
		for (size_t j = 0; j < cellMolecules.size(); j++) {
			assert(cellMolecules[j].id() == molecules[moleculeIndex].id());
			cellMolecules[j] = molecules[moleculeIndex];
			moleculeIndex++;
		}
	}
}


//################################################
//############ PRIVATE METHODS ###################
//################################################


void BlockedReorderedLinkedCells2::initializeCells() {
	_innerCellIndices.clear();
	_boundaryCellIndices.clear();
	_haloCellIndices.clear();
	unsigned long cellIndex;
	for (int iz = 0; iz < _cellsPerDimension[2]; ++iz) {
		for (int iy = 0; iy < _cellsPerDimension[1]; ++iy) {
			for (int ix = 0; ix < _cellsPerDimension[0]; ++ix) {
				cellIndex = cellIndexOf3DIndex(ix, iy, iz);
				if (ix < _haloWidthInNumCells[0] || iy < _haloWidthInNumCells[1] || iz < _haloWidthInNumCells[2] ||
				    ix >= _cellsPerDimension[0]-_haloWidthInNumCells[0] ||
				    iy >= _cellsPerDimension[1]-_haloWidthInNumCells[1] ||
				    iz >= _cellsPerDimension[2]-_haloWidthInNumCells[2]) {
					_cells[cellIndex].assingCellToHaloRegion();
					_haloCellIndices.push_back(cellIndex);
				}
				else if (ix < 2*_haloWidthInNumCells[0] || iy < 2*_haloWidthInNumCells[1] || iz < 2*_haloWidthInNumCells[2] ||
				         ix >= _cellsPerDimension[0]-2*_haloWidthInNumCells[0] ||
				         iy >= _cellsPerDimension[1]-2*_haloWidthInNumCells[1] ||
				         iz >= _cellsPerDimension[2]-2*_haloWidthInNumCells[2]) {
					_cells[cellIndex].assignCellToBoundaryRegion();
					_boundaryCellIndices.push_back(cellIndex);
				}
				else {
					_cells[cellIndex].assignCellToInnerRegion();
					_innerCellIndices.push_back(cellIndex);
				}
			}
		}
	}
}

void BlockedReorderedLinkedCells2::calculateNeighbourIndices() {
	_forwardNeighbourOffsets.clear();
	_backwardNeighbourOffsets.clear();
	double xDistanceSquare;
	double yDistanceSquare;
	double zDistanceSquare;
	int maxNeighbourOffset = 0;
	int minNeighbourOffset = 0;
	double cutoffRadiusSquare = pow(_cutoffRadius, 2);
	for (int zIndex = -_haloWidthInNumCells[2]; zIndex <= _haloWidthInNumCells[2]; zIndex++) {
		// The distance in one dimension is the width of a cell multiplied with the number
		// of cells between the two cells (this is received by substracting one of the
		// absolute difference of the cells, if this difference is not zero)
		if (zIndex != 0) {
			zDistanceSquare = pow((abs(zIndex) - 1) * _cellLength[2], 2);
		}
		else {
			zDistanceSquare = 0;
		}
		for (int yIndex = -_haloWidthInNumCells[1]; yIndex <= _haloWidthInNumCells[1]; yIndex++) {
			if (yIndex != 0) {
				yDistanceSquare = pow((abs(yIndex) - 1) * _cellLength[1], 2);
			}
			else {
				yDistanceSquare = 0;
			}
			for (int xIndex = -_haloWidthInNumCells[0]; xIndex <= _haloWidthInNumCells[0]; xIndex++) {
				if (xIndex != 0) {
					xDistanceSquare = pow((abs(xIndex) - 1) * _cellLength[0], 2);
				}
				else {
					xDistanceSquare = 0;
				}
				if (xDistanceSquare + yDistanceSquare + zDistanceSquare <= cutoffRadiusSquare) {
					long offset = cellIndexOf3DIndex(xIndex, yIndex, zIndex);
					if (offset > 0) {
						_forwardNeighbourOffsets.push_back(offset);
						if (offset > maxNeighbourOffset) {
							maxNeighbourOffset = offset;
						}
					}
					if (offset < 0) {
						_backwardNeighbourOffsets.push_back(offset);
						if (offset < minNeighbourOffset) {
							minNeighbourOffset = offset;
						}
					}
				}
			}
		}
	}

	#ifndef NDEBUG
	global_log->info() << "Neighbour offsets are bounded by "
	<< minNeighbourOffset << ", " << maxNeighbourOffset << endl;
#endif
	_blockTraverse.assignOffsets(_forwardNeighbourOffsets,
			_backwardNeighbourOffsets, maxNeighbourOffset, minNeighbourOffset);
}

unsigned long BlockedReorderedLinkedCells2::getCellIndexOfMolecule(Molecule* molecule) const {
	int cellIndex[3]; // 3D Cell index

	for (int dim = 0; dim < 3; dim++) {
		if (molecule->r(dim) < _haloBoundingBoxMin[dim] || molecule->r(dim) >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "getCellIndexOfMolecule(Molecule* molecule): Molecule is outside of the bounding box" << endl;
		}
		cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);

	}
	return this->cellIndexOf3DIndex( cellIndex[0], cellIndex[1], cellIndex[2] );
}

unsigned long BlockedReorderedLinkedCells2::cellIndexOf3DIndex(int xIndex, int yIndex, int zIndex) const {
	return (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
}


void BlockedReorderedLinkedCells2::deleteMolecule(unsigned long molid, double x, double y, double z) {

	int ix = (int) floor((x - this->_haloBoundingBoxMin[0]) / this->_cellLength[0]);
	int iy = (int) floor((y - this->_haloBoundingBoxMin[1]) / this->_cellLength[1]);
	int iz = (int) floor((z - this->_haloBoundingBoxMin[2]) / this->_cellLength[2]);

	unsigned long hash = this->cellIndexOf3DIndex(ix, iy, iz);
	if (hash >= _cells.size()) {
		global_log->error() << "coordinates for atom deletion lie outside bounding box." << endl;
		exit(1);
	}

	bool found = this->_cells[hash].deleteMolecule(molid);
	if (!found) {
		global_log->error() << "could not delete molecule " << molid << "." << endl;
		exit(1);
	}
}


int BlockedReorderedLinkedCells2::grandcanonicalBalance(DomainDecompBase* comm) {
	comm->collCommInit(1);
	comm->collCommAppendInt(this->_localInsertionsMinusDeletions);
	comm->collCommAllreduceSum();
	int universalInsertionsMinusDeletions = comm->collCommGetInt();
	comm->collCommFinalize();
	return universalInsertionsMinusDeletions;
}

/*****************************************************************/
/************** To be implementended / adapted *******************/
/*****************************************************************/

Molecule* BlockedReorderedLinkedCells2::deleteCurrent() {
#ifndef NDEBUG
	global_log->error() << "TODO: method deleteCurrent() is not implemented / refactored right!" << endl;
	exit(-1);
#endif
	return NULL;
}

double BlockedReorderedLinkedCells2::getEnergy(Molecule* m1) {
#ifndef NDEBUG
	global_log->error() << "TODO: getEnergy() method is not implemented / refactored right!" << endl;
	exit(-1);
#endif
	return 0;
}

/**
 * @todo replace this by a call to component->getNumMolecules() !?
 */
unsigned BlockedReorderedLinkedCells2::countParticles(unsigned int cid) {
	// I won't implement redundant stuff right now...
	#ifndef NDEBUG
		global_log->error() << "TODO: method countParticles(cid) is not implemented / refactored right!" << endl;
		exit(-1);
	#endif
		return 0;
}

/**
 * @todo move this method to the ChemicalPotential, using a call to ParticleContainer::getRegion() !?
 */
unsigned BlockedReorderedLinkedCells2::countParticles(unsigned int cid, double* cbottom, double* ctop) {
	#ifndef NDEBUG
		global_log->error() << "TODO: method countParticles(cid) is not implemented / refactored right!" << endl;
		exit(-1);
	#endif
	return 0;
}

void BlockedReorderedLinkedCells2::grandcanonicalStep(ChemicalPotential* mu, double T) {
#ifndef NDEBUG
	global_log->error() << "TODO: this method is not implemented / refactored right!" << endl;
	exit(-1);
#endif

}
