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
#include <cmath>

#include "particleContainer/LinkedCells.h"

#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "ParticleCell.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "ensemble/GrandCanonical.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


LinkedCells::LinkedCells(
		double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double LJCutoffRadius,
		double cellsInCutoffRadius
)
		: ParticleContainer(bBoxMin, bBoxMax), _cutoffRadius(cutoffRadius),
	      _LJCutoffRadius(LJCutoffRadius), _cellsInCutoffRadius(cellsInCutoffRadius)
{
	global_log->debug() << "cutoff: " << cutoffRadius << endl;
	global_log->debug() << "LJ cutoff:" << LJCutoffRadius << endl;
	global_log->debug() << "# cells in cutoff: " << cellsInCutoffRadius << endl;

	build(bBoxMin, bBoxMax);

	this->_localInsertionsMinusDeletions = 0;
}


LinkedCells::~LinkedCells() {
}



void LinkedCells::build(double bBoxMin[3], double bBoxMax[3]) {
	for (int i = 0; i < 3; i++) {
		this->_boundingBoxMin[i] = bBoxMin[i];
		this->_boundingBoxMax[i] = bBoxMax[i];
	}

	int numberOfCells = 1;
	int boxWidthInNumCells[3];

	for (int d = 0; d < 3; d++) {
		/* first calculate the cell length for this dimension */
		boxWidthInNumCells[d] = floor((_boundingBoxMax[d] - _boundingBoxMin[d]) / _cutoffRadius * _cellsInCutoffRadius);
		// in each dimension at least one layer of (inner+boundary) cells is necessary
		if( boxWidthInNumCells[d] == 0 ) {
			boxWidthInNumCells[d] = 1;
		}
		_cellLength[d] = (_boundingBoxMax[d] - _boundingBoxMin[d]) / boxWidthInNumCells[d];
		_haloWidthInNumCells[d] = ceil(_cellsInCutoffRadius);
		_haloLength[d] = _haloWidthInNumCells[d] * _cellLength[d];
		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _haloLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _haloLength[d];

		_cellsPerDimension[d] = boxWidthInNumCells[d] + 2 * _haloWidthInNumCells[d];

		numberOfCells *= _cellsPerDimension[d];
		assert(numberOfCells > 0);
	}
	global_log->debug() << "Cell size (" << _cellLength[0] << ", " << _cellLength[1] << ", " << _cellLength[2] << ")" << endl;
	global_log->debug() << "Cells per dimension (incl. halo): " << _cellsPerDimension[0] << " x " << _cellsPerDimension[1] << " x " << _cellsPerDimension[2] << endl;

	_cells.resize(numberOfCells);

// If the width of the inner region is less than the width of the halo
	// region a parallelisation is not possible (with the used algorithms).
	// If a particle leaves this box, it would need to be communicated to the two next neighbours.
	if (boxWidthInNumCells[0] < 2* _haloWidthInNumCells[0] ||
			boxWidthInNumCells[1] < 2* _haloWidthInNumCells[1] ||
			boxWidthInNumCells[2] < 2* _haloWidthInNumCells[2]) {
		global_log->error() << "LinkedCells (constructor): bounding box too small for calculated cell length" << endl;
		global_log->error() << "_cellsPerDimension" << _cellsPerDimension[0] << " / " << _cellsPerDimension[1] << " / " << _cellsPerDimension[2] << endl;
		global_log->error() << "_haloWidthInNumCells" << _haloWidthInNumCells[0] << " / " << _haloWidthInNumCells[1] << " / " << _haloWidthInNumCells[2] << endl;
		exit(5);
	}

	initializeCells();
	calculateNeighbourIndices();

	// TODO: We loose particles here as they are not communicated to the new owner
	// WE: I don't think so - in KDDecomposition, particles are communicated and kept
	//     in buffers, then the container is rebuilt, then they are inserted.
	// delete all Particles which are outside of the halo region
	std::list<Molecule>::iterator particleIterator = _particles.begin();
	bool erase_mol;
	while (particleIterator != _particles.end()) {
		erase_mol = false;
		for (unsigned short d = 0; d < 3; ++d) {
			const double& rd = particleIterator->r(d);
			// The molecules has to be within the domain of the process
			// If it is outside in at least one dimension, it has to be
			// erased /
			if (rd < this->_haloBoundingBoxMin[d] || rd >= this->_haloBoundingBoxMax[d])
				erase_mol = true;
		}
		if (erase_mol) {
			particleIterator = _particles.erase(particleIterator);
		}
		else {
			particleIterator++;
		}
	}
	_cellsValid = false;
}


void LinkedCells::update() {
	// clear all Cells
	std::vector<ParticleCell>::iterator celliter;
	for (celliter = (_cells).begin(); celliter != (_cells).end(); ++celliter) {
		(*celliter).removeAllParticles();
	}

	std::list<Molecule>::iterator pos;
	for (pos = _particles.begin(); pos != _particles.end(); ++pos) {
		// determine the cell into which the particle belongs
		Molecule &m = *pos;
		unsigned long index = getCellIndexOfMolecule(&m);
		_cells[index].addParticle(&(*pos));
	}
	_cellsValid = true;
}

void LinkedCells::addParticle(Molecule& particle) {

	double x = particle.r(0);
	double y = particle.r(1);
	double z = particle.r(2);

	if ( ( x >= _haloBoundingBoxMin[0]) && (x < _haloBoundingBoxMax[0]) &&
	     ( y >= _haloBoundingBoxMin[1]) && (y < _haloBoundingBoxMax[1]) &&
	     ( z >= _haloBoundingBoxMin[2]) && (z < _haloBoundingBoxMax[2]) ) {

		_particles.push_front( particle );
		/* TODO: Have a closer look onto this check as there is no warning or error message.
		 *
		 * I (WE) guess this should be a performance optimization: the particle is added into this
		 * container anyway, it is just not sorted into the cells. But as the container is not valid,
		 * update() has to be called anyway.
		 */
		if (_cellsValid) {
			int cellIndex = getCellIndexOfMolecule(&particle);
			_cells[cellIndex].addParticle(&(_particles.front()));
		}
	}
}



/**
 * @todo replace this by a call to component->getNumMolecules() !?
 */
unsigned LinkedCells::countParticles(unsigned int cid) {
	unsigned N = 0;
	std::vector<Molecule*>::iterator molIter1;
	for (unsigned i = 0; i < _cells.size(); i++) {
		ParticleCell& currentCell = _cells[i];
		if( !currentCell.isHaloCell() ) {
			for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {

				if ((*molIter1)->componentid() == cid)
					N++;
			}
		}
	}
	return N;
}

/**
 * @todo move this method to the ChemicalPotential, using a call to ParticleContainer::getRegion() !?
 */
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
	std::vector<Molecule*>::iterator molIter1;
	bool individualCheck;
	int cellid;

	for (cix[0] = minIndex[0]; maxIndex[0] >= cix[0]; (cix[0])++) {
		for (cix[1] = minIndex[1]; maxIndex[1] >= cix[1]; (cix[1])++) {
			for (cix[2] = minIndex[2]; maxIndex[2] >= cix[2]; (cix[2])++) {
				individualCheck = (cix[0] == minIndex[0]) || (cix[0] == minIndex[0] + 1) ||
				                  (cix[0] == maxIndex[0]) || (cix[0] == maxIndex[0] - 1) ||
				                  (cix[1] == minIndex[1]) || (cix[1] == minIndex[1] + 1) ||
				                  (cix[1] == maxIndex[1]) || (cix[1] == maxIndex[1] - 1) ||
				                  (cix[2] == minIndex[2]) || (cix[2] == minIndex[2] + 1) ||
				                  (cix[2] == maxIndex[2]) || (cix[2] == maxIndex[2] - 1);
			cellid = this->cellIndexOf3DIndex(cix[0], cix[1], cix[2]);
				ParticleCell& currentCell = _cells[cellid];
				if (currentCell.isHaloCell())
					continue;
				if (individualCheck) {
					for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
						if (((*molIter1)->r(0) > cbottom[0]) &&
						    ((*molIter1)->r(1) > cbottom[1]) &&
						    ((*molIter1)->r(2) > cbottom[2]) &&
						    ((*molIter1)->r(0) < ctop[0]) &&
						    ((*molIter1)->r(1) < ctop[1]) &&
						    ((*molIter1)->r(2) < ctop[2]) &&
						    ((*molIter1)->componentid() == cid)) {
							N++;
						}
					}
				}
				else {
					for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
						if ((*molIter1)->componentid() == cid)
							N++;
					}
				}
			}
		}
	}

	return N;
}


void LinkedCells::traverseCells(CellProcessor& cellProcessor) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (traversePairs) invalid, call update first" << endl;
		exit(1);
	}

	vector<unsigned long>::iterator neighbourOffsetsIter;

#ifndef NDEBUG
	global_log->debug() << "LinkedCells::traverseCells: Processing pairs and preprocessing Tersoff pairs." << endl;
	global_log->debug() << "_minNeighbourOffset=" << _minNeighbourOffset << "; _maxNeighbourOffset=" << _maxNeighbourOffset<< endl;
#endif

	cellProcessor.initTraversal(_maxNeighbourOffset + _minNeighbourOffset +1);
	// open the window of cells activated
	for (unsigned int cellIndex = 0; cellIndex < _maxNeighbourOffset; cellIndex++) {
		#ifndef NDEBUG
		global_log->debug() << "Open cells window for cell index= " << cellIndex
				<< " numMolecules()="<<_cells[cellIndex].getMoleculeCount() << endl;
		#endif
		cellProcessor.preprocessCell(_cells[cellIndex]);
	}

	// loop over all inner cells and calculate forces to forward neighbours
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		ParticleCell& currentCell = _cells[cellIndex];

		// extend the window of cells with cache activated
		if (cellIndex + _maxNeighbourOffset < _cells.size()) {
			#ifndef NDEBUG
			global_log->debug() << "Opening cells window for cell index=" << (cellIndex + _maxNeighbourOffset)
					<< " with numMolecules()="<< _cells[cellIndex + _maxNeighbourOffset].getMoleculeCount()
					<< " currentCell " << cellIndex << endl;
			#endif
			cellProcessor.preprocessCell(_cells[cellIndex + _maxNeighbourOffset]);
		}

		if (currentCell.isInnerCell()) {
			cellProcessor.processCell(currentCell);
			// loop over all neighbours
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		if (currentCell.isHaloCell()) {
			cellProcessor.processCell(currentCell);
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
				if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) (_cells.size())))
					continue;
				ParticleCell& neighbourCell = _cells[neighbourCellIndex];
				if (!neighbourCell.isHaloCell())
					continue;

				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		// loop over all boundary cells and calculate forces to forward and backward neighbours
		if (currentCell.isBoundaryCell()) {
			cellProcessor.processCell(currentCell);

			// loop over all forward neighbours
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for (neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];
				if (neighbourCell.isHaloCell()) {
					cellProcessor.processCellPair(currentCell, neighbourCell);
				}
			}
		} // if ( isBoundaryCell() )

		// narrow the window of cells activated
		if (cellIndex >= _minNeighbourOffset) {
#ifndef NDEBUG
			global_log->debug() << "Narrowing cells window for cell index=" << (cellIndex - _minNeighbourOffset)
									<< " with size()="<<_cells[cellIndex - _minNeighbourOffset].getMoleculeCount()
									<< " currentCell " << cellIndex << endl;
#endif
			cellProcessor.postprocessCell(_cells[cellIndex - _minNeighbourOffset]);
		}
	} // loop over all cells

	// close the window of cells with cache activated
	for (unsigned int cellIndex = _cells.size() - _minNeighbourOffset; cellIndex < _cells.size(); cellIndex++) {
#ifndef NDEBUG
			global_log->debug() << "Narrowing cells window for cell index=" << cellIndex
					<< " size()="<<_cells[cellIndex].getMoleculeCount() << endl;
#endif
			cellProcessor.postprocessCell(_cells[cellIndex]);
	}
	cellProcessor.endTraversal();
}

unsigned long LinkedCells::getNumberOfParticles() {
	return _particles.size();
}

Molecule* LinkedCells::begin() {
	_particleIter = _particles.begin();
	if (_particleIter != _particles.end()) {
		return &(*_particleIter);
	}
	else {
		return NULL;
	}
}

Molecule* LinkedCells::next() {
	_particleIter++;
	if (_particleIter != _particles.end()) {
		return &(*_particleIter);
	}
	else {
		return NULL;
	}
}

Molecule* LinkedCells::end() {
	return NULL;
}

Molecule* LinkedCells::deleteCurrent() {
	_particleIter = _particles.erase(_particleIter);
	if (_particleIter != _particles.end()) {
		return &(*_particleIter);
	}
	else {
		return NULL;
	}
}

void LinkedCells::deleteOuterParticles() {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (deleteOuterParticles) invalid, call update first" << endl;
		exit(1);
	}

	vector<unsigned long>::iterator cellIndexIter;
	//std::list<Molecule*>::iterator molIter1;
	for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter != _haloCellIndices.end(); cellIndexIter++) {
		ParticleCell& currentCell = _cells[*cellIndexIter];
		currentCell.removeAllParticles();
	}

	std::list<Molecule>::iterator particleIterator = _particles.begin();
	bool erase_mol;
	while (particleIterator != _particles.end()) {
		erase_mol = false;
		for (unsigned short d = 0; d < 3; ++d) {
			const double& rd = particleIterator->r(d);
			// The molecules has to be within the domain of the process
			// If it is outside in at least one dimension, it has to be
			// erased /
			if (rd < this->_boundingBoxMin[d] || rd >= this->_boundingBoxMax[d])
				erase_mol = true;
		}
		if (erase_mol) {
			particleIterator = _particles.erase(particleIterator);
		}
		else {
			particleIterator++;
		}
	}
}

double LinkedCells::get_halo_L(int index) const {
	return _haloLength[index];
}


void LinkedCells::getHaloParticles(list<Molecule*> &haloParticlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (getHaloParticles) invalid, call update first" << endl;
		exit(1);
	}

	std::vector<Molecule*>::iterator particleIter;
	vector<unsigned long>::iterator cellIndexIter;

	// loop over all halo cells
	for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter != _haloCellIndices.end(); cellIndexIter++) {
		ParticleCell& currentCell = _cells[*cellIndexIter];
		// loop over all molecules in the cell
		for (particleIter = currentCell.getParticlePointers().begin(); particleIter != currentCell.getParticlePointers().end(); particleIter++) {
			haloParticlePtrs.push_back(*particleIter);
		}
	}
}

void LinkedCells::getRegion(double lowCorner[3], double highCorner[3], list<Molecule*> &particlePtrs) {
	if (_cellsValid == false) {
		global_log->error() << "Cell structure in LinkedCells (getRegion) invalid, call update first" << endl;
		exit(1);
	}

	int startIndex[3];
	int stopIndex[3];
	int globalCellIndex;
	std::vector<Molecule*>::iterator particleIter;

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
				for (particleIter = _cells[globalCellIndex].getParticlePointers().begin(); particleIter != _cells[globalCellIndex].getParticlePointers().end(); particleIter++) {
					if ((*particleIter)->r(0) >= lowCorner[0] && (*particleIter)->r(0) < highCorner[0] &&
					    (*particleIter)->r(1) >= lowCorner[1] && (*particleIter)->r(1) < highCorner[1] &&
					    (*particleIter)->r(2) >= lowCorner[2] && (*particleIter)->r(2) < highCorner[2]) {
						particlePtrs.push_back(*particleIter);
					}
				}
			}
		}
	}
}

//################################################
//############ PRIVATE METHODS ###################
//################################################


void LinkedCells::initializeCells() {
	_innerCellIndices.clear();
	_boundaryCellIndices.clear();
	_haloCellIndices.clear();
	unsigned long cellIndex;
	for (int iz = 0; iz < _cellsPerDimension[2]; ++iz) {
		for (int iy = 0; iy < _cellsPerDimension[1]; ++iy) {
			for (int ix = 0; ix < _cellsPerDimension[0]; ++ix) {

				cellIndex = cellIndexOf3DIndex(ix, iy, iz);
				_cells[cellIndex].skipCellFromHaloRegion();
				_cells[cellIndex].skipCellFromBoundaryRegion();
				_cells[cellIndex].skipCellFromInnerRegion();

				if (ix < _haloWidthInNumCells[0] || iy < _haloWidthInNumCells[1] || iz < _haloWidthInNumCells[2] ||
				    ix >= _cellsPerDimension[0]-_haloWidthInNumCells[0] ||
				    iy >= _cellsPerDimension[1]-_haloWidthInNumCells[1] ||
				    iz >= _cellsPerDimension[2]-_haloWidthInNumCells[2]) {
					_cells[cellIndex].assignCellToHaloRegion();
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

void LinkedCells::calculateNeighbourIndices() {
	global_log->debug() << "Setting up cell neighbour indice lists." << endl;
	_forwardNeighbourOffsets.clear();
	_backwardNeighbourOffsets.clear();
	_maxNeighbourOffset = 0;
	_minNeighbourOffset = 0;
	double xDistanceSquare;
	double yDistanceSquare;
	double zDistanceSquare;
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
					long int offset = cellIndexOf3DIndex(xIndex, yIndex, zIndex);
					if (offset > 0) {
						_forwardNeighbourOffsets.push_back(abs(offset));
						if (offset > _maxNeighbourOffset) {
							_maxNeighbourOffset = offset;
						}
					}
					if (offset < 0) {
						_backwardNeighbourOffsets.push_back(abs(offset));
						if (abs(offset) > _minNeighbourOffset) {
							_minNeighbourOffset = abs(offset);
						}
					}
				}
			}
		}
	}

	global_log->info() << "Neighbour offsets are bounded by "
			<< _minNeighbourOffset << ", " << _maxNeighbourOffset << endl;
}

unsigned long int LinkedCells::getCellIndexOfMolecule(Molecule* molecule) const {
	int cellIndex[3]; // 3D Cell index

	for (int dim = 0; dim < 3; dim++) {
#ifndef NDEBUG
		if (molecule->r(dim) < _haloBoundingBoxMin[dim] || molecule->r(dim) >= _haloBoundingBoxMax[dim]) {
			global_log->error() << "Molecule is outside of bounding box" << endl;
			global_log->debug() << "Molecule:\n" << *molecule << endl;
		}
#endif
		cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);

	}
	return this->cellIndexOf3DIndex( cellIndex[0], cellIndex[1], cellIndex[2] );
}

long int LinkedCells::cellIndexOf3DIndex(int xIndex, int yIndex, int zIndex) const {
	return (zIndex * _cellsPerDimension[1] + yIndex) * _cellsPerDimension[0] + xIndex;
}


void LinkedCells::deleteMolecule(unsigned long molid, double x, double y, double z) {

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


double LinkedCells::getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1) {

	double u = 0.0;
	std::vector<Molecule*>::iterator molIter2;
	vector<unsigned long>::iterator neighbourOffsetsIter;

	// sqare of the cutoffradius
	double cutoffRadiusSquare = _cutoffRadius * _cutoffRadius;
	double LJCutoffRadiusSquare = _LJCutoffRadius * _LJCutoffRadius;
	double dd;
	double distanceVector[3];

	unsigned long cellIndex = getCellIndexOfMolecule(m1);
	ParticleCell& currentCell = _cells[cellIndex];

	if (m1->numTersoff() > 0) {
		global_log->error() << "The grand canonical ensemble is not implemented for solids." << endl;
		exit(484);
	}
	// molecules in the cell
	for (molIter2 = currentCell.getParticlePointers().begin(); molIter2 != currentCell.getParticlePointers().end(); molIter2++) {
		//cout<<"m1 id "<<m1->id()<<" m2 id"<<(*molIter2)->id()<<endl;
		if (m1->id() == (*molIter2)->id())
			continue;
		dd = (*molIter2)->dist2(*m1, distanceVector);
		if (dd > cutoffRadiusSquare)
			continue;
		
		u += particlePairsHandler->processPair(*m1, **molIter2, distanceVector, MOLECULE_MOLECULE_FLUID, dd, (dd < LJCutoffRadiusSquare));
	}

	// backward and forward neighbours
	for (neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
		if (neighbourOffsetsIter == _backwardNeighbourOffsets.end())
			neighbourOffsetsIter = _forwardNeighbourOffsets.begin();

		ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
		for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2 != neighbourCell.getParticlePointers().end(); molIter2++) {

			dd = (*molIter2)->dist2(*m1, distanceVector);
			//dd = cutoffRadiusSquare - 10;
			if (dd > cutoffRadiusSquare)
				continue;
			u += particlePairsHandler->processPair(*m1, **molIter2, distanceVector, MOLECULE_MOLECULE_FLUID, dd, (dd < LJCutoffRadiusSquare));
		}
	}
	return u;
}

int LinkedCells::grandcanonicalBalance(DomainDecompBase* comm) {
	comm->collCommInit(1);
	comm->collCommAppendInt(this->_localInsertionsMinusDeletions);
	comm->collCommAllreduceSum();
	int universalInsertionsMinusDeletions = comm->collCommGetInt();
	comm->collCommFinalize();
	return universalInsertionsMinusDeletions;
}

void LinkedCells::grandcanonicalStep(ChemicalPotential* mu, double T, Domain* domain) {
	bool accept = true;
	double DeltaUpot;
	Molecule* m;
	ParticlePairs2PotForceAdapter particlePairsHandler(*domain);

	this->_localInsertionsMinusDeletions = 0;

	mu->submitTemperature(T);
	double minco[3];
	double maxco[3];
	for (int d = 0; d < 3; d++) {
		minco[d] = this->getBoundingBoxMin(d);
		maxco[d] = this->getBoundingBoxMax(d);
	}

	bool hasDeletion = true;
	bool hasInsertion = true;
	double ins[3];
	unsigned nextid = 0;
	while (hasDeletion || hasInsertion) {
		if (hasDeletion)
			hasDeletion = mu->getDeletion(this, minco, maxco);
		if (hasDeletion) {
			m = &(*(this->_particleIter));
			DeltaUpot = -1.0 * getEnergy(&particlePairsHandler, m);

			accept = mu->decideDeletion(DeltaUpot / T);
#ifndef NDEBUG
			if(accept) global_log->debug() << "r" << mu->rank() << "d" << m->id() << endl;
			else global_log->debug() << "   (r" << mu->rank() << "-d" << m->id() << ")" << endl;
#endif
			if (accept) {
				m->upd_cache();
				// reset forces and momenta to zero
				{
					double zeroVec[3] = {0.0, 0.0, 0.0};
					m->setF(zeroVec);
					m->setM(zeroVec);
				}
				mu->storeMolecule(*m);
				this->deleteMolecule(m->id(), m->r(0), m->r(1), m->r(2));
				this->_particles.erase(this->_particleIter);
				this->_particleIter = _particles.begin();
				this->_localInsertionsMinusDeletions--;
			}
		}

		if (mu->isWidom()){
			m = &(*(_particles.begin()));
			mu->storeMolecule(*m);
		}
		if (hasInsertion) {
			nextid = mu->getInsertion(ins);
			hasInsertion = (nextid > 0);
		}
		if (hasInsertion) {
			Molecule tmp = mu->loadMolecule();
			for (int d = 0; d < 3; d++)
				tmp.setr(d, ins[d]);
			tmp.setid(nextid);
			this->_particles.push_back(tmp);

			std::list<Molecule>::iterator mit = _particles.end();
			mit--;
			m = &(*mit);
			m->upd_cache();
			// reset forces and momenta to zero
			{
				double zeroVec[3] = {0.0, 0.0, 0.0};
				m->setF(zeroVec);
				m->setM(zeroVec);
			}
			m->check(nextid);
#ifndef NDEBUG
			global_log->debug() << "rank " << mu->rank() << ": insert " << m->id()
			<< " at the reduced position (" << ins[0] << "/" << ins[1] << "/" << ins[2] << ")? " << endl;
#endif

			unsigned long cellid = this->getCellIndexOfMolecule(m);
			this->_cells[cellid].addParticle(m);
			DeltaUpot = getEnergy(&particlePairsHandler, m);
                        domain->submitDU(mu->getComponentID(), DeltaUpot, ins);
			accept = mu->decideInsertion(DeltaUpot / T);

#ifndef NDEBUG
			if(accept) global_log->debug() << "r" << mu->rank() << "i" << mit->id() << ")" << endl;
			else global_log->debug() << "   (r" << mu->rank() << "-i" << mit->id() << ")" << endl;
#endif
			if (accept) {
				this->_localInsertionsMinusDeletions++;
			}
			else {
				// this->deleteMolecule(m->id(), m->r(0), m->r(1), m->r(2));
				this->_cells[cellid].deleteMolecule(m->id());

				mit->check(m->id());
				this->_particles.erase(mit);
			}
		}
	}
	for (m = this->begin(); m != this->end(); m = this->next()) {
#ifndef NDEBUG
		m->check(m->id());
#endif
	}
}
