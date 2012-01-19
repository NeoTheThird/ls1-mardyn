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

#include "BlockedReorderedBlockTraverse.h"
#include "molecules/MoleculeTypes.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "BlockedCell.h"
#include "MemoryManager.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Common.h"

#include <vector>
#include <cmath>

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


BlockedReorderedBlockTraverse::BlockedReorderedBlockTraverse(
		ParticleContainer* moleculeContainer,
		vector<BlockedCell>& cells,
        vector<unsigned long>& innerCellIndices, 
        vector<unsigned long>& boundaryCellIndices,
        vector<unsigned long>& haloCellIndices,
		vector<vector<unsigned long> >& forwardNeighbourOffsets,
        vector<vector<unsigned long> >& backwardNeighbourOffsets
)
		: _moleculeContainer(moleculeContainer),
			_cells(cells),
			_innerCellIndices(innerCellIndices),
			_boundaryCellIndices(boundaryCellIndices),
            _haloCellIndices(haloCellIndices),
			_forwardNeighbourOffsets(&forwardNeighbourOffsets), _backwardNeighbourOffsets(&backwardNeighbourOffsets),
			_allocatedOffsets(false) {
}

BlockedReorderedBlockTraverse::BlockedReorderedBlockTraverse(
		ParticleContainer* moleculeContainer,
		vector<BlockedCell>& cells,
        vector<unsigned long>& innerCellIndices, 
        vector<unsigned long>& boundaryCellIndices, 
        vector<unsigned long>& haloCellIndices
)
		: _moleculeContainer(moleculeContainer),
			_cells(cells),
			_innerCellIndices(innerCellIndices),
			_boundaryCellIndices(boundaryCellIndices),
			_haloCellIndices(haloCellIndices),
			_forwardNeighbourOffsets(0), _backwardNeighbourOffsets(0),
			_allocatedOffsets(true)
{
	_forwardNeighbourOffsets = new vector<vector<unsigned long> >;
	_backwardNeighbourOffsets = new vector<vector<unsigned long> >;
}

BlockedReorderedBlockTraverse::~BlockedReorderedBlockTraverse() {
	if (_allocatedOffsets) {
		delete _forwardNeighbourOffsets;
		delete _backwardNeighbourOffsets;
	}
}

void BlockedReorderedBlockTraverse::assignOffsets(vector<unsigned long>& forwardNeighbourOffsets, vector<unsigned long>& backwardNeighbourOffsets,
		int maxNeighbourOffset, int minNeighbourOffset) {
	_forwardNeighbourOffsets->assign(_cells.size(), forwardNeighbourOffsets);
	_backwardNeighbourOffsets->assign(_cells.size(), backwardNeighbourOffsets);
	_maxNeighbourOffset = maxNeighbourOffset;
	_minNeighbourOffset = abs(minNeighbourOffset);
	global_log->info() << "BlockedReorderedBlockTraverse::assignOffsets() maxNeighbourOffsets=" << maxNeighbourOffset << "; minNeighbourOffsets=" << minNeighbourOffset << endl;
}

void BlockedReorderedBlockTraverse::traversePairs(ParticlePairsHandler* particlePairsHandler) {

	double _cutoffRadius = _moleculeContainer->getCutoff();
	double _LJCutoffRadius = _moleculeContainer->getLJCutoff();
	double _tersoffCutoffRadius = _moleculeContainer->getTersoffCutoff();
	vector<vector<unsigned long> >& forwardNeighbourOffsets = *_forwardNeighbourOffsets;
	vector<vector<unsigned long> >& backwardNeighbourOffsets = *_backwardNeighbourOffsets;

	particlePairsHandler->init();

	// XXX comment
	double distanceVector[3];
	// loop over all cells
	vector<BlockedCell>::iterator cellIter;
	MoleculeArray::iterator molIter1;
	MoleculeArray::iterator molIter2;

#ifndef NDEBUG
	// reset forces and momenta to zero
	global_log->debug() << "Resetting forces and momenta, disconnecting Tersoff pairs." << endl;
#endif
	{
		double zeroVec[3] = {0.0, 0.0, 0.0};

        // TODO: check if the reset is done twice as leaving this part has no difference on the result.
        //Molecule *moleculePtr;
        for (size_t i = 0; i < _cells.size(); i++) {
        	MoleculeArray& particles = _cells[i].getParticles();
        	for (size_t j = 0; j < particles.size(); j++) {
        		particles[j].setF(zeroVec);
        		particles[j].setM(zeroVec);
        		particles[j].clearTersoffNeighbourList();
        	}
        }
	}

	unsigned cellIndex;
	vector<unsigned long>::iterator cellIndexIter;
	vector<unsigned long>::iterator neighbourOffsetsIter;

	// sqare of the cutoff radius
	double cutoffRadiusSquare = _cutoffRadius * _cutoffRadius;
	double LJCutoffRadiusSquare = _LJCutoffRadius * _LJCutoffRadius;
	double tersoffCutoffRadiusSquare = _tersoffCutoffRadius * _tersoffCutoffRadius;

#ifndef NDEBUG
	global_log->debug() << "ReorderedBlockTraversal: Processing pairs and preprocessing Tersoff pairs." << endl;
	global_log->debug() << "_minNeighbourOffset=" << _minNeighbourOffset << "; _maxNeighbourOffset=" << _maxNeighbourOffset<< endl;
#endif


	// open the window of cells with cache activated
	for (unsigned int cellIndex = 0; cellIndex < _maxNeighbourOffset; cellIndex++) {
		#ifndef NDEBUG
		global_log->debug() << "Opening cached cells window for cell index= " << cellIndex
				<< " numMolecules()="<<_cells[cellIndex].getMoleculeCount() << endl;
		#endif
		if ( IsSame<Molecule,HandlerMoleculeType>::Result::value ) {
			_cells[cellIndex].setToHandlerMoleculeType();
		} else {
			_cells[cellIndex].convertToHandlerMoleculeType<Molecule, HandlerMoleculeType>();
		}
	}

	// loop over all inner cells and calculate forces to forward neighbours
	//for (cellIndexIter = _innerCellIndices.begin(); cellIndexIter != _innerCellIndices.end(); cellIndexIter++) {
	for (cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		BlockedCell& currentCell = _cells[cellIndex];

		// extend the window of cells with cache activated
		if (cellIndex + _maxNeighbourOffset < _cells.size()) {
			#ifndef NDEBUG
			global_log->debug() << "Opening cached cells window for cell index=" << (cellIndex + _maxNeighbourOffset)
					<< " with numMolecules()="<< _cells[cellIndex + _maxNeighbourOffset].getMoleculeCount()
					<< " currentCell " << cellIndex << endl;
			#endif
			if ( IsSame<Molecule,HandlerMoleculeType>::Result::value ) {
				_cells[cellIndex + _maxNeighbourOffset].setToHandlerMoleculeType();
			} else {
				_cells[cellIndex + _maxNeighbourOffset].convertToHandlerMoleculeType<Molecule, HandlerMoleculeType>();
			}

#ifdef PROFILE_MEMORY
			CProcessMemoryInformation::collectMaxUsage();
#endif

		}

		HandlerMoleculeTypeArray& currentCellParticles = currentCell.getHandlerTypeParticles();
		int currentParticleCount = currentCellParticles.size();

		// forces between molecules in the cell
		if (currentCell.isInnerCell()) {

#ifdef VECTORIZE
			MemoryManager::fp_memory_type** moleculePositions = currentCell.getMoleculePositons();
			MemoryManager::fp_memory_type** moleculeDistances = BlockedCell::memoryManager.getScratchMemory(currentParticleCount, currentParticleCount);
			global_log->debug() << "Calculate SSE currentParticleCount=" << currentParticleCount << endl;
			calculateDistances(currentParticleCount, moleculePositions, currentParticleCount, moleculePositions, moleculeDistances);
#endif

			for (int i = 0; i < currentParticleCount; i++) {
				HandlerMoleculeType& molecule1 = currentCellParticles[i];
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (int j = i+1; j < currentParticleCount; j++) {
 #ifdef VECTORIZE
					#ifndef NDEBUG
					double dd = currentCellParticles[j].dist2(molecule1, distanceVector);
					assert(fabs(dd - moleculeDistances[i][j]) < 0.0001);
					#endif
					if (moleculeDistances[i][j] < cutoffRadiusSquare) {
						HandlerMoleculeType& molecule2 = currentCellParticles[j];
						assert(&molecule1 != &molecule2);
						for (int iii = 0; iii < 3; iii++) {
							distanceVector[iii] = molecule1.r(iii) - molecule2.r(iii);
						}
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, moleculeDistances[i][j], (moleculeDistances[i][j] < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (moleculeDistances[i][j] < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
#else
					HandlerMoleculeType& molecule2 = currentCellParticles[j];
					assert(&molecule1 != &molecule2);
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
#endif
				}
			}

			// loop over all neighbours
			for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				BlockedCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				HandlerMoleculeTypeArray& neighbourCellParticles = neighbourCell.getHandlerTypeParticles();
				int neighbourParticleCount = neighbourCellParticles.size();

#ifdef VECTORIZE
				MemoryManager::fp_memory_type** currentMoleculePositions = currentCell.getMoleculePositons();
				MemoryManager::fp_memory_type** neighbourMoleculePositions = neighbourCell.getMoleculePositons();
				MemoryManager::fp_memory_type** moleculeDistances = BlockedCell::memoryManager.getScratchMemory(currentParticleCount, neighbourParticleCount);
				global_log->debug() << "Calculate SSE currentParticleCount=" << currentParticleCount << " neighbourParticleCount=" << neighbourParticleCount << endl;
				calculateDistances(currentParticleCount, currentMoleculePositions, neighbourParticleCount, neighbourMoleculePositions, moleculeDistances);
#endif

				// loop over all particles in the cell
				for (int i = 0; i < currentParticleCount; i++) {
					HandlerMoleculeType& molecule1 = currentCellParticles[i];
					unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

					for (int j = 0; j < neighbourParticleCount; j++) {
#ifdef VECTORIZE
					#ifndef NDEBUG
					double dd = neighbourCellParticles[j].dist2(molecule1, distanceVector);
					assert(fabs(dd - moleculeDistances[i][j]) < 0.0001);
					#endif
					if (moleculeDistances[i][j] < cutoffRadiusSquare) {
						HandlerMoleculeType& molecule2 = neighbourCellParticles[j];
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, moleculeDistances[i][j], (moleculeDistances[i][j] < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (moleculeDistances[i][j] < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
#else

						HandlerMoleculeType& molecule2 = neighbourCellParticles[j];
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < cutoffRadiusSquare) {
							particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
							if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
								particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
							}
						}
#endif
					}
				}
			}
		} // if (isInnerCell())

	// loop over halo cells and detect Tersoff neighbours within the halo
	// this is relevant for the angle summation
		if (currentCell.isHaloCell()) {
			for (int i = 0; i < currentParticleCount; i++ ) {
				HandlerMoleculeType& molecule1 = currentCellParticles[i];
				assert(molecule1.numTersoff() == 0);
				if (molecule1.numTersoff() == 0)
					continue;

				for (int j = i+1; j < currentParticleCount; j++) {
					HandlerMoleculeType& molecule2 = currentCellParticles[j];
					assert(&molecule1 != &molecule2);
					if (molecule2.numTersoff() > 0) {
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < tersoffCutoffRadiusSquare)
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
					}
				}

				for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
					int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
					if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) (_cells.size())))
						continue;

					BlockedCell& neighbourCell = _cells[neighbourCellIndex];
					if (!neighbourCell.isHaloCell())
						continue;

					HandlerMoleculeTypeArray& neighbourCellParticles = neighbourCell.getHandlerTypeParticles();
					int neighbourParticleCount = neighbourCellParticles.size();
					for (int j = 0; j < neighbourParticleCount; j++) {
						HandlerMoleculeType& molecule2 = neighbourCellParticles[j];
						if (molecule2.numTersoff() == 0)
							continue;
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < tersoffCutoffRadiusSquare)
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
					}
				}
			}
		}

	// loop over all boundary cells and calculate forces to forward and backward neighbours
		if (currentCell.isBoundaryCell()) {
#ifdef VECTORIZE
			MemoryManager::fp_memory_type** moleculePositions = currentCell.getMoleculePositons();
			MemoryManager::fp_memory_type** moleculeDistances = BlockedCell::memoryManager.getScratchMemory(currentParticleCount, currentParticleCount);
			global_log->debug() << "Calculate SSE currentParticleCount=" << currentParticleCount << endl;
			calculateDistances(currentParticleCount, moleculePositions, currentParticleCount, moleculePositions, moleculeDistances);
#endif
			for (int i = 0; i < currentParticleCount; i++) {
				HandlerMoleculeType& molecule1 = currentCellParticles[i];
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (int j = i+1; j < currentParticleCount; j++) {
#ifdef VECTORIZE
					#ifndef NDEBUG
					double dd = currentCellParticles[j].dist2(molecule1, distanceVector);
					assert(fabs(dd - moleculeDistances[i][j]) < 0.0001);
					#endif
					if (moleculeDistances[i][j] < cutoffRadiusSquare) {
						HandlerMoleculeType& molecule2 = currentCellParticles[j];
						assert(&molecule1 != &molecule2);
						for (int iii = 0; iii < 3; iii++) {
							distanceVector[iii] = molecule1.r(iii) - molecule2.r(iii);
						}
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, moleculeDistances[i][j], (moleculeDistances[i][j] < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (moleculeDistances[i][j] < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
#else
					HandlerMoleculeType& molecule2 = currentCellParticles[j];
					assert(&molecule1 != &molecule2);

					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
#endif
				}
			}

			// loop over all forward neighbours
			for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				BlockedCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				HandlerMoleculeTypeArray& neighbourCellParticles = neighbourCell.getHandlerTypeParticles();
				int neighbourParticleCount = neighbourCellParticles.size();

#ifdef VECTORIZE
				MemoryManager::fp_memory_type** currentMoleculePositions = currentCell.getMoleculePositons();
				MemoryManager::fp_memory_type** neighbourMoleculePositions = neighbourCell.getMoleculePositons();
				MemoryManager::fp_memory_type** moleculeDistances = BlockedCell::memoryManager.getScratchMemory(currentParticleCount, neighbourParticleCount);
				global_log->debug() << "Calculate SSE currentParticleCount=" << currentParticleCount << " neighbourParticleCount=" << neighbourParticleCount << endl;
				calculateDistances(currentParticleCount, currentMoleculePositions, neighbourParticleCount, neighbourMoleculePositions, moleculeDistances);
#endif

				// loop over all particles in the cell
				for (int i = 0; i < currentParticleCount; i++ ) {
					HandlerMoleculeType& molecule1 = currentCellParticles[i];
					unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

					for (int j = 0; j < neighbourParticleCount; j++) {
#ifdef VECTORIZE
					#ifndef NDEBUG
					double dd = neighbourCellParticles[j].dist2(molecule1, distanceVector);
					if (fabs(dd - moleculeDistances[i][j]) > 0.0001) {
						std::cout.precision(8);
						std::cout << "line 381: dd=" << dd << " sse distance=" << moleculeDistances[i][j]  << " i=" << i << " j=" << j << endl;
					}
					assert(fabs(dd - moleculeDistances[i][j]) < 0.0001);
					#endif
					if (moleculeDistances[i][j] < cutoffRadiusSquare) {
						HandlerMoleculeType& molecule2 = neighbourCellParticles[j];

						PairType pairType = MOLECULE_MOLECULE;
						if (neighbourCell.isHaloCell() && ! molecule1.isLessThan(molecule2)) {
							/* Do not sum up values twice. */
							pairType = MOLECULE_HALOMOLECULE;
						}
						for (int iii = 0; iii < 3; iii++) {
							distanceVector[iii] = molecule1.r(iii) - molecule2.r(iii);
						}
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, moleculeDistances[i][j], (moleculeDistances[i][j] < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (moleculeDistances[i][j] < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
						}
					}
#else
					HandlerMoleculeType& molecule2 = neighbourCellParticles[j];

					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						PairType pairType = MOLECULE_MOLECULE;
						if (neighbourCell.isHaloCell() && ! molecule1.isLessThan(molecule2)) {
							/* Do not sum up values twice. */
							pairType = MOLECULE_HALOMOLECULE;
						}
						particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
						}
					}
#endif
					}
				}
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for (neighbourOffsetsIter = backwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != backwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				BlockedCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				HandlerMoleculeTypeArray& neighbourCellParticles = neighbourCell.getHandlerTypeParticles();
				int neighbourParticleCount = neighbourCellParticles.size();


				if (neighbourCell.isHaloCell()) {
#ifdef VECTORIZE
					MemoryManager::fp_memory_type** currentMoleculePositions = currentCell.getMoleculePositons();
					MemoryManager::fp_memory_type** neighbourMoleculePositions = neighbourCell.getMoleculePositons();
					MemoryManager::fp_memory_type** moleculeDistances = BlockedCell::memoryManager.getScratchMemory(currentParticleCount, neighbourParticleCount);
					global_log->debug() << "Calculate SSE currentParticleCount=" << currentParticleCount << " neighbourParticleCount=" << neighbourParticleCount << endl;
					calculateDistances(currentParticleCount, currentMoleculePositions, neighbourParticleCount, neighbourMoleculePositions, moleculeDistances);
#endif
					// loop over all particles in the cell
					for (int i = 0; i < currentParticleCount; i++) {
						HandlerMoleculeType& molecule1 = currentCellParticles[i];
						unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

						for (int j = 0; j < neighbourParticleCount; j++) {
#ifdef VECTORIZE
						#ifndef NDEBUG
							double dd = neighbourCellParticles[j].dist2(molecule1, distanceVector);
							assert(fabs(dd - moleculeDistances[i][j]) < 0.0001);
						#endif

							if (moleculeDistances[i][j] < cutoffRadiusSquare) {
								HandlerMoleculeType& molecule2 = neighbourCellParticles[j];

								PairType pairType = molecule1.isLessThan(molecule2) ? MOLECULE_MOLECULE : MOLECULE_HALOMOLECULE;
								for (int iii = 0; iii < 3; iii++) {
									distanceVector[iii] = molecule1.r(iii) - molecule2.r(iii);
								}
								particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, moleculeDistances[i][j], (moleculeDistances[i][j] < LJCutoffRadiusSquare));
								if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (moleculeDistances[i][j] < tersoffCutoffRadiusSquare)) {
									particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
								}
							}
#else
							HandlerMoleculeType& molecule2 = neighbourCellParticles[j];

							double dd = molecule2.dist2(molecule1, distanceVector);
							if (dd < cutoffRadiusSquare) {
								PairType pairType = molecule1.isLessThan(molecule2) ? MOLECULE_MOLECULE : MOLECULE_HALOMOLECULE;
								particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < LJCutoffRadiusSquare));
								if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
									particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
								}
							}
#endif
						}
					}
				}
			}
		} // if ( isBoundaryCell() )

#ifndef NDEBUG
		global_log->debug() << "processing Tersoff potential." << endl;
#endif
		double params[15];
		double delta_r = 0.;
		bool knowparams = false;

		if (currentCell.isInnerCell() || currentCell.isBoundaryCell()) {
			for (int i = 0; i < currentParticleCount; i++) {
				HandlerMoleculeType& molecule1 = currentCellParticles[i];

				if (molecule1.numTersoff() == 0)
					continue;

				if (!knowparams) {
					delta_r = molecule1.tersoffParameters(params);
					knowparams = true;
				}
				particlePairsHandler->processTersoffAtom(molecule1, params, delta_r);
			}
		}

		// narrow the window of cells with cache activated
		if (cellIndex >= _minNeighbourOffset) {
#ifndef NDEBUG
			global_log->debug() << "Narrowing cached cells window for cell index=" << (cellIndex - _minNeighbourOffset)
					<< " with size()="<<_cells[cellIndex - _minNeighbourOffset].getMoleculeCount()
					<< " currentCell " << cellIndex << endl;
#endif
//			if (applyForces)
//				_cells[cellIndex + _minNeighbourOffset].applyForces();

			if ( IsSame<Molecule,HandlerMoleculeType>::Result::value ) {
				_cells[cellIndex - _minNeighbourOffset].setToMoleculeType();
			} else {
				_cells[cellIndex - _minNeighbourOffset].convertToMoleculeType<Molecule, HandlerMoleculeType>();
			}
		}

	} // for (cellIndex = 0; cellIndex < _cells.size(); cellIndex++)

	// close the window of cells with cache activated
	for (unsigned int cellIndex = _cells.size() - _minNeighbourOffset; cellIndex < _cells.size(); cellIndex++) {
#ifndef NDEBUG
			global_log->debug() << "Narrowing cached cells window for cell index=" << cellIndex
					<< " size()="<<_cells[cellIndex].getMoleculeCount() << endl;
#endif
			if ( IsSame<Molecule,HandlerMoleculeType>::Result::value ) {
				_cells[cellIndex].setToMoleculeType();
			} else {
				_cells[cellIndex].convertToMoleculeType<Molecule, HandlerMoleculeType>();
			}
	}

	particlePairsHandler->finish();
}
//}
