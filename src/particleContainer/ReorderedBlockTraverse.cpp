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

#include "ReorderedBlockTraverse.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "Cell.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"

#include <vector>
#include <cmath>

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


ReorderedBlockTraverse::ReorderedBlockTraverse(
		ParticleContainer* moleculeContainer,
		vector<Cell>& cells, 
        vector<unsigned long>& innerCellIndices, 
        vector<unsigned long>& boundaryCellIndices,
        vector<unsigned long>& haloCellIndices,
		vector<vector<unsigned long> >& forwardNeighbourOffsets,
        vector<vector<unsigned long> >& backwardNeighbourOffsets
)
		: _moleculeContainer(moleculeContainer),
			_particlePairsHandler(moleculeContainer->getPairHandler()),
			_cells(cells),
			_innerCellIndices(innerCellIndices),
			_boundaryCellIndices(boundaryCellIndices),
            _haloCellIndices(haloCellIndices),
			_forwardNeighbourOffsets(&forwardNeighbourOffsets), _backwardNeighbourOffsets(&backwardNeighbourOffsets),
			_allocatedOffsets(false) {
}

ReorderedBlockTraverse::ReorderedBlockTraverse(
		ParticleContainer* moleculeContainer,
		vector<Cell>& cells, 
        vector<unsigned long>& innerCellIndices, 
        vector<unsigned long>& boundaryCellIndices, 
        vector<unsigned long>& haloCellIndices
)
		: _moleculeContainer(moleculeContainer),
			_particlePairsHandler(moleculeContainer->getPairHandler()),
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

ReorderedBlockTraverse::~ReorderedBlockTraverse() {
	if (_allocatedOffsets) {
		delete _forwardNeighbourOffsets;
		delete _backwardNeighbourOffsets;
	}
}

void ReorderedBlockTraverse::assignOffsets(vector<unsigned long>& forwardNeighbourOffsets, vector<unsigned long>& backwardNeighbourOffsets,
		int maxNeighbourOffset, int minNeighbourOffset) {
	_forwardNeighbourOffsets->assign(_cells.size(), forwardNeighbourOffsets);
	_backwardNeighbourOffsets->assign(_cells.size(), backwardNeighbourOffsets);
	_maxNeighbourOffset = maxNeighbourOffset;
	_minNeighbourOffset = minNeighbourOffset;
}

void ReorderedBlockTraverse::traversePairs() {

	if (!IsSame<Molecule, HandlerMoleculeType>::Result::value) {
		global_log->error() << "For this implementation, Molecule and MoleculeHandlerType have to be " <<
				"Typedefs for the same class type!" << std::endl;
		exit(-1);
	}

	_particlePairsHandler = _moleculeContainer->getPairHandler();
	double _cutoffRadius = _moleculeContainer->getCutoff();
	double _LJCutoffRadius = _moleculeContainer->getLJCutoff();
	double _tersoffCutoffRadius = _moleculeContainer->getTersoffCutoff();
	vector<vector<unsigned long> >& forwardNeighbourOffsets = *_forwardNeighbourOffsets;
	vector<vector<unsigned long> >& backwardNeighbourOffsets = *_backwardNeighbourOffsets;

	_particlePairsHandler->init();

	// XXX comment
	double distanceVector[3];
	// loop over all cells
	vector<Cell>::iterator cellIter;
	vector<Molecule*>::iterator molIter1;
	vector<Molecule*>::iterator molIter2;

#ifndef NDEBUG
	// reset forces and momenta to zero
	global_log->debug() << "Resetting forces and momenta, disconnecting Tersoff pairs." << endl;
#endif
	{
		double zeroVec[3] = {0.0, 0.0, 0.0};

        // TODO: check if the reset is done twice as leaving this part has no difference on the result.
        //Molecule *moleculePtr;
        for (cellIter = _cells.begin(); cellIter != _cells.end(); ++cellIter) {
        	for (molIter1 = cellIter->getParticlePointers().begin(); molIter1 != cellIter->getParticlePointers().end(); molIter1++) {
        		Molecule& molecule1 = **molIter1;
        		molecule1.setF(zeroVec);
        		molecule1.setM(zeroVec);
        		molecule1.clearTersoffNeighbourList();
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
#endif

	// loop over all inner cells and calculate forces to forward neighbours
	for (cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		Cell& currentCell = _cells[cellIndex];
		std::vector<Molecule*>& currentCellParticles = currentCell.getParticlePointers();
		int currentParticleCount = currentCellParticles.size();

		// forces between molecules in the cell
		if (currentCell.isInnerCell()) {
			for (int i = 0; i < currentParticleCount; i++) {
				HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (int j = i+1; j < currentParticleCount; j++) {
					HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[j]);
					assert(&molecule1 != &molecule2);
					double dd = molecule2.dist2(molecule1, distanceVector);

					if (dd < cutoffRadiusSquare) {
						_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
				}
			}

			// loop over all neighbours
			for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				Cell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				std::vector<Molecule*>& neighbourCellParticles = neighbourCell.getParticlePointers();
				int neighbourParticleCount = neighbourCellParticles.size();

				// loop over all particles in the cell
				for (int i = 0; i < currentParticleCount; i++) {
					HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);
					unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

					for (int j = 0; j < neighbourParticleCount; j++) {
						HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(neighbourCellParticles[j]);
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < cutoffRadiusSquare) {
							_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
							if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
								_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
							}
						}
					}

				}
			}
		} // if (isInnerCell())

	// loop over halo cells and detect Tersoff neighbours within the halo
	// this is relevant for the angle summation
		if (currentCell.isHaloCell()) {
			for (int i = 0; i < currentParticleCount; i++) {
				HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);
				if (molecule1.numTersoff() == 0)
					continue;

				for (int j = i+1; j < currentParticleCount; j++) {
					HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[j]);
					assert(&molecule1 != &molecule2);
					if (molecule2.numTersoff() > 0) {
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < tersoffCutoffRadiusSquare)
							_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
					}
				}

				for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
					int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
					if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) (_cells.size())))
						continue;
					Cell& neighbourCell = _cells[neighbourCellIndex];
					if (!neighbourCell.isHaloCell())
						continue;

					std::vector<Molecule*>& neighbourCellParticles = neighbourCell.getParticlePointers();
					int neighbourParticleCount = neighbourCellParticles.size();
					for (int j = 0; j < neighbourParticleCount; j++) {
						HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(neighbourCellParticles[j]);
						if (molecule2.numTersoff() == 0)
							continue;
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < tersoffCutoffRadiusSquare)
							_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
					}
				}
			}
		}

	// loop over all boundary cells and calculate forces to forward and backward neighbours
		if (currentCell.isBoundaryCell()) {
			for (int i = 0; i < currentParticleCount; i++) {
				HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (int j = i+1; j < currentParticleCount; j++) {
					HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[j]);
					assert(&molecule1 != &molecule2);

					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
				}
			}

			// loop over all forward neighbours
			for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				Cell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				std::vector<Molecule*>& neighbourCellParticles = neighbourCell.getParticlePointers();
				int neighbourParticleCount = neighbourCellParticles.size();

				// loop over all particles in the cell
				for (int i = 0; i < currentParticleCount; i++) {
					HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);
					unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

					for (int j = 0; j < neighbourParticleCount; j++) {
						HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(neighbourCellParticles[j]);

						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < cutoffRadiusSquare) {
							PairType pairType = MOLECULE_MOLECULE;
							if (neighbourCell.isHaloCell() && ! molecule1.isLessThan(molecule2)) {
								/* Do not sum up values twice. */
								pairType = MOLECULE_HALOMOLECULE;
							}
							_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < LJCutoffRadiusSquare));
							if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
								_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
							}
						}
					}
				}
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for (neighbourOffsetsIter = backwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != backwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				Cell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				std::vector<Molecule*>& neighbourCellParticles = neighbourCell.getParticlePointers();
				int neighbourParticleCount = neighbourCellParticles.size();

				if (neighbourCell.isHaloCell()) {
					// loop over all particles in the cell
					for (int i = 0; i < currentParticleCount; i++) {
						HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);
						unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

						for (int j = 0; j < neighbourParticleCount; j++) {
							HandlerMoleculeType& molecule2 = *reinterpret_cast<HandlerMoleculeType*>(neighbourCellParticles[j]);

							double dd = molecule2.dist2(molecule1, distanceVector);
							if (dd < cutoffRadiusSquare) {
								PairType pairType = molecule1.isLessThan(molecule2) ? MOLECULE_MOLECULE : MOLECULE_HALOMOLECULE;
								_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < LJCutoffRadiusSquare));
								if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
									_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
								}
							}
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
				HandlerMoleculeType& molecule1 = *reinterpret_cast<HandlerMoleculeType*>(currentCellParticles[i]);

				if (molecule1.numTersoff() == 0)
					continue;

				if (!knowparams) {
					delta_r = molecule1.tersoffParameters(params);
					knowparams = true;
				}
				_particlePairsHandler->processTersoffAtom(molecule1, params, delta_r);
			}
		}
	} // for (cellIndex = 0; cellIndex < _cells.size(); cellIndex++)

	_particlePairsHandler->finish();
}
