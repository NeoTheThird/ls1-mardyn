#ifndef CELL_H_
#define CELL_H_

#include "utils/DynamicArray.h"
#include "molecules/MoleculeTypes.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"

#include "MemoryManager.h"
#include "Simulation.h"
#include "Domain.h"


class HybridCell {

public:

	static MemoryManager memoryManager;

private:

	//! the basic molecules contained in this cell
	std::vector<Molecule*> _particlePointers;

	//! the molecules as the handler expects them
	HandlerMoleculeTypeArray* _handlerParticles;

	// contains the x, y and z coordinates of the molecules, if converted to HandlerMoleculeType.
	MemoryManager::fp_memory_type* _moleculePositions[3];

	//! true when the cell is in the halo region
	bool _haloCellState;
	//! true when the cell is in the boundary region
	bool _boundaryCellState;
	//! true when the cell is in the inner region
	bool _innerCellState;



#ifndef NDEBUG
	//! indicate what is the current internal representation of the particles
	//! this will be useful for debugging
	enum CurrentMoleculeType {BasicMolecule, HandlerMolecule};

	CurrentMoleculeType _currentMoleculeType;
#endif

public:
	//! The constructor sets all cell-states to false 
	HybridCell();

	HybridCell(const HybridCell& other);

	~HybridCell();

	//! removes all elements from the list molecules
	void removeAllParticles();

	//! insert a single molecule into this cell
	void addParticle(Molecule* particle);

	//! return a reference to the list of molecules in this cell
	std::vector<Molecule*>& getParticlePointers();

	//! return a reference to the list of molecules in this cell
	//! in the representation for the ParticlePairsHandler
	HandlerMoleculeTypeArray& getHandlerTypeParticles();

	//! @todo Return type bool neccessary!?
	bool deleteMolecule(unsigned long molid);

	//! Set the flag for a Halo Cell
	void assingCellToHaloRegion();

	//! Set the flag for a Boundary Cell
	void assignCellToBoundaryRegion();

	//! Set the flag for a Inner Cell
	void assignCellToInnerRegion();

	//! returns true, if the cell is a Halo Cell, otherwise false
	bool isHaloCell() const;

	//! returns true, if the cell is a Boundary Cell, otherwise false
	bool isBoundaryCell() const;

	//! returns true, if the cell is a Inner Cell, otherwise false
	bool isInnerCell() const;

	//! return the number of molecules contained in this cell
	int getMoleculeCount() const;

	/**
	 * convert the internal representation of the molecules to the handlerMoleculeType.
	 * After the conversion, only that type should be used.
	 */
	template <class Molecule, class HandlerMoleculeType>
	void convertToHandlerMoleculeType();

	/**
	 * convert the internal representation of the molecules to simple molecules.
	 * After the conversion, only that type should be used.
	 */
	template <class Molecule, class HandlerMolecule>
	void convertToMoleculeType();


	MemoryManager::fp_memory_type** getMoleculePositons() {
#ifndef VECTORIZE
		assert(false);
#endif
		assert(_currentMoleculeType == HandlerMolecule);
		return _moleculePositions;
	}
};



template <class Molecule, class HandlerMoleculeType>
void HybridCell::convertToHandlerMoleculeType() {

	#ifndef NDEBUG
	assert(_currentMoleculeType == BasicMolecule);
	_currentMoleculeType = HandlerMolecule;
	#endif

	std::vector<Component>& components = global_simulation->getDomain()->getComponents();

	_handlerParticles = memoryManager.getMoleculeArray();
#ifdef VECTORIZE
	_moleculePositions[0] = memoryManager.getFPMemory(_particlePointers.size());
	_moleculePositions[1] = memoryManager.getFPMemory(_particlePointers.size());
	_moleculePositions[2] = memoryManager.getFPMemory(_particlePointers.size());
#endif

	for (size_t i = 0; i < _particlePointers.size(); i++) {
		//_handlerParticles->push_back(HandlerMoleculeType((*_particles)[i]));
		_handlerParticles->push_back(HandlerMoleculeType (
						_particlePointers[i]->id(),
						_particlePointers[i]->componentid(),
						_particlePointers[i]->r(0),
						_particlePointers[i]->r(1),
						_particlePointers[i]->r(2),
						_particlePointers[i]->v(0),
						_particlePointers[i]->v(1),
						_particlePointers[i]->v(2),
						_particlePointers[i]->q().qw(),
						_particlePointers[i]->q().qx(),
						_particlePointers[i]->q().qy(),
						_particlePointers[i]->q().qz(),
						_particlePointers[i]->D(0),
						_particlePointers[i]->D(1),
						_particlePointers[i]->D(2),
						&components)
				);

#ifdef VECTORIZE
		_moleculePositions[0][i] = _particlePointers[i]->r(0);
		_moleculePositions[1][i] = _particlePointers[i]->r(1);
		_moleculePositions[2][i] = _particlePointers[i]->r(2);
#endif
	}

	for (size_t i = 0; i < _particlePointers.size(); i++) {
		(*_handlerParticles)[i].upd_cache();
	}
}


template <class Molecule, class HandlerMolecule>
void HybridCell::convertToMoleculeType() {

	#ifndef NDEBUG
	assert(_currentMoleculeType == HandlerMolecule);
	_currentMoleculeType = BasicMolecule;
	#endif

	for (size_t i = 0; i < _particlePointers.size(); i++) {
		assert(_particlePointers[i]->id() == (*_handlerParticles)[i].id());
		assert(_particlePointers[i]->componentid() == (*_handlerParticles)[i].componentid());
		assert(_particlePointers[i]->r(0) == (*_handlerParticles)[i].r(0));
		assert(_particlePointers[i]->r(1) == (*_handlerParticles)[i].r(1));
		assert(_particlePointers[i]->r(2) == (*_handlerParticles)[i].r(2));
		assert(_particlePointers[i]->v(0) == (*_handlerParticles)[i].v(0));
		assert(_particlePointers[i]->v(1) == (*_handlerParticles)[i].v(1));
		assert(_particlePointers[i]->v(2) == (*_handlerParticles)[i].v(2));


		(*_handlerParticles)[i].calcFM();

		double force[3];
		force[0] = (*_handlerParticles)[i].F(0);
		force[1] = (*_handlerParticles)[i].F(1);
		force[2] = (*_handlerParticles)[i].F(2);
		_particlePointers[i]->setF(force);
		double M[3];
		M[0] = (*_handlerParticles)[i].M(0);
		M[1] = (*_handlerParticles)[i].M(1);
		M[2] = (*_handlerParticles)[i].M(2);
		_particlePointers[i]->setM(M);
	}

	memoryManager.releaseMoleculeArray(_handlerParticles);
#ifdef VECTORIZE
	memoryManager.releaseFPMemory(_moleculePositions[0]);
	memoryManager.releaseFPMemory(_moleculePositions[1]);
	memoryManager.releaseFPMemory(_moleculePositions[2]);
#endif
	_handlerParticles = NULL;
}

#endif /*CELL_H_*/
