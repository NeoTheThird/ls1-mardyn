#ifndef CELL_H_
#define CELL_H_

#include "utils/DynamicArray.h"
#include "molecules/MoleculeTypes.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"

#include "MemoryManager.h"
#include "Simulation.h"
#include "Domain.h"

//! @brief Cell data structure.
//! @author Martin Buchholz
//!
//! A Cell represents a small cuboid area of the domain and stores a list of 
//! pointers to the molecules in that area. Depending on the actual position
//! of the cell, it belongs to one of four different regions: \n
//! - completele outside (more than the cutoffradius away from all cells that
//!                       belong directly to the MoleculeContainer)
//!                       Such a cell shouldn't exist!!!
//! - halo region (not belonging directly to the MoleculeContainer, but within
//!                       the cutoffradius of at least one cell of the MoleculeContainer)
//! - boundary region (belonging directly to the MoleculeContainer, but not more than
//!                           the cutoffradius away from the boundary)
//! - inner region (more than the cutoffradius away from the boundary)
//! 
//! There are three boolean member variables for the last three regions. \n
//! If more than one of them is true, there must be an error in the code \n
//! If none of them is true, the cell wasn't assigned to any region yet.
//! A cell which is completely outside shouldn't exist, as it completely useless.
//! 
//! For cuboid domains, a cell doesn't need much information
//! (you could even avoid a seperate class completely)
//! Most information (position, neighbours,...) can be calculated
//! from the cell's index in the linked cell vector.
//! e.g. for a domain of 10x10x10 cells, the cell with index (counting from zero)
//! 561 has spacial index (counting from zero)  (1,6,5) and the cell indices of
//! the neighbours are e.g. 450, 451, 452, 460, 461, ... 
//! But for either unregular domains or unregular spacial decompositions
//! in the parallelisation (e.g. space filling curves), it's not easy to store that
//! information in the linked cell vector (if it is a vector at all)
//! That's why this class was introduced, so additional information which is
//! specific to a cell should be stored here.

class BlockedCell {


private:

	static MemoryManager memoryManager;

	//! the basic molecules contained in this cell
	utils::DynamicArray<Molecule, true, false>* _particles;

	//! the molecules as the handler expects them
	utils::DynamicArray<HandlerMoleculeType, true, false>* _handlerParticles;

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
	BlockedCell();

	BlockedCell(const BlockedCell& other);

	~BlockedCell();

	//! removes all elements from the list molecules
	void removeAllParticles();

	//! insert a single molecule into this cell
	void addParticle(const Molecule& particle);

	//! return a reference to the list of molecules in this cell
	utils::DynamicArray<Molecule, true, false>& getParticles();

	//! return a reference to the list of molecules in this cell
	//! in the representation for the ParticlePairsHandler
	utils::DynamicArray<HandlerMoleculeType, true, false>& getHandlerTypeParticles();

	//! @todo Return type bool neccessary!?
	bool deleteMolecule(unsigned long molid);

	//! delete molecule at iterator position.
	//! the iterator is set to point to the next particle
	utils::DynamicArray<Molecule, true, false>::iterator& deleteMolecule(utils::DynamicArray<Molecule, true, false>::iterator& it);

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

};



template <class Molecule, class HandlerMoleculeType>
void BlockedCell::convertToHandlerMoleculeType() {

	#ifndef NDEBUG
	assert(_currentMoleculeType == BasicMolecule);
	_currentMoleculeType = HandlerMolecule;
	#endif

	std::vector<Component>& components = global_simulation->getDomain()->getComponents();

	_handlerParticles = memoryManager.getMoleculeArray();

	for (size_t i = 0; i < _particles->size(); i++) {
		//_handlerParticles->push_back(HandlerMoleculeType((*_particles)[i]));
		_handlerParticles->push_back(HandlerMoleculeType (
						(*_particles)[i].id(),
						(*_particles)[i].componentid(),
						(*_particles)[i].r(0),
						(*_particles)[i].r(1),
						(*_particles)[i].r(2),
						(*_particles)[i].v(0),
						(*_particles)[i].v(1),
						(*_particles)[i].v(2),
						(*_particles)[i].q().qw(),
						(*_particles)[i].q().qx(),
						(*_particles)[i].q().qy(),
						(*_particles)[i].q().qz(),
						(*_particles)[i].D(0),
						(*_particles)[i].D(1),
						(*_particles)[i].D(2),
						&components)
				);
	}

	for (size_t i = 0; i < _particles->size(); i++) {
		(*_handlerParticles)[i].upd_cache();
	}
}


template <class Molecule, class HandlerMolecule>
void BlockedCell::convertToMoleculeType() {

	#ifndef NDEBUG
	assert(_currentMoleculeType == HandlerMolecule);
	_currentMoleculeType = BasicMolecule;
	#endif

	for (size_t i = 0; i < _particles->size(); i++) {
		assert((*_particles)[i].id() == (*_handlerParticles)[i].id());
		assert((*_particles)[i].componentid() == (*_handlerParticles)[i].componentid());
		assert((*_particles)[i].r(0) == (*_handlerParticles)[i].r(0));
		assert((*_particles)[i].r(1) == (*_handlerParticles)[i].r(1));
		assert((*_particles)[i].r(2) == (*_handlerParticles)[i].r(2));
		assert((*_particles)[i].v(0) == (*_handlerParticles)[i].v(0));
		assert((*_particles)[i].v(1) == (*_handlerParticles)[i].v(1));
		assert((*_particles)[i].v(2) == (*_handlerParticles)[i].v(2));


		(*_handlerParticles)[i].calcFM();

		double force[3];
		force[0] = (*_handlerParticles)[i].F(0);
		force[1] = (*_handlerParticles)[i].F(1);
		force[2] = (*_handlerParticles)[i].F(2);
		(*_particles)[i].setF(force);
		double M[3];
		M[0] = (*_handlerParticles)[i].M(0);
		M[1] = (*_handlerParticles)[i].M(1);
		M[2] = (*_handlerParticles)[i].M(2);
		(*_particles)[i].setM(M);
	}

	memoryManager.releaseMoleculeArray(_handlerParticles);

}

#endif /*CELL_H_*/
