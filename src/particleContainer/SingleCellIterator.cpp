/*
 * SingleCellIterator.cpp
 *
 *  Created on: 27 Sep 2017
 *      Author: tchipevn
 */

#include "SingleCellIterator.h"
#include "ParticleCellBase.h"

SingleCellIterator::SingleCellIterator(ParticleCellBase * cell_arg,
		size_t index_arg) : _cell(cell_arg), _mol_index(index_arg), _currentParticleDeleted(false) {
}

Molecule& SingleCellIterator:: operator * () const {
	// .at method performs automatically an out-of-bounds check
	Molecule *moleculePtr = nullptr;

#ifdef ENABLE_REDUCED_MEMORY_MODE
	moleculePtr = const_cast<Molecule *>(& _AoSMoleculeReservoir);
#endif

	_cell->moleculesAtNew(_mol_index, moleculePtr);

	return *moleculePtr;
}

void SingleCellIterator :: operator ++ () {
	// if the "current" particle was deleted, then there is a new particle at _mol_index
	// and the _mol_index value should not be incremented.
	_mol_index += _currentParticleDeleted ? 0 : 1;

	_currentParticleDeleted = false;
}

void SingleCellIterator :: deleteCurrentParticle () {
	_cell->deleteMoleculeByIndex(_mol_index);
	_currentParticleDeleted = true;
}

bool SingleCellIterator::isValid() const {
	return _cell != nullptr and _mol_index < static_cast<size_t>(_cell->getMoleculeCount());
}
