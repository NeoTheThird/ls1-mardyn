#include "particleContainer/BlockedCell.h"
#include "molecules/Molecule.h"

using namespace std;

BlockedCell::BlockedCell() :	_particles(NULL), _handlerParticles(NULL)
{
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
#ifndef NDEBUG
	_currentMoleculeType = BasicMolecule;
#endif
	_particles = new utils::DynamicArray<Molecule, true, false>(2);
	_handlerParticles = _particles;
	//std::cout << "Constructor BlockedCell: _particles at " << _particles << std::endl;
}


BlockedCell::BlockedCell(const BlockedCell& other) {
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
#ifndef NDEBUG
	_currentMoleculeType = BasicMolecule;
#endif
	_particles = new utils::DynamicArray<Molecule, true, false>(*(other._particles));
	_handlerParticles = _particles;
	//std::cout << "CopyConstructor BlockedCell: _particles at " << _particles << std::endl;
}


BlockedCell::~BlockedCell() {
	delete _particles;
}

void BlockedCell::removeAllParticles() {
	assert(_currentMoleculeType == BasicMolecule);
	_particles->clear();
}

void BlockedCell::addParticle(const Molecule& particle) {
	assert(_currentMoleculeType == BasicMolecule);
	_particles->push_back(particle);
}

utils::DynamicArray<Molecule, true, false>& BlockedCell::getParticles() {
	assert(_currentMoleculeType == BasicMolecule);
	return *(this->_particles);
}

utils::DynamicArray<HandlerMoleculeType, true, false>& BlockedCell::getHandlerTypeParticles() {
	assert(_currentMoleculeType == HandlerMolecule);
	return *(this->_handlerParticles);
}

void BlockedCell::assingCellToHaloRegion() {
	this->_haloCellState = true;
}

void BlockedCell::assignCellToBoundaryRegion() {
	this->_boundaryCellState = true;
}

void BlockedCell::assignCellToInnerRegion() {
	this->_innerCellState = true;
}

bool BlockedCell::isHaloCell() const {
	return _haloCellState;
}

bool BlockedCell::isBoundaryCell() const {
	return _boundaryCellState;
}

bool BlockedCell::isInnerCell() const {
	return _innerCellState;
}

int BlockedCell::getMoleculeCount() const {
	return _particles->size();
}

bool BlockedCell::deleteMolecule(unsigned long molid) {
	assert(_currentMoleculeType == BasicMolecule);
	bool found = false;
	utils::DynamicArray<Molecule, true, false>::iterator cellit = _particles->begin();

	while (cellit != _particles->end()) {
		if ((*cellit).id() == molid) {
			found = true;
			_particles->erase(cellit);
			break;
		}
		++cellit;
	}
	return found;
}


utils::DynamicArray<Molecule, true, false>::iterator& BlockedCell::deleteMolecule(utils::DynamicArray<Molecule, true, false>::iterator& it) {
	assert(_currentMoleculeType == BasicMolecule);
	return _particles->erase(it);
}
