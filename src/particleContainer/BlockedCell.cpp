#include "particleContainer/BlockedCell.h"
#include "molecules/MoleculeTypes.h"

using namespace std;

MemoryManager BlockedCell::memoryManager;

BlockedCell::BlockedCell() :	_particles(NULL), _handlerParticles(NULL)
{
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
#ifndef NDEBUG
	_currentMoleculeType = BasicMolecule;
#endif
	_particles = new MoleculeArray(2);
	_handlerParticles = NULL;
	//std::cout << "Constructor BlockedCell: _particles at " << _particles << std::endl;
}


BlockedCell::BlockedCell(const BlockedCell& other) {
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
#ifndef NDEBUG
	_currentMoleculeType = BasicMolecule;
#endif
	_particles = new MoleculeArray(*(other._particles));
	_handlerParticles = NULL;
	//std::cout << "CopyConstructor BlockedCell: _particles at " << _particles << std::endl;
}


BlockedCell::~BlockedCell() {
	assert(_currentMoleculeType == BasicMolecule);
	assert(_handlerParticles == NULL);

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

MoleculeArray& BlockedCell::getParticles() {
	assert(_currentMoleculeType == BasicMolecule);
	return *(this->_particles);
}

HandlerMoleculeTypeArray& BlockedCell::getHandlerTypeParticles() {
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
	MoleculeArray::iterator cellit = _particles->begin();

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


MoleculeArray::iterator& BlockedCell::deleteMolecule(MoleculeArray::iterator& it) {
	assert(_currentMoleculeType == BasicMolecule);
	return _particles->erase(it);
}


void BlockedCell::setToHandlerMoleculeType() {
	#ifndef NDEBUG
	assert(_currentMoleculeType == BasicMolecule);
	_currentMoleculeType = HandlerMolecule;
	std::cout << "Set ptr to HandlerMoleculeArray." << std::endl;
	#endif
	_handlerParticles = reinterpret_cast<HandlerMoleculeTypeArray*>(_particles);
}


void BlockedCell::setToMoleculeType() {
	#ifndef NDEBUG
	assert(_currentMoleculeType == HandlerMolecule);
	_currentMoleculeType = BasicMolecule;
	std::cout << "Release ptr to HandlerMoleculeArray." << std::endl;
	#endif

	_handlerParticles = NULL;
}
