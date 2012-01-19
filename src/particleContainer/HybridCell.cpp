#include "particleContainer/HybridCell.h"
#include "molecules/MoleculeTypes.h"

using namespace std;

MemoryManager HybridCell::memoryManager;

HybridCell::HybridCell() : _handlerParticles(NULL)
{
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
#ifndef NDEBUG
	_currentMoleculeType = BasicMolecule;
#endif
	std::cout << "Created HybridCell." << endl;
}


HybridCell::HybridCell(const HybridCell& other) {
	this->_haloCellState = false;
	this->_boundaryCellState = false;
	this->_innerCellState = false;
#ifndef NDEBUG
	_currentMoleculeType = BasicMolecule;
#endif
	//_particles = new MoleculeArray(*(other._particles));
	_handlerParticles = NULL;
	//std::cout << "CopyConstructor BlockedCell: _particles at " << _particles << std::endl;
}


HybridCell::~HybridCell() {
	assert(_currentMoleculeType == BasicMolecule);
	assert(_handlerParticles == NULL);
}

void HybridCell::removeAllParticles() {
	assert(_currentMoleculeType == BasicMolecule);
	this->_particlePointers.clear();
}

void HybridCell::addParticle(Molecule* particle) {
	assert(_currentMoleculeType == BasicMolecule);
	_particlePointers.push_back(particle);
}

vector<Molecule*>& HybridCell::getParticlePointers() {
	return this->_particlePointers;
}

HandlerMoleculeTypeArray& HybridCell::getHandlerTypeParticles() {
	assert(_currentMoleculeType == HandlerMolecule);
	return *(this->_handlerParticles);
}

void HybridCell::assingCellToHaloRegion() {
	this->_haloCellState = true;
}

void HybridCell::assignCellToBoundaryRegion() {
	this->_boundaryCellState = true;
}

void HybridCell::assignCellToInnerRegion() {
	this->_innerCellState = true;
}

bool HybridCell::isHaloCell() const {
	return _haloCellState;
}

bool HybridCell::isBoundaryCell() const {
	return _boundaryCellState;
}

bool HybridCell::isInnerCell() const {
	return _innerCellState;
}

int HybridCell::getMoleculeCount() const {
	return _particlePointers.size();
}

bool HybridCell::deleteMolecule(unsigned long molid) {
	assert(_currentMoleculeType == BasicMolecule);
	bool found = false;
	vector<Molecule*>::iterator cellit;

	for (cellit = _particlePointers.begin(); cellit != _particlePointers.end(); cellit++) {
		if ((*cellit)->id() == molid) {
			found = true;
			_particlePointers.erase(cellit);
			break;
		}
	}
	return found;
}

