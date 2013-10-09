#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"

using namespace std;

ParticleCell::ParticleCell() :
		_cellDataSoA(0) {
}

ParticleCell::~ParticleCell() {
	assert(!_cellDataSoA);
}

void ParticleCell::removeAllParticles() {
	molecules.clear();
}

void ParticleCell::addParticle(Molecule* particle_ptr) {
	molecules.push_back(particle_ptr);
}


void ParticleCell::update() {
	removeAllParticles();

	for (int i = 0; i < NBUCKET*NBUCKET; i++) {
		std::vector<Molecule*>& currentBucket = _buckets[i];
		int nmols = currentBucket.size();
		int kmax = nmols;
		bool swapped;
		do {
			swapped = false;
			for (int k = 0; k < kmax-1; k++) {
				if (currentBucket[k+1]->r(2) < currentBucket[k]->r(2)){
					Molecule* tmp = currentBucket[k];
					currentBucket[k] = currentBucket[k+1];
					currentBucket[k+1] = tmp;
					swapped = true;
				}
			}
			kmax = kmax-1;
		} while (swapped == true);

		molecules.insert(molecules.end(), currentBucket.begin(), currentBucket.end());
		currentBucket.clear();
	}
}


vector<Molecule*>& ParticleCell::getParticlePointers() {
	return molecules;
}

int ParticleCell::getMoleculeCount() const {
	return molecules.size();
}

bool ParticleCell::deleteMolecule(unsigned long molecule_id) {
	bool found = false;
	vector<Molecule*>::iterator molecule_iter;

	for (molecule_iter = molecules.begin(); molecule_iter != molecules.end(); molecule_iter++) {
		Molecule *molecule = *molecule_iter;
		if (molecule->id() == molecule_id) {
			found = true;
			molecules.erase(molecule_iter);
			break;
		}
	}
	return found;
}
