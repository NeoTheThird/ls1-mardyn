#include <fstream>
#include <ctime>

#include "DomainDecompDummy.h"

#include "molecules/MoleculeTypes.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/DynamicArray.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

DomainDecompDummy::DomainDecompDummy() {
}

DomainDecompDummy::~DomainDecompDummy() {
}

void DomainDecompDummy::exchangeMolecules(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {

	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double halo_L[3]; // width of the halo strip //ENABLE_MPI
	for (int i = 0; i < 3; i++) {
		rmin[i] = moleculeContainer->getBoundingBoxMin(i);
		rmax[i] = moleculeContainer->getBoundingBoxMax(i);
		halo_L[i] = moleculeContainer->get_halo_L(i);
	}

	Molecule* currentMolecule;
	// molecules that have to be copied (because of halo), get a new position
	double new_position[3];

	double phaseSpaceSize[3];

	double low_limit; // particles below this limit have to be copied or moved to the lower process
	double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process

	// 8192 should be ok as a compromise: should use 1.2 MB for BasicMolecules
	// in Double Precision and about 4MB for CachingMolecules
	utils::DynamicArray<Molecule, false, false> moleculesToAdd(8192, 1.5);

	for (unsigned short d = 0; d < 3; ++d) {
		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit = rmin[d] + halo_L[d];
		high_limit = rmax[d] - halo_L[d];
		currentMolecule = moleculeContainer->begin();

		//cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << endl;
		//cout << "halo_L: " << halo_L[0] << " / " << halo_L[1] << " / " << halo_L[2] << endl;
		//cout << "proc_domain_L: " << proc_domain_L[0] << " / " << proc_domain_L[1] << " / " << proc_domain_L[2] << endl;
		while (currentMolecule != moleculeContainer->end()) {
			const double& rd = currentMolecule->r(d);
			if (rd < low_limit) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						new_position[d2] = rd + phaseSpaceSize[d2];
					else
						new_position[d2] = currentMolecule->r(d2);
				}
				Molecule m1 = Molecule(currentMolecule->id(),currentMolecule->componentid(),
				                       new_position[0], new_position[1], new_position[2],
				                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
				                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
				                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2), &components);

				//moleculeContainer->addParticle(m1);
				moleculesToAdd.push_back(m1);

				currentMolecule = moleculeContainer->next();
			}
			else if (rd >= high_limit) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						new_position[d2] = rd - phaseSpaceSize[d2];
					else
						new_position[d2] = currentMolecule->r(d2);
				}
				Molecule m1 = Molecule(currentMolecule->id(),currentMolecule->componentid(),
				                       new_position[0], new_position[1], new_position[2],
				                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
				                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
				                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2), &components);

				//moleculeContainer->addParticle(m1);
				moleculesToAdd.push_back(m1);

				currentMolecule = moleculeContainer->next();
			}
			else
				currentMolecule = moleculeContainer->next();
		}
	}

	for (size_t i =  0; i < moleculesToAdd.size(); i++) {
		moleculeContainer->addParticle(moleculesToAdd[i]);
	}
	global_log->info() << "DomainDecompDummy::exchangeMolecules() Memory of tmpArray: " << (moleculesToAdd.capacity() * sizeof(Molecule) / 1024) << " kB "
			<< " capacity=" << moleculesToAdd.capacity() << " size=" << moleculesToAdd.size() << endl;
}

void DomainDecompDummy::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {
	exchangeMolecules(moleculeContainer, components, domain);
}

unsigned long DomainDecompDummy::countMolecules(ParticleContainer* moleculeContainer, vector<unsigned long> &compCount) {
	for (unsigned i = 0; i < compCount.size(); i++) {
		compCount[i] = 0;
	}
	Molecule* tempMolecule;
	for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()) {
		compCount[tempMolecule->componentid()] += 1;
	}
	int numMolecules = 0;
	for (unsigned i = 0; i < compCount.size(); i++) {
		numMolecules += compCount[i];
	}
	return numMolecules;
}

double DomainDecompDummy::getBoundingBoxMin(int dimension, Domain* domain) {
	return 0.0;
}
double DomainDecompDummy::getBoundingBoxMax(int dimension, Domain* domain) {
	return domain->getGlobalLength(dimension);
}

void DomainDecompDummy::writeMoleculesToFile(string filename, ParticleContainer* moleculeContainer) {

	ofstream checkpointfilestream(filename.c_str(), ios::app);

	Molecule* tempMolecule;
	for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()) {
		tempMolecule->write(checkpointfilestream);
	}
	checkpointfilestream.close();
}

double DomainDecompDummy::getTime() {
	return double(clock()) / CLOCKS_PER_SEC;
}

unsigned DomainDecompDummy::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	*minrnd = 0.0;
	*maxrnd = 1.0;
	return localN;
}

void DomainDecompDummy::assertIntIdentity(int IX) {
}

void DomainDecompDummy::assertDisjunctivity(TMoleculeContainer* mm) {
}

void DomainDecompDummy::printDecomp(std::string filename, Domain* domain) {
	global_log->warning() << "printDecomp useless in serial mode" << endl;
}
