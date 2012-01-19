#include "DomainDecomposition.h"

#include "molecules/MoleculeTypes.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/BlockedReorderedLinkedCells.h"
#include "Domain.h"
#include "parallel/ParticleData.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;
using namespace utils;

#define SEND_RAWBYTES

DomainDecomposition::DomainDecomposition() {

	int period[DIM]; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder; // 1(true) if the ranking may be reordered by MPI_Cart_create
	int num_procs; // Number of processes

	// We create a torus topology, so all boundary conditions are periodic
	for (int d = 0; d < DIM; d++)
		period[d] = 1;
	// Allow reordering of process ranks
	reorder = 1;
	// Find out appropriate grid dimensions
	MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &num_procs) );
	setGridSize(num_procs);
	// Create the communicator
	MPI_CHECK( MPI_Cart_create(MPI_COMM_WORLD, DIM, _gridSize, period, reorder, &_comm) );

	// introduce coordinates
	MPI_CHECK( MPI_Comm_rank(_comm, &_rank) );
	MPI_CHECK( MPI_Cart_coords(_comm, _rank, DIM, _coords) );
	// find lower and higher neighbours:
	for (int d = 0; d < DIM; d++) {
		MPI_CHECK( MPI_Cart_shift(_comm, d, 1, &_neighbours[d][LOWER], &_neighbours[d][HIGHER]) );
	}
	// Initialize MPI Dataype for the particle exchange once at the beginning.
	ParticleData::setMPIType(_mpi_Particle_data);
	ParticleData::setMPITypeForBasicMolecule(_mpi_basic_molecule_data);
}

DomainDecomposition::~DomainDecomposition() {
}

void DomainDecomposition::exchangeMolecules(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {
	BlockedReorderedLinkedCells* blockedReorderedLinkedCells = dynamic_cast<BlockedReorderedLinkedCells*>(moleculeContainer);
	if (IsSame<Molecule, BasicMolecule>::Result::value && blockedReorderedLinkedCells != NULL) {
		exchangeBasicMolecules(blockedReorderedLinkedCells, components, domain);
	} else {
		exchangeParticleData(moleculeContainer, components, domain);
	}
}


void DomainDecomposition::exchangeBasicMolecules(BlockedReorderedLinkedCells* moleculeContainer, const vector<Component>& components, Domain* domain) {

	global_log->info() << "Performing Exchange of BasicMolecules!" << endl;

	// corners of the process-specific domain
	double rmin[DIM]; // lower corner
	double rmax[DIM]; // higher corner
	double halo_L[DIM]; // width of the halo strip

	for (int d = 0; d < DIM; d++) {
		rmin[d] = moleculeContainer->getBoundingBoxMin(d);
		rmax[d] = moleculeContainer->getBoundingBoxMax(d);
		halo_L[d] = moleculeContainer->get_halo_L(d);
	}

#ifndef NDEBUG
	vector<unsigned long> compCount;
	compCount.resize(1);
	countMolecules(moleculeContainer, compCount);
	global_log->debug() << "Num molecules: " << compCount[0] << endl;
#endif

	// temporal data for the particle exchange
	unsigned int numPartsToSend[DIM][2];
	unsigned int numPartsToRecv[DIM][2];
	DynamicArray<BasicMolecule, true, false>* particlesSendBuffs[DIM][2];
	BasicMolecule* particlesRecvBuffs[DIM][2];

	// MPI communication status and requests
	MPI_Status status;
	MPI_Status send_statuses[DIM][2];
	MPI_Status recv_statuses[DIM][2];
	MPI_Request send_requests[DIM][2];
	MPI_Request recv_requests[DIM][2];

	int direction;

	for (unsigned short d = 0; d < DIM; d++) {
		// when moving a particle across a periodic boundary, the molecule position has to change
		// these offset specify for each dimension (x, y and z) and each direction ("left"/lower
		// neighbour and "right"/higher neighbour, how the paritcle coordinates have to be changed.
		// e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
		// moving to the left get the length of the whole domain added to their x-value
		double offsetLower[DIM];
		double offsetHigher[DIM];
		offsetLower[d] = 0.0;
		offsetHigher[d] = 0.0;

		// process on the left boundary
		if (_coords[d] == 0)
			offsetLower[d] = domain->getGlobalLength(d);
		// process on the right boundary
		if (_coords[d] == _gridSize[d] - 1)
			offsetHigher[d] = -domain->getGlobalLength(d);

		double regToSendLow[DIM]; // Region that belongs to a neighbouring process
		double regToSendHigh[DIM]; // -> regToSendLow
		for (direction = LOWER; direction <= HIGHER; direction++) {
			// find the region that each neighbour will get
			for (int i = 0; i < DIM; i++) {
				regToSendLow[i] = rmin[i] - halo_L[i];
				regToSendHigh[i] = rmax[i] + halo_L[i];
			}
			switch (direction) {
			case LOWER:
				regToSendHigh[d] = rmin[d] + halo_L[d];
				break;
			case HIGHER:
				regToSendLow[d] = rmax[d] - halo_L[d];
				break;
			}

			particlesSendBuffs[d][direction] = new DynamicArray<BasicMolecule, true, false>;
			moleculeContainer->getRegion(regToSendLow, regToSendHigh, *particlesSendBuffs[d][direction]);

			// initialize send buffer
			numPartsToSend[d][direction] = particlesSendBuffs[d][direction]->size();
			//particlesSendBuffs[d][direction] = new ParticleData[numPartsToSend[d][direction]];

			//std::list<Molecule*>::iterator particlePtrIter;
			unsigned int partCount = 0;
			double shift = 0.0;
			if (direction == LOWER)
				shift = offsetLower[d];
			if (direction == HIGHER)
				shift = offsetHigher[d];

			for (size_t particleIter = 0; particleIter < numPartsToSend[d][direction]; particleIter++) {
			//for (particlePtrIter = particlePtrsToSend.begin(); particlePtrIter != particlePtrsToSend.end(); particlePtrIter++) {
				// copy relevant data from the Molecule to ParticleData type
				//ParticleData::MoleculeToParticleData(particlesSendBuffs[d][direction][partCount], **particlePtrIter);
				// add offsets for particles transfered over the periodic boundary
				DynamicArray<BasicMolecule, true, false>* arrayPtr = particlesSendBuffs[d][direction];
				arrayPtr->operator [](partCount).setr( d, arrayPtr->operator [](partCount).r(d) + shift);
				//std::cout << "[At " << _rank << "] Sending Particle index=" << particlesSendBuffs[d][direction]->operator [](partCount).id() << " to [" << _neighbours[d][direction] << std::endl;
				partCount++;
			}
		}

		// Communicate to lower and higher neighbour
		for (direction = LOWER; direction <= HIGHER; direction++) {
			// Send to lower, receive from upper
			// Send number of values that have to be sent

#ifdef SEND_RAWBYTES
			int numsend = numPartsToSend[d][direction] * sizeof(BasicMolecule);
			int numrecv;
			// Send values to lower/upper and receive values from upper/lower
			MPI_CHECK( MPI_Isend(&((*particlesSendBuffs[d][direction])[0]), numsend, MPI_BYTE, _neighbours[d][direction], 99, _comm, &send_requests[d][direction]) );
			MPI_CHECK( MPI_Probe(_neighbours[d][(direction + 1) % 2], 99, _comm, &status) );
			MPI_CHECK( MPI_Get_count(&status, MPI_BYTE, &numrecv) );
			assert(numrecv % sizeof(BasicMolecule) == 0);
			numrecv = numrecv / sizeof(BasicMolecule);
			// initialize receive buffer
			particlesRecvBuffs[d][direction] = new BasicMolecule[numrecv];
			numPartsToRecv[d][direction] = numrecv;
			MPI_CHECK( MPI_Irecv(particlesRecvBuffs[d][direction], numrecv * sizeof (BasicMolecule), MPI_BYTE, _neighbours[d][(direction + 1) % 2], 99, _comm, &recv_requests[d][direction]) );
#else
			int numsend = numPartsToSend[d][direction];
			int numrecv;
			// Send values to lower/upper and receive values from upper/lower
			MPI_CHECK( MPI_Isend(&((*particlesSendBuffs[d][direction])[0]), numsend, _mpi_basic_molecule_data, _neighbours[d][direction], 99, _comm, &send_requests[d][direction]) );
			MPI_CHECK( MPI_Probe(_neighbours[d][(direction + 1) % 2], 99, _comm, &status) );
			MPI_CHECK( MPI_Get_count(&status, _mpi_basic_molecule_data, &numrecv) );
			// initialize receive buffer
			particlesRecvBuffs[d][direction] = new BasicMolecule[numrecv];
			numPartsToRecv[d][direction] = numrecv;
			MPI_CHECK( MPI_Irecv(particlesRecvBuffs[d][direction], numrecv, _mpi_basic_molecule_data, _neighbours[d][(direction + 1) % 2], 99, _comm, &recv_requests[d][direction]) );
#endif
		}

		// warning because of the reinterpret-cast below
		if (!IsSame<Molecule, BasicMolecule>::Result::value) {
			global_log->error() << "DomainDecomposition::exchangeBasicMolecules() must only be called for BasicMolecules!" << endl;
			return;
		}

		// Insert molecules into domain
		for (direction = LOWER; direction <= HIGHER; direction++) {
			int numrecv = numPartsToRecv[d][direction];
			MPI_CHECK( MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]) );
			MPI_CHECK( MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]) );
			// insert received molecules into list of molecules
			for (int i = 0; i < numrecv; i++) {
				moleculeContainer->addParticle(reinterpret_cast<Molecule*>(particlesRecvBuffs[d][direction])[i]);
				//std::cout << "[At " << _rank << "] Receive Particle index=" << particlesSendBuffs[d][direction]->operator [](i).id() << " from [" << _neighbours[d][direction] << std::endl;
			}
			// free memory
			delete particlesSendBuffs[d][direction];
			delete[] particlesRecvBuffs[d][direction];
		}
	}
}


void DomainDecomposition::exchangeParticleData(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {

	// corners of the process-specific domain
	double rmin[DIM]; // lower corner
	double rmax[DIM]; // higher corner
	double halo_L[DIM]; // width of the halo strip

	for (int d = 0; d < DIM; d++) {
		rmin[d] = moleculeContainer->getBoundingBoxMin(d);
		rmax[d] = moleculeContainer->getBoundingBoxMax(d);
		halo_L[d] = moleculeContainer->get_halo_L(d);
	}

#ifndef NDEBUG
	vector<unsigned long> compCount;
	compCount.resize(1);
	countMolecules(moleculeContainer, compCount);
	global_log->debug() << "Num molecules: " << compCount[0] << endl;
#endif

	// temporal data for the particle exchange
	int numPartsToSend[DIM][2];
	int numPartsToRecv[DIM][2];
	ParticleData* particlesSendBuffs[DIM][2];
	ParticleData* particlesRecvBuffs[DIM][2];

	// MPI communication status and requests
	MPI_Status status;
	MPI_Status send_statuses[DIM][2];
	MPI_Status recv_statuses[DIM][2];
	MPI_Request send_requests[DIM][2];
	MPI_Request recv_requests[DIM][2];

	int direction;

	for (unsigned short d = 0; d < DIM; d++) {
		// when moving a particle across a periodic boundary, the molecule position has to change
		// these offset specify for each dimension (x, y and z) and each direction ("left"/lower
		// neighbour and "right"/higher neighbour, how the paritcle coordinates have to be changed.
		// e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
		// moving to the left get the length of the whole domain added to their x-value
		double offsetLower[DIM];
		double offsetHigher[DIM];
		offsetLower[d] = 0.0;
		offsetHigher[d] = 0.0;

		// process on the left boundary
		if (_coords[d] == 0)
			offsetLower[d] = domain->getGlobalLength(d);
		// process on the right boundary
		if (_coords[d] == _gridSize[d] - 1)
			offsetHigher[d] = -domain->getGlobalLength(d);

		double regToSendLow[DIM]; // Region that belongs to a neighbouring process
		double regToSendHigh[DIM]; // -> regToSendLow
		for (direction = LOWER; direction <= HIGHER; direction++) {
			// find the region that each neighbour will get
			for (int i = 0; i < DIM; i++) {
				regToSendLow[i] = rmin[i] - halo_L[i];
				regToSendHigh[i] = rmax[i] + halo_L[i];
			}
			switch (direction) {
			case LOWER:
				regToSendHigh[d] = rmin[d] + halo_L[d];
				break;
			case HIGHER:
				regToSendLow[d] = rmax[d] - halo_L[d];
				break;
			}

			list<Molecule*> particlePtrsToSend;
			moleculeContainer->getRegion(regToSendLow, regToSendHigh, particlePtrsToSend);

			// initialize send buffer
			numPartsToSend[d][direction] = particlePtrsToSend.size();
			particlesSendBuffs[d][direction] = new ParticleData[numPartsToSend[d][direction]];

			std::list<Molecule*>::iterator particlePtrIter;
			long partCount = 0;
			double shift = 0.0;
			if (direction == LOWER)
				shift = offsetLower[d];
			if (direction == HIGHER)
				shift = offsetHigher[d];

			for (particlePtrIter = particlePtrsToSend.begin(); particlePtrIter != particlePtrsToSend.end(); particlePtrIter++) {
				// copy relevant data from the Molecule to ParticleData type
				ParticleData::MoleculeToParticleData(particlesSendBuffs[d][direction][partCount], **particlePtrIter);
				// add offsets for particles transfered over the periodic boundary
				particlesSendBuffs[d][direction][partCount].r[d] += shift;
				partCount++;
				//std::cout << "[At " << _rank << "] Sending Particle index=" << particlesSendBuffs[d][direction][partCount].id << " to [" << _neighbours[d][direction] << std::endl;
			}
		}

		// Communicate to lower and higher neighbour
		for (direction = LOWER; direction <= HIGHER; direction++) {
			// Send to lower, receive from upper
			// Send number of values that have to be sent
			int numsend = numPartsToSend[d][direction];
			int numrecv;

			// Send values to lower/upper and receive values from upper/lower
			MPI_CHECK( MPI_Isend(particlesSendBuffs[d][direction], numsend, _mpi_Particle_data, _neighbours[d][direction], 99, _comm, &send_requests[d][direction]) );
			MPI_CHECK( MPI_Probe(_neighbours[d][(direction + 1) % 2], 99, _comm, &status) );
			MPI_CHECK( MPI_Get_count(&status, _mpi_Particle_data, &numrecv) );
			// initialize receive buffer
			particlesRecvBuffs[d][direction] = new ParticleData[numrecv];
			numPartsToRecv[d][direction] = numrecv;
			MPI_CHECK( MPI_Irecv(particlesRecvBuffs[d][direction], numrecv, _mpi_Particle_data, _neighbours[d][(direction + 1) % 2], 99, _comm, &recv_requests[d][direction]) );
		}

		// Insert molecules into domain
		for (direction = LOWER; direction <= HIGHER; direction++) {
			int numrecv = numPartsToRecv[d][direction];
			MPI_CHECK( MPI_Wait(&send_requests[d][direction], &send_statuses[d][direction]) );
			MPI_CHECK( MPI_Wait(&recv_requests[d][direction], &recv_statuses[d][direction]) );
			// insert received molecules into list of molecules
			for (int i = 0; i < numrecv; i++) {
				Molecule *m;
				ParticleData::ParticleDataToMolecule(particlesRecvBuffs[d][direction][i], &m, &components);
				moleculeContainer->addParticle(*m);
				//std::cout << "[At " << _rank << "] Received Particle index=" << particlesSendBuffs[d][direction][i].id << " from [" << _neighbours[d][direction] << std::endl;
				delete m;
			}
			// free memory
			delete[] particlesRecvBuffs[d][direction];
			delete[] particlesSendBuffs[d][direction];
		}
	}
}

void DomainDecomposition::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {
	exchangeMolecules(moleculeContainer, components, domain);
}

bool DomainDecomposition::procOwnsPos(double x, double y, double z, Domain* domain) {
	if (x < getBoundingBoxMin(0, domain) || x >= getBoundingBoxMax(0, domain))
		return false;
	else if (y < getBoundingBoxMin(1, domain) || y >= getBoundingBoxMax(1, domain))
		return false;
	else if (z < getBoundingBoxMin(2, domain) || z >= getBoundingBoxMax(2, domain))
		return false;
	else
		return true;
}

double DomainDecomposition::guaranteedDistance(double x, double y, double z, Domain* domain) {
	double xdist = 0;
	double ydist = 0;
	double zdist = 0;

	if (x < getBoundingBoxMin(0, domain))
		xdist = getBoundingBoxMin(0, domain) - x;
	else if (x >= getBoundingBoxMax(0, domain))
		xdist = x - getBoundingBoxMax(0, domain);

	if (y < getBoundingBoxMin(1, domain))
		ydist = getBoundingBoxMin(1, domain) - y;
	else if (y >= getBoundingBoxMax(1, domain))
		ydist = y - getBoundingBoxMax(1, domain);

	if (z < getBoundingBoxMin(2, domain))
		zdist = getBoundingBoxMin(2, domain) - z;
	else if (z >= getBoundingBoxMax(2, domain))
		zdist = z - getBoundingBoxMax(2, domain);

	return sqrt(xdist * xdist + ydist * ydist + zdist * zdist);
}

unsigned long DomainDecomposition::countMolecules(ParticleContainer* moleculeContainer, vector<unsigned long> &compCount) {
	vector<int> localCompCount;
	localCompCount.resize(compCount.size());
	for (unsigned int i = 0; i < localCompCount.size(); i++)
		localCompCount[i] = 0;

	Molecule* tempMolecule;
	for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()) {
		localCompCount[tempMolecule->componentid()] += 1;
	}
	int numMolecules = 0;
	for (unsigned int i = 0; i < localCompCount.size(); i++) {
		MPI_CHECK( MPI_Allreduce(&localCompCount[i], &compCount[i], 1, MPI_INT, MPI_SUM, _comm) );
		numMolecules += compCount[i];
	}
	return numMolecules;
}

double DomainDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	return _coords[dimension] * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

double DomainDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	return (_coords[dimension] + 1) * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

void DomainDecomposition::printDecomp(string filename, Domain* domain) {
	int numprocs;
	MPI_CHECK( MPI_Comm_size(_comm, &numprocs) );

	if (_rank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
		povcfgstrm << "cells " << _gridSize[0] << " " << _gridSize[1] << " " << _gridSize[2] << endl;
		povcfgstrm << "procs " << numprocs << endl;
		povcfgstrm << "data DomainDecomp" << endl;
		povcfgstrm.close();
	}

	for (int process = 0; process < numprocs; process++) {
		if (_rank == process) {
			ofstream povcfgstrm(filename.c_str(), ios::app);
			povcfgstrm << _coords[2] * _gridSize[0] * _gridSize[1] + _coords[1] * _gridSize[0] + _coords[0] << " " << _rank << endl;
			povcfgstrm.close();
		}
		barrier();
	}
}

void DomainDecomposition::writeMoleculesToFile(string filename, ParticleContainer* moleculeContainer) {

	int numprocs;
	MPI_CHECK( MPI_Comm_size(_comm, &numprocs) );

	for (int process = 0; process < numprocs; process++) {
		if (_rank == process) {
			ofstream checkpointfilestream(filename.c_str(), ios::app);
			Molecule* tempMolecule;
			for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next())
				tempMolecule->write(checkpointfilestream);
			checkpointfilestream.close();
		}
		barrier();
	}
}

inline int DomainDecomposition::getRank(int x, int y, int z) {
	int neigh_coords[DIM]; // Array for the coordinates
	int neigh_rank; // Rank of the neighbour
	neigh_coords[0] = x;
	neigh_coords[1] = y;
	neigh_coords[2] = z;
	MPI_CHECK( MPI_Cart_rank(_comm, neigh_coords, &neigh_rank) );
	return (neigh_rank);
}

int DomainDecomposition::getNumProcs() {
	int numProcs;
	MPI_CHECK( MPI_Comm_size(_comm, &numProcs) );
	return numProcs;
}

double DomainDecomposition::getTime() {
	return MPI_Wtime();
}

void DomainDecomposition::setGridSize(int num_procs) {
	for( int i = 0; i < DIM; i++ )
	_gridSize[i] = 0;
	MPI_CHECK( MPI_Dims_create( num_procs, DIM, (int *) &_gridSize ) );
}

unsigned DomainDecomposition::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	int num_procs;
	MPI_CHECK( MPI_Comm_size(_comm, &num_procs) );
	unsigned* moldistribution = new unsigned[num_procs];
	MPI_CHECK( MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, _comm) );
	unsigned globalN = 0;
	for (int r = 0; r < _rank; r++)
		globalN += moldistribution[r];
	unsigned localNbottom = globalN;
	globalN += moldistribution[_rank];
	unsigned localNtop = globalN;
	for (int r = _rank + 1; r < num_procs; r++)
		globalN += moldistribution[r];
	delete[] moldistribution;
	*minrnd = (float) localNbottom / globalN;
	*maxrnd = (float) localNtop / globalN;
	return globalN;
}

void DomainDecomposition::assertIntIdentity(int IX) {
	if (_rank)
		MPI_CHECK( MPI_Send(&IX, 1, MPI_INT, 0, 2 * _rank + 17, _comm) );
	else {
		int recv;
		int num_procs;
		MPI_CHECK( MPI_Comm_size(_comm, &num_procs) );
		MPI_Status s;
		for (int i = 1; i < num_procs; i++) {
			MPI_CHECK( MPI_Recv(&recv, 1, MPI_INT, i, 2 * i + 17, _comm, &s) );
			if (recv != IX) {
				global_log->error() << "IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
				MPI_Abort(MPI_COMM_WORLD, 911);
			}
		}
		global_log->info() << "IX = " << recv << " for all " << num_procs << " ranks.\n";
	}
}

void DomainDecomposition::assertDisjunctivity(TMoleculeContainer* mm) {
	Molecule* m;

	if (_rank) {
		int num_molecules = mm->getNumberOfParticles();
		unsigned long *tids;
		tids = new unsigned long[num_molecules];

		int i = 0;
		for (m = mm->begin(); m != mm->end(); m = mm->next()) {
			tids[i] = m->id();
			i++;
		}
		MPI_CHECK( MPI_Send(tids, num_molecules, MPI_UNSIGNED_LONG, 0, 2674 + _rank, _comm) );
		delete[] tids;
		global_log->info() << "Data consistency checked: for results see rank 0." << endl;
	}
	else {
		map<unsigned long, int> check;
		int num_procs;
		MPI_CHECK( MPI_Comm_size(_comm, &num_procs) );

		for (m = mm->begin(); m != mm->end(); m = mm->next())
			check[m->id()] = 0;

		MPI_Status status;
		for (int i = 1; i < num_procs; i++) {
			int num_recv = 0;
			unsigned long *recv;
			MPI_CHECK( MPI_Probe(i, 2674 + i, _comm, &status) );
			MPI_CHECK( MPI_Get_count(&status, MPI_UNSIGNED_LONG, &num_recv) );
			recv = new unsigned long[num_recv];

			MPI_CHECK( MPI_Recv(recv, num_recv, MPI_UNSIGNED_LONG, i, 2674 + i, _comm, &status) );
			for (int j = 0; j < num_recv; j++) {
				if (check.find(recv[j]) != check.end()) {
					global_log->error() << "Ranks " << check[recv[j]] << " and " << i << " both propagate ID " << recv << endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
					check[recv[j]] = i;
			}
			delete[] recv;
		}
		global_log->info() << "Data consistency checked: No duplicate IDs detected among " << check.size() << " entries." << endl;
	}
}
