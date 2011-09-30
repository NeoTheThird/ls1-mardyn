#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "molecules/MoleculeTypes.h"

#include <mpi.h>
#include "utils/Logger.h"

#define DIM 3

//! @brief class to represent that particle data that is necessary for the exchange between processes
//! @author Martin Buchholz
//!
//! The data that has to be transfered for each particle consists of several values of different
//! types. As MPI-commands are expensive, it is desirable to transfer all data with a single
//! command. Therefore, two things have to be done:
//! - The data has to be "packed together", so all data at succeeding memory locations
//! - A MPI-datatype has to be defined.
//! This class achieves both. It is a class without constructor and destructor and only static
//! methods, which means that only the member variables use memory, so it can be used like a
//! C stuct. The static method setMPIType sets the given type to represent all the member variables.
class ParticleData {
public:

	static void setMPITypeForBasicMolecule(MPI_Datatype &sendPartType) {
		int blocklengths[4] = { 1, 1, 13, 1 }; // 1 unsLong value (id), 1 int value (cid), 13 double values (3r, 3v, 4q, 3D)
		MPI_Datatype types[3];
		if ( IsSame<BasicMolecule::fp_type,float>::Result::value ) {
			types[0] = MPI_UNSIGNED_LONG;
			types[1] = MPI_INT;
			types[2] = MPI_FLOAT;
			types[3] = MPI_UB;
		} else {
			types[0] = MPI_UNSIGNED_LONG;
			types[1] = MPI_INT;
			types[2] = MPI_DOUBLE;
			types[3] = MPI_UB;
		}

		MPI_Aint base;
		MPI_Aint displacements[4];
		BasicMolecule dummyMolecules[2];
		MPI_CHECK( MPI_Address(&dummyMolecules[0], &base) );
		MPI_CHECK( MPI_Address(&dummyMolecules[0]._id, displacements) );
		MPI_CHECK( MPI_Address(&dummyMolecules[0]._componentid, displacements + 1) );
		MPI_CHECK( MPI_Address(&dummyMolecules[0]._r, displacements + 2) ); // Note: it is probably potentially dangerous to
		                                                                 // assume that the quaternion just takes 4 double values
		                                                                 // (the compiler might add some vtable pointers, alignment, etc...)
		MPI_CHECK( MPI_Address(&dummyMolecules[1], displacements + 3) );

		for (int i = 0; i < 4; i++) {
			displacements[i] -= base;
		}

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
		MPI_CHECK( MPI_Type_create_struct(4, blocklengths, displacements, types, &sendPartType) );
#else
		MPI_CHECK( MPI_Type_struct(3, blocklengths, displacements, types, &sendPartType) );
#endif
		MPI_CHECK( MPI_Type_commit(&sendPartType) );
	}


	//! @brief defines a MPI datatype which can be used to transfer a MacroscopicData object
	static void setMPIType(MPI_Datatype &sendPartType) {
		int blocklengths[] = { 1, 1, 13 }; // 1 unsLong value (id), 1 int value (cid), 13 double values (3r, 3v, 4q, 3D)
		MPI_Datatype types[] = { MPI_UNSIGNED_LONG, MPI_INT, MPI_DOUBLE };

		MPI_Aint displacements[3];
		ParticleData pdata_dummy;
		MPI_CHECK( MPI_Address(&pdata_dummy, displacements) );
		MPI_CHECK( MPI_Address(&pdata_dummy.cid, displacements + 1) );
		MPI_CHECK( MPI_Address(&pdata_dummy.r[0], displacements + 2) );
		MPI_Aint base = displacements[0];
		for (int i = 0; i < 3; i++)
			displacements[i] -= base;

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
		MPI_CHECK( MPI_Type_create_struct(3, blocklengths, displacements, types, &sendPartType) );
#else
	MPI_CHECK( MPI_Type_struct(3, blocklengths, displacements, types, &sendPartType) );
#endif
		MPI_CHECK( MPI_Type_commit(&sendPartType) );
	}

	//! @brief copy data from object of class Molecule to object of class ParticleData
	static void MoleculeToParticleData(ParticleData &particleStruct, Molecule &molecule) {
		particleStruct.id = molecule.id();
		particleStruct.cid = molecule.componentid();
		particleStruct.r[0] = molecule.r(0);
		particleStruct.r[1] = molecule.r(1);
		particleStruct.r[2] = molecule.r(2);
		particleStruct.v[0] = molecule.v(0);
		particleStruct.v[1] = molecule.v(1);
		particleStruct.v[2] = molecule.v(2);
		particleStruct.q[0] = molecule.q().qw();
		particleStruct.q[1] = molecule.q().qx();
		particleStruct.q[2] = molecule.q().qy();
		particleStruct.q[3] = molecule.q().qz();
		particleStruct.D[0] = molecule.D(0);
		particleStruct.D[1] = molecule.D(1);
		particleStruct.D[2] = molecule.D(2);
	}

	//! @brief copy data from object of class class ParticleData to object of class Molecule
	static void ParticleDataToMolecule(ParticleData &particleStruct, Molecule **molecule, const std::vector<Component>* components = NULL) {
		*molecule = new Molecule(particleStruct.id, particleStruct.cid,
		                         particleStruct.r[0], particleStruct.r[1], particleStruct.r[2],
		                         particleStruct.v[0], particleStruct.v[1], particleStruct.v[2],
		                         particleStruct.q[0], particleStruct.q[1], particleStruct.q[2], particleStruct.q[3],
		                         particleStruct.D[0], particleStruct.D[1], particleStruct.D[2],
		                         components
		);
	}

	unsigned long id;
	int cid;
	double r[3];
	double v[3];
	double q[4];
	double D[3];
};

#endif /*PARTICLE_H_*/
