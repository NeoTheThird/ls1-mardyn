#ifndef COLLECTIVECOMMUNICATION_H_
#define COLLECTIVECOMMUNICATION_H_

#include "utils/Logger.h"
#include "CollectiveCommBase.h"
#include <mpi.h>
#include <cassert>

/* Enable agglomerated reduce operations. This will store all values in one array and apply a
 * user defined reduce operation so that the MPI reduce operation is only called once. */
#define ENABLE_AGGLOMERATED_REDUCE 1

//! @brief This class is used to transfer several values of different types with a single command
//! @author Nikola Tchipev, Martin Buchholz
//!
//! At different points in the simulation process, collective communication commands are
//! necessary. One thing that always has to be done is a reduce command to get global
//! macroscopic values from the local ones. But depending on the application, other commands
//! might also be possible (also broadcast commands). The number of values to be transferred
//! from/to each proc is usually very small (about ten or even less). The main costs of
//! those communications are therefore not caused by the amount of data, but by the MPI
//! overhead caused by each call of a communication command. The purpose of this class is
//! to use a single MPI command to transfer several values of possible different types.
//! Currently supported commands are:
//! - broadcast
//! - reduce using add as reduce operation
//!
//! Currently supported datatypes are:
//! - MPI_INT
//! - MPI_UNSIGNED_LONG
//! - MPI_FLOAT
//! - MPI_DOUBLE
//! - MPI_LONG_DOUBLE
//!
//! Further datatypes and reduce operations could be added very easily.
//! A typical usage of this class could look like this:
//! @code
//!   // create object
//!   CollectiveCommunication collComm;
//!   // Set number of values to be sent
//!   collComm.init(comm, 4);
//!
//!   // store values to be sent
//!   collComm.appendInt(5);
//!   collComm.appendInt(8);
//!   collComm.appendDouble(1.3);
//!   collComm.appendUnsLong(9);
//!
//!   // execute collective communication
//!   collComm.allreduceSum();
//!
//!   // read values (IMPORTANT: same order as for storing)
//!   int i1 = collComm.getInt();
//!   int i2 = collComm.getInt();
//!   double d = collComm.getDouble();
//!   unsigned long l = collComm.getUnsLong();
//!
//!   // finalize the communication (important for deleting memory)
//!   collComm.finalize();
//! @endcode
class CollectiveCommunication: public CollectiveCommBase {
public:
	CollectiveCommunication() {
		_communicator = 0;
		_agglomeratedType = MPI_DATATYPE_NULL;
	}

	virtual ~CollectiveCommunication() {
		assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}

	//! @brief allocate memory for the values to be sent, initialize counters
	//! @param numValues number of values that shall be communicated
	void init(MPI_Comm communicator, int numValues) {
		CollectiveCommBase::init(numValues);

		_communicator = communicator;
		_types.reserve(numValues);
	}

	// documentation in base class
	void finalize() {
		CollectiveCommBase::finalize();
		_types.clear();

		assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}

	// documentation in base class
	void appendInt(int intValue) {
		CollectiveCommBase::appendInt(intValue);
		_types.push_back(MPI_INT);
	}

	// documentation in base class
	void appendUnsLong(unsigned long unsLongValue) {
		CollectiveCommBase::appendUnsLong(unsLongValue);
		_types.push_back(MPI_UNSIGNED_LONG);
	}

	// documentation in base class
	void appendFloat(float floatValue) {
		CollectiveCommBase::appendFloat(floatValue);
		_types.push_back(MPI_FLOAT);
	}

	// documentation in base class
	void appendDouble(double doubleValue) {
		CollectiveCommBase::appendDouble(doubleValue);
		_types.push_back(MPI_DOUBLE);
	}

	// documentation in base class
	void appendLongDouble(long double longDoubleValue) {
		CollectiveCommBase::appendLongDouble(longDoubleValue);
		_types.push_back(MPI_LONG_DOUBLE);
	}

	//! Get the MPI communicator
	//! @return MPI communicator
	MPI_Comm getTopology() {
		return _communicator;
	}

	// Getters don't need to be overridden, see parent class

	// documentation in base class
	void broadcast(int root = 0) {
		setMPIType();
		valType * startOfValues = &(_values[0]);
		MPI_CHECK(
				MPI_Bcast(startOfValues, 1, _agglomeratedType, root,
						_communicator));
		MPI_CHECK(MPI_Type_free(&_agglomeratedType));
	}

	// documentation in base class
	void allreduceSum() {
#if ENABLE_AGGLOMERATED_REDUCE
		setMPIType();
		MPI_Op agglomeratedTypeAddOperator;
		const int commutative = 1;
		valType * startOfValues = &(_values[0]);
		MPI_CHECK(
				MPI_Op_create(
						(MPI_User_function * ) CollectiveCommunication::add,
						commutative, &agglomeratedTypeAddOperator));
		MPI_CHECK(
				MPI_Allreduce(MPI_IN_PLACE, startOfValues, 1, _agglomeratedType, agglomeratedTypeAddOperator, _communicator));
		MPI_CHECK(MPI_Op_free(&agglomeratedTypeAddOperator));
		MPI_CHECK(MPI_Type_free(&_agglomeratedType));
#else
		for( int i = 0; i < _numValues; i++ ) {
			MPI_CHECK( MPI_Allreduce( MPI_IN_PLACE, &(_values[i]), 1, _types[i], MPI_SUM, _communicator ) );
		}
#endif
	}

private:
	//! @brief defines a MPI datatype which can be used to transfer a CollectiveCommunication object
	//!
	//! before this method is called, init has to be called and all values to be
	//! communicated have to be added with the append<type> methods. Then this method
	//! will construct a MPI-datatype which represents all the added values. The
	//! datatype is stored in the member variable _valuesType;
	void setMPIType() {
		int numblocks = _values.size();
		int blocklengths[numblocks];
		MPI_Aint disps[numblocks];
		int disp = 0;
		for (int i = 0; i < numblocks; i++) {
			blocklengths[i] = 1;
			disps[i] = disp;
			disp += sizeof(valType);
		}
		MPI_Datatype * startOfTypes = &(_types[0]);
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
		MPI_CHECK(
				MPI_Type_create_struct(numblocks, blocklengths, disps,
						startOfTypes, &_agglomeratedType));
#else
		MPI_CHECK( MPI_Type_struct(numblocks, blocklengths, disps, startOfTypes, &_agglomeratedType) );
#endif
		MPI_CHECK(MPI_Type_commit(&_agglomeratedType));
	}

	//! @brief method used by MPI to add variables of this type
	//!
	//! For collective reduce commands, an operation has to be specified
	//! which should be used to combine the values from the different processes.
	//! As with this class, special datatypes are sent, the build-in
	//! MPI operations don't work on those datatypes. Therefore, operations have
	//! to be defined which work on those datatypes. MPI allows to create own
	//! operations (type MPI_Op) by specifying a function which will be used
	//! in the reduce operation. MPI specifies the signature of such functions
	//! This methods checks from which basic datatypes the given datatype
	//! was constructed and performs an add operation for each of the basic types.
	static void add(valType *invec, valType *inoutvec, int */*len*/,
			MPI_Datatype *dtype) {
		int numints;
		int numaddr;
		int numtypes;
		int combiner;

		MPI_CHECK(
				MPI_Type_get_envelope(*dtype, &numints, &numaddr, &numtypes,
						&combiner));

		int arrayInts[numints];
		MPI_Aint arrayAddr[numaddr];
		MPI_Datatype arrayTypes[numtypes];

		MPI_CHECK(
				MPI_Type_get_contents(*dtype, numints, numaddr, numtypes,
						arrayInts, arrayAddr, arrayTypes));

		for (int i = 0; i < numtypes; i++) {
			if (arrayTypes[i] == MPI_INT) {
				inoutvec[i].v_int += invec[i].v_int;
			}
			else if (arrayTypes[i] == MPI_UNSIGNED_LONG) {
				inoutvec[i].v_unsLong += invec[i].v_unsLong;
			}
			else if (arrayTypes[i] == MPI_FLOAT) {
				inoutvec[i].v_float += invec[i].v_float;
			}
			else if (arrayTypes[i] == MPI_DOUBLE) {
				inoutvec[i].v_double += invec[i].v_double;
			}
			else if (arrayTypes[i] == MPI_LONG_DOUBLE) {
				inoutvec[i].v_longDouble += invec[i].v_longDouble;
			}
		}
	}

	//! Vector of the corresponding MPI types for the values stored in _values
	std::vector<MPI_Datatype> _types;

	//! MPI_Datatype which will be used in the communication command and which
	//! represents all values
	MPI_Datatype _agglomeratedType;

	//! Communicator to be used by the communication commands
	MPI_Comm _communicator;

};

#endif /* COLLECTIVECOMMUNICATION_H_ */
