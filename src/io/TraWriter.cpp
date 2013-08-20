// TraWriter.cpp

#include "io/TraWriter.h"
#include "Common.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include <list>

#ifdef ENABLE_MPI
#include <mpi.h>
#include <sys/resource.h>
#endif

using namespace std;

TraWriter::TraWriter(unsigned long writeFrequency, std::string outputPrefix, bool incremental) {
	// _filename = filename;
        _filename = outputPrefix;
        _outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	// _numberOfTimesteps = numberOfTimesteps;

	_filenameisdate = false;
}

TraWriter::~TraWriter(){}

void TraWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain){
	string filename = _filename + ".tra";
	ofstream fileout(filename.c_str(), ios::out);
	fileout.close();
	_wroteTra = false;
}

void TraWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;

}

void TraWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu){
	if (simstep % _writeFrequency == 0) {

		  struct rlimit limit;
		  /* Set the stack limit in bytes. */
		  //
		  limit.rlim_cur = 1099511627776;
		  limit.rlim_max = 1099511627776;
		  if (setrlimit(RLIMIT_STACK, &limit) != 0) {
		    cout << "setrlimit() failed.";
		    exit(1);
		  }



		stringstream filenamestream, outputstream;
		filenamestream << _filename;
		filenamestream << ".tra";
		char filename[filenamestream.str().size()+1];
		strcpy(filename,filenamestream.str().c_str());
		
		//ofstream visittfstrm(filenamestream.str().c_str());
#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		int numprocs = domainDecomp->getNumProcs();
		if (rank== 0){
#endif
			if (!_wroteTra){
				outputstream << "      ZS          x          y          z          \n";
				_wroteTra = true;
			}
			else
			{
				// outputstream << "#" << endl;
			}
#ifdef ENABLE_MPI
		}
#endif

/*
		std::list<Molecule*> moleculeList;
		double dLowerCorner[3];
		double dUpperCorner[3];

		dLowerCorner[0] = _dMidpoint[0] - _dBoxLength[0] / 2;
		dLowerCorner[1] = _dMidpoint[1] - _dBoxLength[1] / 2;
		dLowerCorner[2] = _dMidpoint[2] - _dBoxLength[2] / 2;

		dUpperCorner[0] = _dMidpoint[0] + _dBoxLength[0] / 2;
		dUpperCorner[1] = _dMidpoint[1] + _dBoxLength[1] / 2;
		dUpperCorner[2] = _dMidpoint[2] + _dBoxLength[2] / 2;

		particleContainer->getRegion(dLowerCorner, dUpperCorner, moleculeList);
*/

		outputstream << simstep;
		Molecule* pos = particleContainer->begin();


		// for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				outputstream << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(2)
				            << pos->componentid() << setprecision(3);
				for (unsigned short d = 0; d < 3; d++)
				{
					outputstream << setw(11) << pos->r(d);
				}
				outputstream << "   ";

			}
		// }
		outputstream << endl;

		long outputsize = outputstream.str().size();
		//cout << "rank: " << rank << "; step: " << simstep << "; outputsize: " << outputsize << endl;
		char output[outputsize+1];
		strcpy(output,outputstream.str().c_str());
#ifdef ENABLE_MPI
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		for (int dest = rank+1; dest < numprocs; dest++){
			int sendcount = 1;
		    int sendtag = 0;
		    MPI_Request request;
		    MPI_Isend(&outputsize, sendcount, MPI_LONG, dest, sendtag, MPI_COMM_WORLD, &request);
		}
		MPI_Status status;
		long offset = 0;
		long outputsize_get;
		for (int source = 0; source < rank; source++){
			int recvcount = 1;
		    int recvtag = 0;
		    MPI_Recv(&outputsize_get, recvcount, MPI_LONG, source, recvtag, MPI_COMM_WORLD, &status);
		    offset += outputsize_get;
		}

		//cout << "rank: " << rank << "; step: " << simstep << "; offset: " << offset << endl;
		MPI_File_seek(fh, offset, MPI_SEEK_END);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_File_write(fh, output, outputsize, MPI_CHAR, &status);
		MPI_File_close(&fh);
#else
		ofstream fileout(filename, ios::out|ios::app);
		fileout << output;
		fileout.close();
#endif
	}
}

void TraWriter::finishOutput(ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp, Domain* domain) {
}
