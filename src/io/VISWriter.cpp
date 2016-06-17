#include "io/VISWriter.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "Common.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

using Log::global_log;
using namespace std;

VISWriter::VISWriter(unsigned long writeFrequency, unsigned long writeFrequencyFile, string outputPrefix) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_writeFrequencyFile = writeFrequencyFile;
	_wroteVIS = false;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

VISWriter::~VISWriter(){}

void VISWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void VISWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain) {
	/*
    string filename = _outputPrefix + ".vis";
	ofstream fileout(filename.c_str(), ios::out);
	fileout.close();
	*/
}

void VISWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu) {
	if (simstep % _writeFrequency == 0) {
		stringstream filenamestream, outputstream;
		filenamestream << _outputPrefix;

		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}
		filenamestream << "_" << std::setw(10) << std::setfill('0') << (((simstep-_writeFrequency)/_writeFrequencyFile)+1)*_writeFrequencyFile << ".vis";
		
		char filename[filenamestream.str().size()+1];
		strcpy(filename,filenamestream.str().c_str());

		// cout << "filename = " << filename << endl;

#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		int numprocs = domainDecomp->getNumProcs();
		if (rank== 0){
#endif
			//if (!_wroteVIS){
			if ((simstep-_writeFrequency)%_writeFrequencyFile == 0){
				//outputstream << "      id t          x          y          z     q0     q1     q2     q3        c\n";
				outputstream << "     id t      x      y      z     q0     q1     q2     q3\n";
				_wroteVIS = true;
			}
			else
				outputstream << "#" << endl;
#ifdef ENABLE_MPI
		}
#endif

		// originally VIS files had a fixed width of 8 (and no t), here I use 12 (with 2 for t)
		//ostrm << "t           x           y           z          q0          q1          q2          q3" << endl;
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
//			if (!halo) {
//				outputstream << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(2)
//				            << pos->componentid() << setprecision(3);
//				for (unsigned short d = 0; d < 3; d++) outputstream << setw(11) << pos->r(d);
//				outputstream << setprecision(3) << setw(7) << pos->q().qw() << setw(7) << pos->q().qx()
//				            << setw(7) << pos->q().qy()<< setw(7) << pos->q().qz()
//				            << setw(9) << right << 0 << "\n";
//			}
			if (!halo) {
				outputstream << setiosflags(ios::fixed) << setw(7) << pos->id() << setw(2)
				            << pos->componentid() << setprecision(1);
				for (unsigned short d = 0; d < 3; d++) outputstream << setw(7) << pos->r(d);
				outputstream << setprecision(3) << setw(7) << pos->q().qw() << setw(7) << pos->q().qx()
				            << setw(7) << pos->q().qy()<< setw(7) << pos->q().qz() << endl;
			}
		}
		long outputsize = outputstream.str().size();

		char output[outputsize+1];
		strcpy(output,outputstream.str().c_str());
#ifdef ENABLE_MPI
		MPI_File fh;

		if ((simstep-_writeFrequency)%_writeFrequencyFile == 0)
		{
			MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
		}
		else
		{
			MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh);
		}

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

		MPI_File_seek(fh, offset, MPI_SEEK_END);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_File_write(fh, output, outputsize, MPI_CHAR, &status);
		MPI_File_close(&fh);
#else
		if ((simstep-_writeFrequency)%1000 == 0)
		{
			ofstream fileout(filename, ios::out);
			fileout << output;
			fileout.close();
		}
		else
		{
			ofstream fileout(filename, ios::out|ios::app);
			fileout << output;
			fileout.close();
		}
#endif
	}
}

void VISWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}
