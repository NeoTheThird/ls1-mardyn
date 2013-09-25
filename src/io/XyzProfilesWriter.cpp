// XyzProfilesWriter.cpp

#include "io/XyzProfilesWriter.h"
#include "io/xyVal.h"
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

XyzProfilesWriter::XyzProfilesWriter( unsigned long writeFrequency, unsigned long updateFrequency, double dSliceWidth, unsigned int nMoveSteps,
                                      unsigned int* nRecSwitches, bool incremental )
{
	_writeFrequency = writeFrequency;
	_updateFrequency = updateFrequency;

	_incremental = incremental;
	_filenameisdate = false;
	_appendTimestamp = false;

	_dSliceWidth = dSliceWidth;
	_nMoveSteps = nMoveSteps;

	// init number of molecules help variables;
	_nNumMolsSlice = 0;
	_nNumMolsSliceComp = 0;

	// init profile record switches
	for(int i = 0; i < 6; i++)
	{
		_nRecSwitches[i] = nRecSwitches[i];
	}
}

XyzProfilesWriter::~XyzProfilesWriter(){}

void XyzProfilesWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain)
{
	/*
	string filename = _filename + ".tra";
	ofstream fileout(filename.c_str(), ios::out);
	fileout.close();
	_wroteTra = false;
	*/
}

void XyzProfilesWriter::readXML(XMLfileUnits& xmlconfig)
{
	/*
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
*/
}

void XyzProfilesWriter::doOutput( ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu)
{
	std::vector<xyVal*>::iterator it;
	double dPos, dC;

	// update profiles with respect to update frequency
	if ( simstep % _updateFrequency == 0 )
	{
		if(_nRecSwitches[0] == 1)
		{
			// calculate new concentration profile, dimension: X
			CalcConcentrProfile( particleContainer, domainDecomp, domain, 1, XYZPFD_DIMENSION_X);
		}

		if(_nRecSwitches[1] == 1)
		{
			// calculate new concentration profile, dimension: Y
			CalcConcentrProfile( particleContainer, domainDecomp, domain, 1, XYZPFD_DIMENSION_Y);
		}

		if(_nRecSwitches[2] == 1)
		{
			// calculate new concentration profile, dimension: Z
			CalcConcentrProfile( particleContainer, domainDecomp, domain, 1, XYZPFD_DIMENSION_Z);
		}
	}

	// write out profiles with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;

	struct rlimit limit;
	/* Set the stack limit in bytes. */
	//
	limit.rlim_cur = 1099511627776;
	limit.rlim_max = 1099511627776;
	if (setrlimit(RLIMIT_STACK, &limit) != 0)
	{
		cout << "setrlimit() failed.";
		exit(1);
	}

	// writing .prf-files
	std::stringstream outputstream;
	std::stringstream filenamestream;
	filenamestream << "profile_" << simstep << ".prf";
	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

	#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank== 0)
		{
	#endif
			outputstream << "x           c           \n";

			// global_log->info() << "concentrAtX size: " << _Cx.size() << endl;

			for(it = _Cxyz[0].begin(); it != _Cxyz[0].end(); it++)
			{
				dPos = (*it)->GetXvalue();
				dC = (*it)->GetYvalue();
				outputstream << std::setw(11) << std::setprecision(6) << dPos;
				outputstream << std::setw(11) << std::setprecision(6) << dC << endl;
			}

			// typumwandlung
			long outputsize = outputstream.str().size();
			//cout << "rank: " << rank << "; step: " << simstep << "; outputsize: " << outputsize << endl;
			char output[outputsize+1];
			strcpy(output,outputstream.str().c_str());

			// Datei zum schreiben öffnen, daten schreiben
			ofstream fileout(filename, ios::out|ios::app);
			fileout << output;
			fileout.close();
	#ifdef ENABLE_MPI
		}
	#endif
}

void XyzProfilesWriter::finishOutput( ParticleContainer* particleContainer,
                                      DomainDecompBase* domainDecomp, Domain* domain)
{
}

void XyzProfilesWriter::CalcConcentrProfile( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
		                                     unsigned int nCompID, unsigned int nDimension)
{
	unsigned long nNumSlices;
	list<Molecule*> particlePtrs;
	list<Molecule*>::iterator it;
	double dNumMolsSlice;
	double dNumMolsSliceComp;
	double dMoveSteps = (double) _nMoveSteps;
	double dPos, dC;

	global_log->info() << "dMoveSteps" << dMoveSteps << endl;

	nCompID--; // cid in input file starts with 1, in ls1 with 0

	double dBoxLength[3];
	dBoxLength[0] = domain->getGlobalLength(0);
	dBoxLength[1] = domain->getGlobalLength(1);
	dBoxLength[2] = domain->getGlobalLength(2);

	// clear Cx-vector
	this->ClearCxyzVector(nDimension);

	// Anzahl der Slices
	nNumSlices = dBoxLength[nDimension] / _dSliceWidth;

	double dLowerCorner[3];
	double dUpperCorner[3];
	dLowerCorner[0] = 0.0;
	dLowerCorner[1] = 0.0;
	dLowerCorner[2] = 0.0;
	dUpperCorner[0] = dBoxLength[0];
	dUpperCorner[1] = dBoxLength[1];
	dUpperCorner[2] = dBoxLength[2];

	dUpperCorner[nDimension] = _dSliceWidth;  // only one dimension divided into slices

	global_log->info() << "calculating concentration profile..." << endl;
//	global_log->info() << "_dDeltaX: " << _dDeltaX << endl;

	for(unsigned long i=0; i < nNumSlices * _nMoveSteps; i++)
	{
		particlePtrs.clear();
		particleContainer->getRegion(dLowerCorner, dUpperCorner, particlePtrs);

		// global_log->info() << "particlePtrs - size: " << particlePtrs.size() << endl;

		_nNumMolsSlice = 0;
		_nNumMolsSliceComp = 0;

		for(it = particlePtrs.begin(); it != particlePtrs.end(); ++it)
		{
			_nNumMolsSlice++;
			unsigned int cid = (*it)->componentid();

			// global_log->info() << "cid: " << cid << endl;

			if(cid == nCompID)
				_nNumMolsSliceComp++;
		}

		// Wert akkumulieren mit anderen Prozessen
		domainDecomp->collCommInit(2);
		domainDecomp->collCommAppendInt(_nNumMolsSlice);
		domainDecomp->collCommAppendInt(_nNumMolsSliceComp);
		domainDecomp->collCommAllreduceSum();
		_nNumMolsSlice = domainDecomp->collCommGetInt();
		_nNumMolsSliceComp = domainDecomp->collCommGetInt();
		domainDecomp->collCommFinalize();

//		global_log->info() << "calculating x ..." << endl;
		// dX = _dDeltaX / 2.0 + _dDeltaX * (double)(i);
		dPos = ( dUpperCorner[nDimension] + dLowerCorner[nDimension] ) / 2.0;
//		global_log->info() << "x: " << dX << endl;

//		global_log->info() << "calculating c ..." << endl;
//		global_log->info() << "_nNumMolsSlice: " << _nNumMolsSlice << endl;
//		global_log->info() << "_nNumMolsSliceComp: " << _nNumMolsSliceComp << endl;
		dNumMolsSlice = (double) _nNumMolsSlice;
		dNumMolsSliceComp = (double) _nNumMolsSliceComp;
		dC = dNumMolsSliceComp / dNumMolsSlice;
//		global_log->info() << "c: " << dC << endl;
//		global_log->info() << "dNumMolsSlice: " << dNumMolsSlice << endl;
//		global_log->info() << "dNumMolsSliceComp: " << dNumMolsSliceComp << endl;

//		global_log->info() << "adding value ..." << endl;
		_Cxyz[nDimension].push_back(new xyVal(dPos, dC) );

		// Begrenzung nächster Slice
		dLowerCorner[nDimension] += _dSliceWidth / dMoveSteps;
		dUpperCorner[nDimension] += _dSliceWidth / dMoveSteps;
	}
}

void XyzProfilesWriter::ClearCxyzVector(unsigned int nDimension)
{
	std::vector<xyVal*>::iterator it;

	for(it = _Cxyz[nDimension].begin(); it != _Cxyz[nDimension].end(); it++)
	{
		delete (*it);
		(*it) = NULL;
	}

	_Cxyz[nDimension].clear();
}

















