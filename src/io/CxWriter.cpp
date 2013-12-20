// CxWriter.cpp

#include "io/CxWriter.h"
#include "io/xyVal.h"
#include "io/Region.h"
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

CxWriter::CxWriter( unsigned long writeFrequency, unsigned long updateFrequency, double dDeltaX, unsigned int nDeltaXSteps,
		            bool incremental)
{
	_writeFrequency = writeFrequency;
	_updateFrequency = updateFrequency;

	_incremental = incremental;
	_filenameisdate = false;
	_appendTimestamp = false;

	_dDeltaX = dDeltaX;
	_nDeltaXSteps = nDeltaXSteps;

	// init number of molecules help variables;
	_nNumMolsSlice = 0;
	_nNumMolsSliceComp = 0;
}

CxWriter::~CxWriter(){}

void CxWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain)
{
	/*
	string filename = _filename + ".tra";
	ofstream fileout(filename.c_str(), ios::out);
	fileout.close();
	_wroteTra = false;
	*/
}

void CxWriter::readXML(XMLfileUnits& xmlconfig)
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

void CxWriter::doOutput( ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu)
{
	std::vector<xyVal*>::iterator it;
	double dX, dC;

	// update cx-profile with respect to update frequency
	if ( simstep % _updateFrequency == 0 )
	{
		// calculate new concentration profile
		CalcConcentrAtX( particleContainer, domainDecomp, domain, 1);
	}

	// write out cx-profile with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;

#ifdef ENABLE_MPI
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
#endif

	// writing .cox-file
	std::stringstream outputstream;
	std::stringstream filenamestream;
	filenamestream << "cx-profile_" << simstep << ".cox";
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

			for(it = _Cx.begin(); it != _Cx.end(); it++)
			{
				dX = (*it)->GetXvalue();
				dC = (*it)->GetYvalue();
				outputstream << std::setw(11) << std::setprecision(6) << dX;
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

void CxWriter::finishOutput( ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp, Domain* domain)
{
}

void CxWriter::CalcConcentrAtX( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned int nCompID)
{
	unsigned long nNumSlices;
	list<Molecule*> particlePtrs;
	list<Molecule*>::iterator it;
	double dNumMolsSlice;
	double dNumMolsSliceComp;
	double dDeltaXSteps = (double) _nDeltaXSteps;
	double dX, dC;

	global_log->info() << "dDeltaXSteps" << dDeltaXSteps << endl;

	nCompID--; // cid in input file starts with 1, in ls1 with 0

	double dBoxLength[3];
	dBoxLength[0] = domain->getGlobalLength(0);
	dBoxLength[1] = domain->getGlobalLength(1);
	dBoxLength[2] = domain->getGlobalLength(2);

	// clear Cx-vector
	this->ClearCxVector();

	// Anzahl der Slices
	nNumSlices = dBoxLength[0] / _dDeltaX;

	double dLowerCorner[3];
	double dUpperCorner[3];
	dLowerCorner[0] = 0.0;
	dLowerCorner[1] = 0.0;
	dLowerCorner[2] = 0.0;
	dUpperCorner[0] = _dDeltaX;
	dUpperCorner[1] = dBoxLength[1];
	dUpperCorner[2] = dBoxLength[2];

	global_log->info() << "calculating concentration profile..." << endl;
//	global_log->info() << "_dDeltaX: " << _dDeltaX << endl;

	for(unsigned long i=0; i < nNumSlices * _nDeltaXSteps; i++)
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
		dX = ( dUpperCorner[0] + dLowerCorner[0] ) / 2.0;
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
		_Cx.push_back(new xyVal(dX, dC) );

		// Begrenzung nächster Slice
		dLowerCorner[0] += _dDeltaX / dDeltaXSteps;
		dUpperCorner[0] += _dDeltaX / dDeltaXSteps;
	}

	// Update region list with respect to new concentration profile
	this->UpdateRegionList(domain);
}

void CxWriter::ClearCxVector(void)
{
	std::vector<xyVal*>::iterator it;

	for(it = _Cx.begin(); it != _Cx.end(); it++)
	{
		delete (*it);
		(*it) = NULL;
	}

	_Cx.clear();
}

int CxWriter::UpdateRegionList(Domain* domain)
{
	int nNumRegions = 0;
	unsigned int nNumValues;
	std::vector<xyVal*>::iterator it;
	std::vector<xyVal*>::iterator itLast;
	double dX, dC;
	double dLowerCorner[3];
	double dUpperCorner[3];
	int nPresentPhasesState = PPS_UNKNOWN;
	int nPresentPhasesStateLast = PPS_UNKNOWN;
	int UpdateRegionListState = URLS_INITIAL_STEP;
	double dDeltaXSteps = (double) _nDeltaXSteps;
	double dBoxLengthX = domain->getGlobalLength(0);

	// last element
	itLast = _Cx.end();
	itLast--;

	// init y, z length
	dLowerCorner[1] = 0.0;
	dLowerCorner[2] = 0.0;
	dUpperCorner[1] = domain->getGlobalLength(1);
	dUpperCorner[2] = domain->getGlobalLength(2);

	global_log->info() << "CxWriter::UpdateRegionList - updating region list." << endl;

	// clear list
	this->ClearRegionList();

	unsigned int i = 0;
	// unsigned int nSize = _Cx.size();

	for(it = _Cx.begin(); it != _Cx.end(); it++)
	{
		dX = (*it)->GetXvalue();
		dC = (*it)->GetYvalue();

		// check if end of domain reached
		if(it == itLast)
			UpdateRegionListState = URLS_LAST_ELEMENT;

		// check present phases state
		if(dC == 1.0)
			nPresentPhasesState = PPS_COMPONENT_ONE_ONLY;
		else if(dC == 0.0)
			nPresentPhasesState = PPS_COMPONENT_TWO_ONLY;
		else if(dC > 0.0 && dC < 1.0)
			nPresentPhasesState = PPS_COMPONENT_ONE_AND_TWO_MIXED;
		else
			nPresentPhasesState = PPS_UNKNOWN;

		switch(UpdateRegionListState)
		{
		case URLS_INITIAL_STEP:
//			global_log->info() << "state: URLS_INITIAL_STEP" << endl;
			dLowerCorner[0] = 0.0;
			nNumValues = 1;
			nPresentPhasesStateLast = nPresentPhasesState;
			UpdateRegionListState = URLS_LEFT_BOUNDARY_SET;
			break;
//		case URLS_LEFT_BOUNDARY_SET:
//			global_log->info() << "state: URLS_LEFT_BOUNDARY_SET" << endl;
//			global_log->info() << "nPresentPhasesStateLast: " << nPresentPhasesStateLast << endl;
//			global_log->info() << "nPresentPhasesState: " << nPresentPhasesState << endl;
//			if(nPresentPhasesStateLast == nPresentPhasesState && nPresentPhasesState != PPS_UNKNOWN)
//				nNumValues++;
//			else
//				UpdateRegionListState = URLS_INITIAL_STEP;
//
//			if(nNumValues == nNumValuesMin)
//				UpdateRegionListState = URLS_NUM_VALUES_MIN_REACHED;
//			break;
		case URLS_LEFT_BOUNDARY_SET:
//			global_log->info() << "state: URLS_NUM_VALUES_MIN_REACHED" << endl;
//			global_log->info() << "nPresentPhasesStateLast: " << nPresentPhasesStateLast << endl;
//			global_log->info() << "nPresentPhasesState: " << nPresentPhasesState << endl;
			if( nPresentPhasesStateLast == nPresentPhasesState && nPresentPhasesState != PPS_UNKNOWN)
				nNumValues++;
			else
			{
				dUpperCorner[0] = dX - _dDeltaX / dDeltaXSteps * 0.5;
				_regionList.push_back( new Region(dLowerCorner, dUpperCorner, nPresentPhasesStateLast) );
				nNumRegions++;
				global_log->info() << "added region nr: " << nNumRegions << endl;

				global_log->info() << "LowerCorner x: " << dLowerCorner[0] <<  ", y: " << dLowerCorner[1] <<  ", z: " << dLowerCorner[2] << endl;
				global_log->info() << "UpperCorner x: " << dUpperCorner[0] <<  ", y: " << dUpperCorner[1] <<  ", z: " << dUpperCorner[2] << endl;

				UpdateRegionListState = URLS_RIGHT_BOUNDARY_SET;
			}
			break;
		case URLS_RIGHT_BOUNDARY_SET:
//			global_log->info() << "state: URLS_INITIAL_STEP" << endl;
			dLowerCorner[0] = dX - _dDeltaX / dDeltaXSteps * 1.5;
			nNumValues = 1;
			nPresentPhasesStateLast = nPresentPhasesState;
			UpdateRegionListState = URLS_LEFT_BOUNDARY_SET;
			break;

		case URLS_LAST_ELEMENT:
			dUpperCorner[0] = dBoxLengthX;
			_regionList.push_back( new Region(dLowerCorner, dUpperCorner, nPresentPhasesStateLast) );
			nNumRegions++;
			global_log->info() << "added region nr: " << nNumRegions << endl;

			global_log->info() << "LowerCorner x: " << dLowerCorner[0] <<  ", y: " << dLowerCorner[1] <<  ", z: " << dLowerCorner[2] << endl;
			global_log->info() << "UpperCorner x: " << dUpperCorner[0] <<  ", y: " << dUpperCorner[1] <<  ", z: " << dUpperCorner[2] << endl;

			UpdateRegionListState = URLS_END_OF_DOMAIN_REACHED;
			break;
		case URLS_END_OF_DOMAIN_REACHED:
			break;
		}

		i++;
	}

	// update domain region list
	domain->UpdateRegionList(_regionList);

	return nNumRegions;
}

void CxWriter::ClearRegionList(void)
{
	std::vector<Region*>::iterator it;

	for(it = _regionList.begin(); it != _regionList.end(); it++)
	{
		delete (*it);
		(*it) = NULL;
	}

	_regionList.clear();
}




















