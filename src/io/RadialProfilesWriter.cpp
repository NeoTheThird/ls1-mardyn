// RadialProfilesWriter.cpp

#include "io/RadialProfilesWriter.h"
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

RadialProfilesWriter::RadialProfilesWriter( unsigned long writeFrequency, unsigned long updateFrequency, unsigned long writeFreqAverage,
									        unsigned int nNumShells, unsigned int nNumDiscretSteps, unsigned long nNumDetectionTimesteps,
										    bool incremental )
{
	_writeFrequency   = writeFrequency;
	_updateFrequency  = updateFrequency;
	_writeFreqAverage = writeFreqAverage;

	_nNumShells = nNumShells;
	_nNumDiscretSteps = nNumDiscretSteps;
	_nNumDetectionTimesteps = nNumDetectionTimesteps;
	_nNumTimestepsDetected  = 0;
	_dVeloMax = 0;  // max value has to be determined first

	_incremental = incremental;
	_bMaxVeloDetermined = false;
	_bDiscretisationDone = false;
	_bWroteHeaderVmax = false;

	_veloDistrMatrix = NULL;
	_veloDistrMatrixAve = NULL;

	_appendTimestamp = false;
}

RadialProfilesWriter::~RadialProfilesWriter()
{
	// TODO: Speicher von dynamischen arrays freigeben
}

void RadialProfilesWriter::initOutput( ParticleContainer* particleContainer,
                                    DomainDecompBase* domainDecomp, Domain* domain)
{
	/*
	string filename = _filename + ".tra";
	ofstream fileout(filename.c_str(), ios::out);
	fileout.close();
	_wroteTra = false;
	*/
}

void RadialProfilesWriter::readXML(XMLfileUnits& xmlconfig)
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

void RadialProfilesWriter::doOutput( ParticleContainer* particleContainer,
									 DomainDecompBase* domainDecomp, Domain* domain,
									 unsigned long simstep, list<ChemicalPotential>* lmu)
{
	double dVeloMax;

	// max value
	if(_bMaxVeloDetermined == false)
		dVeloMax = CalcMaxVelocity( particleContainer, domainDecomp, domain );

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

	// writing .max-files
	std::stringstream outputstream;
	std::stringstream filenamestream;
	filenamestream << "velocity_quare.max";
	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

	#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank== 0)
		{
	#endif
			if( _bWroteHeaderVmax == false)
			{
				outputstream << "TS           v2_max           \n";
				_bWroteHeaderVmax = true;
			}

			outputstream << simstep << std::setw(11) << std::setprecision(6) << dVeloMax << std::endl;

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

	// update profiles with respect to update frequency
	if ( simstep % _updateFrequency == 0 )
	{
		CalcVelocityDistributaion( particleContainer, domainDecomp, domain);
	}

	// write out profiles with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;


}

void RadialProfilesWriter::finishOutput( ParticleContainer* particleContainer,
                                         DomainDecompBase* domainDecomp, Domain* domain)
{
}

void RadialProfilesWriter::CalcVelocityDistributaion( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain)
{
	// implementation
}

double RadialProfilesWriter::CalcMaxVelocity( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain )
{
	double dVelo, dVeloMaxLocal, dVeloMaxGlobal;
	dVeloMaxLocal = 0;
	dVeloMaxGlobal = 0;

	for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next() )
	{
		dVelo = pos->v2();
		dVeloMaxLocal = dVeloMaxLocal + (dVelo - dVeloMaxLocal) * (dVelo > dVeloMaxLocal);

		// global_log->info() << "dVelo: " << dVelo << endl;
		// global_log->info() << "_dVeloMax: " << _dVeloMax << endl;
		// global_log->info() << "dVeloMaxLocal: " << dVeloMaxLocal << endl;
	}

	#ifdef ENABLE_MPI

		MPI_Allreduce( &dVeloMaxLocal, &dVeloMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if(_bMaxVeloDetermined == false)
			_dVeloMax = _dVeloMax + (dVeloMaxGlobal - _dVeloMax) * (dVeloMaxGlobal > _dVeloMax);

	#else

		if(_bMaxVeloDetermined == false)
			_dVeloMax = _dVeloMax + (dVeloMaxLocal - _dVeloMax) * (dVeloMaxLocal > _dVeloMax);

	#endif

	_nNumTimestepsDetected++;
	_bMaxVeloDetermined = _nNumTimestepsDetected >= _nNumDetectionTimesteps;

	#ifdef ENABLE_MPI

		return dVeloMaxGlobal;

	#else

		return dVeloMaxLocal;

	#endif
}














