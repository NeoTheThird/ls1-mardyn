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

//	_veloDistrMatrix = NULL;
//	_veloDistrMatrixAve = NULL;

	_appendTimestamp = false;

	_dAverageFactor = 1.0;
}

RadialProfilesWriter::~RadialProfilesWriter()
{
	// delete all dynamic allocated arrays
	delete[] _dDiscreteRadiusValues;
	delete[] _dDiscreteVelocityValues;
	delete[] _dDensityProfile;
	delete[] _dShellVolumesReduced;

	delete[]  _nNumMoleculesInsideShellLocal;
#ifdef ENABLE_MPI
	delete[] _nNumMoleculesInsideShellGlobal;
#endif

	// delete all dynamic allocated 2d-arrays
	for(int i=0; i < _nNumShells; i++)
	{
		delete[] _veloDistrMatrix[i];
		delete[] _veloDistrMatrixAve[i];
	}
	delete[] _veloDistrMatrix;
	delete[] _veloDistrMatrixAve;

	// set all pointers to NULL
	_dDiscreteRadiusValues = NULL;
	_dDiscreteVelocityValues = NULL;
	_dDensityProfile = NULL;
	_dShellVolumesReduced = NULL;

	_nNumMoleculesInsideShellLocal = NULL;
#ifdef ENABLE_MPI
	_nNumMoleculesInsideShellGlobal = NULL;
#endif

	_veloDistrMatrix = NULL;
	_veloDistrMatrixAve = NULL;
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


	// max value
	if(_bMaxVeloDetermined == false)
	{
		dVeloMax = CalcMaxVelocity( particleContainer, domainDecomp, domain );

	// writing .max-files
	std::stringstream outputstream;
	std::stringstream filenamestream;
	filenamestream << "velocity_square.max";
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

	}  // if(_bMaxVeloDetermined == false)


	if(_bDiscretisationDone == false && _bMaxVeloDetermined == true)
	{
		DoDiscretisation(domain);  // initialise discretisation data structure
		global_log->info() << "Discretisationn done!" << endl;
	}


	if(_bDiscretisationDone == false)
		return;

	// update profiles with respect to update frequency
	if ( !(simstep % _updateFrequency == 0) )
		return;

	// start averaging after first simstep
	if(simstep > 0)
		_dAverageFactor = 0.5;

	// calculate radial profiles
	global_log->info() << "calculating radial profiles..." << endl;

	CalcRadialProfiles( particleContainer, domainDecomp, domain );

	global_log->info() << "... calculation done!" << endl;



	// write out profiles with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;

	global_log->info() << "writing radial density profile..." << endl;

	// writing .rpf-files
	std::stringstream outputstream, outputstreamVelo;
	std::stringstream filenamestream, filenamestreamVelo;
	filenamestream << "radial_density_profile_averaged_TS" << simstep << ".rpf";
	filenamestreamVelo << "velocity_distribution_TS" << simstep << ".rpf";
	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());
	char filenameVelo[filenamestreamVelo.str().size()+1];
	strcpy(filenameVelo,filenamestreamVelo.str().c_str());

	#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank== 0)
		{
	#endif

			// radial density profile
			outputstream << "r           rho           \n";

			for(int i = 0; i < _nNumShells; i++)
			{
				outputstream << _dDiscreteRadiusValues[i] << std::setw(11) << std::setprecision(6) << _dDensityProfile[i] << std::endl;
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
			// radial density profile


			global_log->info() << "writing velocity distribution..." << endl;


			// radial velocity distribution
			outputstreamVelo << "v/R       ";

			// first line - discrete radius values
			for(int i = 0; i < _nNumShells; i++)
			{
				outputstreamVelo << _dDiscreteRadiusValues[i] << std::setw(11) << std::setprecision(6);
			}
			outputstreamVelo << endl;

			// velocity distribution matrix
			for(int i = 0; i < _nNumDiscretSteps; i++)
			{
				outputstreamVelo << _dDiscreteVelocityValues[i] << std::setw(11) << std::setprecision(6);

				for(int j = 0; j < _nNumShells; j++)
				{
					outputstreamVelo << _veloDistrMatrix[j][i] << std::setw(11) << std::setprecision(6);
				}
				outputstreamVelo << endl;
			}

			// typumwandlung
			long outputsize2 = outputstreamVelo.str().size();
			//cout << "rank: " << rank << "; step: " << simstep << "; outputsize: " << outputsize << endl;
			char output2[outputsize2+1];
			strcpy(output2,outputstreamVelo.str().c_str());

			// Datei zum schreiben öffnen, daten schreiben
			ofstream fileout2(filenameVelo, ios::out|ios::app);
			fileout2 << output2;
			fileout2.close();
			// radial velocity distribution

	#ifdef ENABLE_MPI
		}
	#endif

}

void RadialProfilesWriter::finishOutput( ParticleContainer* particleContainer,
                                         DomainDecompBase* domainDecomp, Domain* domain)
{
}

void RadialProfilesWriter::CalcRadialProfiles( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain )
{
	// calculate radius vector
	double dRadiusLength;
	double dShellWidth;
	unsigned int nRadiusIndex;
	unsigned int nVelocityIndex;
	unsigned int nIndexMax = _nNumShells-1;
	unsigned int nIndexMaxVelo = _nNumDiscretSteps-1;
	double dSystemMidpoint[3];
	double dRadiusVector[3];
	double dMaxVelo = _dVeloMax * 1.2;  // velocity discretisation, highest value with safety factor


	global_log->info() << "dMaxVelo: " << dMaxVelo << endl;


	// init local density profile array
	for(int i = 0; i < _nNumShells; i++)
	{
		_nNumMoleculesInsideShellLocal[i] = 0;
	#ifdef ENABLE_MPI
		_nNumMoleculesInsideShellGlobal[i] = 0;
	#endif
	}
	// radial density profile


	// radial velocity distribution
	unsigned long** veloDistrMatrixLocal = new unsigned long*[_nNumShells];

	for(int i = 0; i < _nNumShells; i++)
	{
		veloDistrMatrixLocal[i] = new unsigned long[_nNumDiscretSteps];
	}

	// init array
	for(int i = 0; i < _nNumShells; i++)
	{
		for(int j = 0; j < _nNumDiscretSteps; j++)
		{
			veloDistrMatrixLocal[i][j] = 0;
		}
	}
	// radial velocity distribution


	// calc shell width
	dShellWidth = domain->getGlobalLength(0) * 0.5 / _nNumShells;

	dSystemMidpoint[0] = domain->getGlobalLength(0) * 0.5;
	dSystemMidpoint[1] = domain->getGlobalLength(1) * 0.5;
	dSystemMidpoint[2] = domain->getGlobalLength(2) * 0.5;

//	global_log->info() << "dSystemMidpoint[0] " << dSystemMidpoint[0] << endl;
//	global_log->info() << "dSystemMidpoint[1] " << dSystemMidpoint[1] << endl;
//	global_log->info() << "dSystemMidpoint[2] " << dSystemMidpoint[2] << endl;
//
//	global_log->info() << "preparing density profile calculation..." << endl;

	/*
	// prepare to calc average with new time step
	for(int i = 0; i < _nNumShells; i++)
	{
		_dDensityProfile[i] = _dDensityProfile[i] * 0.5;
	}
	*/
//	global_log->info() << "...preparation done!" << endl;

	for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next() )
	{
		// calc radius vector
		dRadiusVector[0] = pos->r(0) - dSystemMidpoint[0];
		dRadiusVector[1] = pos->r(1) - dSystemMidpoint[1];
		dRadiusVector[2] = pos->r(2) - dSystemMidpoint[2];

//		global_log->info() << "pos->r(0) :" << pos->r(0) << endl;
//		global_log->info() << "pos->r(1) :" << pos->r(1) << endl;
//		global_log->info() << "pos->r(2) :" << pos->r(2) << endl;
//
//		global_log->info() << "dRadiusVector[0] " << dRadiusVector[0] << endl;
//		global_log->info() << "dRadiusVector[1] " << dRadiusVector[1] << endl;
//		global_log->info() << "dRadiusVector[2] " << dRadiusVector[2] << endl;

		// calc length of radius ^2
		dRadiusLength = sqrt(dRadiusVector[0] * dRadiusVector[0] + dRadiusVector[1] * dRadiusVector[1] + dRadiusVector[2] * dRadiusVector[2]);

//		global_log->info() << "dRadiusLength " << dRadiusLength << endl;
//		global_log->info() << "dShellWidth " << dShellWidth << endl;

		// radius^2 --> radiusIndex
		nRadiusIndex = (unsigned int) (dRadiusLength / dShellWidth);

//		global_log->info() << "nRadiusIndex " << pos->id() << ": " << nRadiusIndex << endl;

		/*
		// calc new average density profile
		if(nRadiusIndex < nIndexMax)
		{
			_dDensityProfile[nRadiusIndex] += 0.5 / _dShellVolumesReduced[nRadiusIndex];
		}
		*/

		// radial density profile
		if(nRadiusIndex < nIndexMax)  // respect finite resolution of radius
		{
			_nNumMoleculesInsideShellLocal[nRadiusIndex]++;
		}
		// radial density profile


		// radial velocity distribution
		nVelocityIndex = (unsigned int) (pos->v2() / dMaxVelo * _nNumDiscretSteps);

		// global_log->info() << "nVelocityIndex: " << nVelocityIndex << endl;

		if(nVelocityIndex < nIndexMaxVelo && nRadiusIndex < nIndexMax)  // respect finite resolution of radius and velocity
		{
			veloDistrMatrixLocal[nRadiusIndex][nVelocityIndex]++;
		}
		// radial velocity distribution

	}  // foreach molecule

	#ifdef ENABLE_MPI

		// radial density profile

		//MPI_Reduce( dDensityProfileLocal, _dDensityProfile, _nNumShells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _nNumMoleculesInsideShellLocal, _nNumMoleculesInsideShellGlobal, _nNumShells, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		for(int i = 0; i < _nNumShells; i++)
		{
			_dDensityProfile[i] = (_dDensityProfile[i] + _nNumMoleculesInsideShellGlobal[i] / _dShellVolumesReduced[i] ) * _dAverageFactor;  // _dAverageFactor == 1.0 or 0.5
		}
		// radial density profile


		// radial velocity distribution

		for(int i = 0; i < _nNumShells; i++)
		{
			MPI_Reduce( veloDistrMatrixLocal[i], _veloDistrMatrix[i], _nNumShells, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		// radial velocity distribution

	#else

		// radial density profile
		for(int i = 0; i < _nNumShells; i++)
		{
			_dDensityProfile[i] = (_dDensityProfile[i] + _nNumMoleculesInsideShellLocal[i] / _dShellVolumesReduced[i] ) * _dAverageFactor;  // _dAverageFactor == 1.0 or 0.5
		}
		// radial density profile


		// radial velocity distribution
		for(int i = 0; i < _nNumShells; i++)
		{
			for(int j = 0; j < _nNumDiscretSteps; j++)
			{
				_veloDistrMatrix[i][j] = veloDistrMatrixLocal[i][j];
			}
		}
		// radial velocity distribution

	#endif

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


void RadialProfilesWriter::DoDiscretisation(Domain* domain)
{
	if(_bDiscretisationDone == true)  // if allready done -> return
		return;

	double pi = 3.141592654;

	double dVeloMax = _dVeloMax * 1.2;  // velocity discretisation, highest value with safety factor
	double dRadiusMax = domain->getGlobalLength(0) * 0.5;
	double dShellWidth = dRadiusMax / _nNumShells;

	double dRi, dRa;

	// velocity distribution
	_veloDistrMatrix = new unsigned long* [_nNumShells];

	for(int i = 0; i < _nNumShells; i++)
	{
		_veloDistrMatrix[i] = new unsigned long[_nNumDiscretSteps];
	}

	// init matrix
	for(int i = 0; i < _nNumShells; i++)
	{
		for(int j = 0; j < _nNumDiscretSteps; j++)
		{
			_veloDistrMatrix[i][j] = 0;
		}
	}

	// radial density profile
	_nNumMoleculesInsideShellLocal = new unsigned long[_nNumShells];
#ifdef ENABLE_MPI
	_nNumMoleculesInsideShellGlobal = new unsigned long[_nNumShells];
#endif

	// densityProfile
	_dDensityProfile = new double[_nNumShells];

	for(int i = 0; i < _nNumShells; i++)
	{
		_dDensityProfile[i] = 0.0;
		_nNumMoleculesInsideShellLocal[i] = 0;

#ifdef ENABLE_MPI
		_nNumMoleculesInsideShellGlobal[i] = 0;
#endif
	}


	// calc discrete radius values
	_dDiscreteRadiusValues = new double[_nNumShells];

	for(int i = 0; i < _nNumShells; i++)
	{
		_dDiscreteRadiusValues[i] = (i + 0.5) * dShellWidth;
	}


	// calc reduced shell volumes
	_dShellVolumesReduced = new double[_nNumShells];

	for(int i = 0; i < _nNumShells; i++)
	{
		dRi = dShellWidth * i;
		dRa = dShellWidth * (i+1);
		_dShellVolumesReduced[i] = pi * 4. / 3. * (dRa*dRa*dRa - dRi*dRi*dRi);  // V_shell = 4/3 pi (Ra^3 - Ri^3)
	}


	// calc discrete velocity values
	_dDiscreteVelocityValues = new double[_nNumDiscretSteps];
	double dDeltaVelo = dVeloMax / _nNumDiscretSteps;

	for(int i = 0; i < _nNumDiscretSteps; i++)
	{
		_dDiscreteVelocityValues[i] = dDeltaVelo * (i + 0.5);
	}

	_bDiscretisationDone = true;
}











