// RadialProfilesWriter.cpp

#include "io/RadialProfilesWriter.h"
#include "io/xyVal.h"
#include "utils/Vector3d.h"
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

	_appendTimestamp = false;

	_dAverageFactor = 1.0;


	// radial velocity distribution
	_veloDistrMatrixLocal    = new unsigned long*[_nNumShells];
	_veloDistrMatrixLocal_r  = new unsigned long*[_nNumShells];
	_veloDistrMatrixLocal_t1 = new unsigned long*[_nNumShells];
	_veloDistrMatrixLocal_t2 = new unsigned long*[_nNumShells];

	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		_veloDistrMatrixLocal[i]    = new unsigned long[_nNumDiscretSteps];
		_veloDistrMatrixLocal_r[i]  = new unsigned long[_nNumDiscretSteps];
		_veloDistrMatrixLocal_t1[i] = new unsigned long[_nNumDiscretSteps];
		_veloDistrMatrixLocal_t2[i] = new unsigned long[_nNumDiscretSteps];
	}

#ifdef ENABLE_MPI
	_veloDistrMatrixGlobal    = new unsigned long*[_nNumShells];
	_veloDistrMatrixGlobal_r  = new unsigned long*[_nNumShells];
	_veloDistrMatrixGlobal_t1 = new unsigned long*[_nNumShells];
	_veloDistrMatrixGlobal_t2 = new unsigned long*[_nNumShells];

	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		_veloDistrMatrixGlobal[i]    = new unsigned long[_nNumDiscretSteps];
		_veloDistrMatrixGlobal_r[i]  = new unsigned long[_nNumDiscretSteps];
		_veloDistrMatrixGlobal_t1[i] = new unsigned long[_nNumDiscretSteps];
		_veloDistrMatrixGlobal_t2[i] = new unsigned long[_nNumDiscretSteps];
	}
#endif
	// radial velocity distribution


	// radial density profile
	_nNumMoleculesInsideShellLocal = new unsigned long[_nNumShells];
#ifdef ENABLE_MPI
	_nNumMoleculesInsideShellGlobal = new unsigned long[_nNumShells];
#endif

	_dDensityProfile = new double[_nNumShells];
	// radial density profile
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
	for(unsigned int i=0; i < _nNumShells; i++)
	{
		delete[] _veloDistrMatrixLocal[i];
		delete[] _veloDistrMatrixLocal_r[i];
		delete[] _veloDistrMatrixLocal_t1[i];
		delete[] _veloDistrMatrixLocal_t2[i];

#ifdef ENABLE_MPI
		delete[] _veloDistrMatrixGlobal[i];
		delete[] _veloDistrMatrixGlobal_r[i];
		delete[] _veloDistrMatrixGlobal_t1[i];
		delete[] _veloDistrMatrixGlobal_t2[i];
#endif
	}

	delete[] _veloDistrMatrixLocal;
	delete[] _veloDistrMatrixLocal_r;
	delete[] _veloDistrMatrixLocal_t1;
	delete[] _veloDistrMatrixLocal_t2;

#ifdef ENABLE_MPI
	delete[] _veloDistrMatrixGlobal;
	delete[] _veloDistrMatrixGlobal_r;
	delete[] _veloDistrMatrixGlobal_t1;
	delete[] _veloDistrMatrixGlobal_t2;
#endif

	// set all pointers to NULL
	_dDiscreteRadiusValues = NULL;
	_dDiscreteVelocityValues = NULL;
	_dDensityProfile = NULL;
	_dShellVolumesReduced = NULL;

	_nNumMoleculesInsideShellLocal = NULL;
#ifdef ENABLE_MPI
	_nNumMoleculesInsideShellGlobal = NULL;
#endif

	_veloDistrMatrixLocal = NULL;
	_veloDistrMatrixLocal_r = NULL;
	_veloDistrMatrixLocal_t1 = NULL;
	_veloDistrMatrixLocal_t2 = NULL;

#ifdef ENABLE_MPI
	_veloDistrMatrixGlobal = NULL;
	_veloDistrMatrixGlobal_r = NULL;
	_veloDistrMatrixGlobal_t1 = NULL;
	_veloDistrMatrixGlobal_t2 = NULL;
#endif

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

	/*
	// update profiles with respect to update frequency
	if ( !(simstep % _updateFrequency == 0) )
		return;
	*/

	// start averaging after first simstep
	if(simstep > 0)
		_dAverageFactor = 0.5;

	// calculate radial profiles
	global_log->info() << "calculating radial profiles..." << endl;

	CalcRadialProfiles( particleContainer, domainDecomp, domain, simstep );

	global_log->info() << "... calculation done!" << endl;



	// write out profiles with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;

	global_log->info() << "writing radial properties files..." << endl;

	// writing .rpf-files
	std::stringstream outputstream;
	std::stringstream outputstreamVelo, outputstreamVelo_r, outputstreamVelo_t1, outputstreamVelo_t2;
	std::stringstream filenamestream;
	std::stringstream filenamestreamVelo, filenamestreamVelo_r, filenamestreamVelo_t1, filenamestreamVelo_t2;
	filenamestream << "radial_density_profile_averaged_TS" << simstep << ".rpf";
	filenamestreamVelo << "velocity_distribution_TS" << simstep << ".rpf";
	filenamestreamVelo_r << "velocity_distribution_TS" << simstep << "_r.rpf";
	filenamestreamVelo_t1 << "velocity_distribution_TS" << simstep << "_t1.rpf";
	filenamestreamVelo_t2 << "velocity_distribution_TS" << simstep << "_t2.rpf";
	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());
	char filenameVelo[filenamestreamVelo.str().size()+1];
	strcpy(filenameVelo,filenamestreamVelo.str().c_str());

	char filenameVelo_r[filenamestreamVelo_r.str().size()+1];
	strcpy(filenameVelo_r,filenamestreamVelo_r.str().c_str());
	char filenameVelo_t1[filenamestreamVelo_t1.str().size()+1];
	strcpy(filenameVelo_t1,filenamestreamVelo_t1.str().c_str());
	char filenameVelo_t2[filenamestreamVelo_t2.str().size()+1];
	strcpy(filenameVelo_t2,filenamestreamVelo_t2.str().c_str());

	#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank== 0)
		{
	#endif

			// radial density profile
			outputstream << "r           rho           \n";

			for(unsigned int i = 0; i < _nNumShells; i++)
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
			outputstreamVelo_r << "v/R       ";
			outputstreamVelo_t1 << "v/R       ";
			outputstreamVelo_t2 << "v/R       ";

			// first line - discrete radius values
			for(unsigned int i = 0; i < _nNumShells; i++)
			{
				outputstreamVelo << _dDiscreteRadiusValues[i] << std::setw(11) << std::setprecision(6);
				outputstreamVelo_r << _dDiscreteRadiusValues[i] << std::setw(11) << std::setprecision(6);
				outputstreamVelo_t1 << _dDiscreteRadiusValues[i] << std::setw(11) << std::setprecision(6);
				outputstreamVelo_t2 << _dDiscreteRadiusValues[i] << std::setw(11) << std::setprecision(6);
			}
			outputstreamVelo << endl;
			outputstreamVelo_r << endl;
			outputstreamVelo_t1 << endl;
			outputstreamVelo_t2 << endl;

			global_log->info() << "header done." << endl;

			// velocity distribution matrix
			for(unsigned int i = 0; i < _nNumDiscretSteps; i++)
			{
				outputstreamVelo << _dDiscreteVelocityValues[i] << std::setw(11) << std::setprecision(6);
				outputstreamVelo_r << _dDiscreteVelocityValues[i] << std::setw(11) << std::setprecision(6);
				outputstreamVelo_t1 << _dDiscreteVelocityValues[i] << std::setw(11) << std::setprecision(6);
				outputstreamVelo_t2 << _dDiscreteVelocityValues[i] << std::setw(11) << std::setprecision(6);

				for(unsigned int j = 0; j < _nNumShells; j++)
				{
					#ifdef ENABLE_MPI

					outputstreamVelo << _veloDistrMatrixGlobal[j][i] << std::setw(11) << std::setprecision(6);
					outputstreamVelo_r << _veloDistrMatrixGlobal_r[j][i] << std::setw(11) << std::setprecision(6);
					outputstreamVelo_t1 << _veloDistrMatrixGlobal_t1[j][i] << std::setw(11) << std::setprecision(6);
					outputstreamVelo_t2 << _veloDistrMatrixGlobal_t2[j][i] << std::setw(11) << std::setprecision(6);

					#else

					outputstreamVelo << _veloDistrMatrixLocal[j][i] << std::setw(11) << std::setprecision(6);
					outputstreamVelo_r << _veloDistrMatrixLocal_r[j][i] << std::setw(11) << std::setprecision(6);
					outputstreamVelo_t1 << _veloDistrMatrixLocal_t1[j][i] << std::setw(11) << std::setprecision(6);
					outputstreamVelo_t2 << _veloDistrMatrixLocal_t2[j][i] << std::setw(11) << std::setprecision(6);

					#endif
				}
				outputstreamVelo << endl;
				outputstreamVelo_r << endl;
				outputstreamVelo_t1 << endl;
				outputstreamVelo_t2 << endl;
			}

			global_log->info() << "matrix done." << endl;
			global_log->info() << "typumwandlung ..." << endl;

			// typumwandlung
			long outputsize2 = outputstreamVelo.str().size();
			long outputsize2_r = outputstreamVelo_r.str().size();
			long outputsize2_t1 = outputstreamVelo_t1.str().size();
			long outputsize2_t2 = outputstreamVelo_t2.str().size();
			//cout << "rank: " << rank << "; step: " << simstep << "; outputsize: " << outputsize << endl;
			char output2[outputsize2+1];
			char output2_r[outputsize2_r+1];
			char output2_t1[outputsize2_t1+1];
			char output2_t2[outputsize2_t2+1];
			strcpy(output2,outputstreamVelo.str().c_str());
			strcpy(output2_r,outputstreamVelo_r.str().c_str());
			strcpy(output2_t1,outputstreamVelo_t1.str().c_str());
			strcpy(output2_t2,outputstreamVelo_t2.str().c_str());

			global_log->info() << "typumwandlung ..." << endl;
			global_log->info() << "opening files ..." << endl;

			// Datei zum schreiben öffnen, daten schreiben
			ofstream fileout2(filenameVelo, ios::out|ios::app);
			fileout2 << output2;
			fileout2.close();

			ofstream fileout2_r(filenameVelo_r, ios::out|ios::app);
			fileout2_r << output2_r;
			fileout2_r.close();

			ofstream fileout2_t1(filenameVelo_t1, ios::out|ios::app);
			fileout2_t1 << output2_t1;
			fileout2_t1.close();

			ofstream fileout2_t2(filenameVelo_t2, ios::out|ios::app);
			fileout2_t2 << output2_t2;
			fileout2_t2.close();

			global_log->info() << "files cclosed." << endl;

			// radial velocity distribution

	#ifdef ENABLE_MPI
		}
	#endif

}

void RadialProfilesWriter::finishOutput( ParticleContainer* particleContainer,
                                         DomainDecompBase* domainDecomp, Domain* domain)
{
}

void RadialProfilesWriter::CalcRadialProfiles( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep )
{
	// calculate radius vector
	double dRadiusLength;
	double dShellWidth;
	unsigned int nRadiusIndex;

	double dVelocity;
	double dVelocity_r;
	double dVelocity_t1;
	double dVelocity_t2;

	unsigned int nVelocityIndex;
	unsigned int nVelocityIndex_r;
	unsigned int nVelocityIndex_t1;
	unsigned int nVelocityIndex_t2;

	unsigned int nIndexMax = _nNumShells-1;
	unsigned int nIndexMaxVelo = _nNumDiscretSteps-1;
	double dSystemMidpoint[3];
	double dRadiusVector[3];
	double dMaxVelo = _dVeloMax * 1.1;  // velocity discretisation, highest value with safety factor


//	global_log->info() << "dMaxVelo: " << dMaxVelo << endl;


	// init local density profile array
	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		_nNumMoleculesInsideShellLocal[i] = 0;
	#ifdef ENABLE_MPI
		_nNumMoleculesInsideShellGlobal[i] = 0;
	#endif
	}
	// radial density profile


/*
	if (_writeFreqAverage % simstep == 0)
	{
		// reset local radial velocity profile array
		for(int i = 0; i < _nNumShells; i++)
		{
			for(int j = 0; j < _nNumDiscretSteps; j++)
			{
				_veloDistrMatrixLocal[i][j] = 0;
			}
		}
		// radial velocity distribution
	}
*/

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
		dVelocity = sqrt( pos->v2() );
		global_log->info() << "dVelocity: " << dVelocity << endl;

		Vector3d vec3dMidpoint( dSystemMidpoint[0], dSystemMidpoint[1], dSystemMidpoint[2] );
		Vector3d vec3dPosition( pos->r(0), pos->r(1), pos->r(2) );
		Vector3d vec3dVelocity( pos->v(0), pos->v(1), pos->v(2) );

		// in Komponenten zerlegen
		vec3dVelocity.CalcRadialTangentialComponentLength( vec3dMidpoint, vec3dPosition, dVelocity_r, dVelocity_t1, dVelocity_t2);

		nVelocityIndex    = (unsigned int) (dVelocity    / dMaxVelo * _nNumDiscretSteps);
		nVelocityIndex_r  = (unsigned int) (dVelocity_r  / dMaxVelo * _nNumDiscretSteps);
		nVelocityIndex_t1 = (unsigned int) (dVelocity_t1 / dMaxVelo * _nNumDiscretSteps);
		nVelocityIndex_t2 = (unsigned int) (dVelocity_t2 / dMaxVelo * _nNumDiscretSteps);

//		global_log->info() << "nVelocityIndex: " << nVelocityIndex << endl;
//		global_log->info() << "nVelocityIndex_r: " << nVelocityIndex_r << endl;
//		global_log->info() << "nVelocityIndex_t1: " << nVelocityIndex_t1 << endl;
//		global_log->info() << "nVelocityIndex_t2: " << nVelocityIndex_t2 << endl;

		if(nVelocityIndex < nIndexMaxVelo && nRadiusIndex < nIndexMax)  // respect finite resolution of radius and velocity
		{
			_veloDistrMatrixLocal[nRadiusIndex][nVelocityIndex]++;
		}
		if(nVelocityIndex_r < nIndexMaxVelo && nRadiusIndex < nIndexMax)  // respect finite resolution of radius and velocity
		{
			_veloDistrMatrixLocal_r[nRadiusIndex][nVelocityIndex_r]++;
		}
		if(nVelocityIndex_t1 < nIndexMaxVelo && nRadiusIndex < nIndexMax)  // respect finite resolution of radius and velocity
		{
			_veloDistrMatrixLocal_t1[nRadiusIndex][nVelocityIndex_t1]++;
		}
		if(nVelocityIndex_t2 < nIndexMaxVelo && nRadiusIndex < nIndexMax)  // respect finite resolution of radius and velocity
		{
			_veloDistrMatrixLocal_t2[nRadiusIndex][nVelocityIndex_t2]++;
		}
		// radial velocity distribution

//		global_log->info() << "calculation of velocity distribution done." << endl;

	}  // foreach molecule

	#ifdef ENABLE_MPI

		// radial density profile

		//MPI_Reduce( dDensityProfileLocal, _dDensityProfile, _nNumShells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _nNumMoleculesInsideShellLocal, _nNumMoleculesInsideShellGlobal, _nNumShells, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		for(unsigned int i = 0; i < _nNumShells; i++)
		{
			_dDensityProfile[i] = (_dDensityProfile[i] + _nNumMoleculesInsideShellGlobal[i] / _dShellVolumesReduced[i] ) * _dAverageFactor;  // _dAverageFactor == 1.0 or 0.5
		}
		// radial density profile


		global_log->info() << "MPI - reducing..." << endl;


		// radial velocity distribution

		for(unsigned int i = 0; i < _nNumShells; i++)
		{
			MPI_Reduce( _veloDistrMatrixLocal[i],    _veloDistrMatrixGlobal[i],    _nNumDiscretSteps, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
			// global_log->info() << "i = " << i << ", ges - reduction done." << endl;
		}

		for(unsigned int i = 0; i < _nNumShells; i++)
		{
			MPI_Reduce( _veloDistrMatrixLocal_r[i],  _veloDistrMatrixGlobal_r[i],  _nNumDiscretSteps, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
			// global_log->info() << "i = " << i << ", r - reduction done." << endl;
		}

		for(unsigned int i = 0; i < _nNumShells; i++)
		{
			MPI_Reduce( _veloDistrMatrixLocal_t1[i], _veloDistrMatrixGlobal_t1[i], _nNumDiscretSteps, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
			global_log->info() << "i = " << i << ", t1 - reduction done." << endl;
		}

		for(unsigned int i = 0; i < _nNumShells; i++)
		{
			MPI_Reduce( _veloDistrMatrixLocal_t2[i], _veloDistrMatrixGlobal_t2[i], _nNumDiscretSteps, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
			global_log->info() << "i = " << i << ", t2 - reduction done." << endl;
		}
		// radial velocity distribution

		global_log->info() << "MPI - reduction done." << endl;

	#else

		// radial density profile
		for(int i = 0; i < _nNumShells; i++)
		{
			_dDensityProfile[i] = (_dDensityProfile[i] + _nNumMoleculesInsideShellLocal[i] / _dShellVolumesReduced[i] ) * _dAverageFactor;  // _dAverageFactor == 1.0 or 0.5
		}
		// radial density profile

	#endif

}

double RadialProfilesWriter::CalcMaxVelocity( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain )
{
	double dVelo, dVeloMaxLocal, dVeloMaxGlobal;
	dVeloMaxLocal = 0;
	dVeloMaxGlobal = 0;

	for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next() )
	{
		dVelo = sqrt( pos->v2() );
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

	double dVeloMax = _dVeloMax * 1.1;  // velocity discretisation, highest value with safety factor
	double dRadiusMax = domain->getGlobalLength(0) * 0.5;
	double dShellWidth = dRadiusMax / _nNumShells;

	double dRi, dRa;


	// init local radial velocity distribution matrix
	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		for(unsigned int j = 0; j < _nNumDiscretSteps; j++)
		{
			_veloDistrMatrixLocal[i][j]    = 0;
			_veloDistrMatrixLocal_r[i][j]  = 0;
			_veloDistrMatrixLocal_t1[i][j] = 0;
			_veloDistrMatrixLocal_t2[i][j] = 0;
		}
	}


	// init local molecules array / density profile
	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		_dDensityProfile[i] = 0.0;
		_nNumMoleculesInsideShellLocal[i] = 0;

#ifdef ENABLE_MPI
		_nNumMoleculesInsideShellGlobal[i] = 0;
#endif
	}


	// calc discrete radius values
	_dDiscreteRadiusValues = new double[_nNumShells];

	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		_dDiscreteRadiusValues[i] = (i + 0.5) * dShellWidth;
	}


	// calc reduced shell volumes
	_dShellVolumesReduced = new double[_nNumShells];

	for(unsigned int i = 0; i < _nNumShells; i++)
	{
		dRi = dShellWidth * i;
		dRa = dShellWidth * (i+1);
		_dShellVolumesReduced[i] = pi * 4. / 3. * (dRa*dRa*dRa - dRi*dRi*dRi);  // V_shell = 4/3 pi (Ra^3 - Ri^3)
	}


	// calc discrete velocity values
	_dDiscreteVelocityValues = new double[_nNumDiscretSteps];
	double dDeltaVelo = dVeloMax / _nNumDiscretSteps;

	for(unsigned int i = 0; i < _nNumDiscretSteps; i++)
	{
		_dDiscreteVelocityValues[i] = dDeltaVelo * (i + 0.5);
	}

	_bDiscretisationDone = true;
}











