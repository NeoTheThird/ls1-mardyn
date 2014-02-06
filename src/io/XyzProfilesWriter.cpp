// XyzProfilesWriter.cpp

#include "io/XyzProfilesWriter.h"
// #include "io/xyVal.h"  // <-- not used any more
#include "Common.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include <list>
#include <cmath>

#ifdef ENABLE_MPI
#include <mpi.h>
#include <sys/resource.h>
#endif

using namespace std;

XyzProfilesWriter::XyzProfilesWriter( unsigned long writeFrequency, unsigned long updateFrequency, unsigned long nNumAverageTimesteps,
		                              unsigned int* nNumShells, unsigned int* nMoveSteps, Domain* ptrDomain)
{
	// pointer to domain
	_ptrDomain = ptrDomain;

	unsigned int nNumComponents;
	nNumComponents = _ptrDomain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

	_writeFrequency  = writeFrequency;
	_updateFrequency = updateFrequency;

    // number of shells that divide the domain
	_nNumShellsX = nNumShells[0];
	_nNumShellsY = nNumShells[1];
	_nNumShellsZ = nNumShells[2];

    // moving shell, steps to next shell --> increasing datapoints
	_nMoveStepsX = nMoveSteps[0];
	_nMoveStepsY = nMoveSteps[1];
	_nMoveStepsZ = nMoveSteps[2];

	// shell width
	_dShellWidth = new double[3];  // 3 dimensions: x, y, z
	_dShellWidth[0] = ptrDomain->getGlobalLength(0) / _nNumShellsX;
	_dShellWidth[1] = ptrDomain->getGlobalLength(1) / _nNumShellsY;
	_dShellWidth[2] = ptrDomain->getGlobalLength(2) / _nNumShellsZ;

	// shell volumes
	_dShellVolumeX = 0.0;
	_dShellVolumeY = 0.0;
	_dShellVolumeZ = 0.0;

	CalcShellVolumes(_ptrDomain);

	// number of coordinates
	_nNumMidpointPositions = new unsigned int[3];  // 3 dimensions: x, y, z
	_nNumMidpointPositions[0] = _nMoveStepsX * (_nNumShellsX-1) + 1;
	_nNumMidpointPositions[1] = _nMoveStepsY * (_nNumShellsY-1) + 1;
	_nNumMidpointPositions[2] = _nMoveStepsZ * (_nNumShellsZ-1) + 1;

	// allocate number of molecules data structure
	_nNumMoleculesInsideShells = new unsigned long**[3];  // 3 dimensions: x, y, z

	// number of components
	for(unsigned long d=0; d<3; d++)
	{
		_nNumMoleculesInsideShells[d] = new unsigned long*[nNumComponents];
	}

	// number of coordinates
	for(unsigned long d=0; d<3; d++)
	{
		for(unsigned long c=0; c < nNumComponents; c++)
		{
			_nNumMoleculesInsideShells[d][c] = new unsigned long[ _nNumMidpointPositions[d] ];
		}
	}

	// allocate mean squared velocity data structure
	_dKineticEnergyInsideShells = new long double**[3];  // 3 dimensions: x, y, z

	// number of components
	for(unsigned long d=0; d<3; d++)
	{
		_dKineticEnergyInsideShells[d] = new long double*[nNumComponents];
	}

	// number of coordinates
	for(unsigned long d=0; d<3; d++)
	{
		for(unsigned long c=0; c < nNumComponents; c++)
		{
			_dKineticEnergyInsideShells[d][c] = new long double[ _nNumMidpointPositions[d] ];
		}
	}

	// values of coordinates
	_dShellMidpointPositionsX = new double[ _nNumMidpointPositions[0] ];
	_dShellMidpointPositionsY = new double[ _nNumMidpointPositions[1] ];
	_dShellMidpointPositionsZ = new double[ _nNumMidpointPositions[2] ];

	CalcShellMidpointPositions(_ptrDomain);  // coordinates

	// density profiles
	_densityProfileX = new double[ _nNumMidpointPositions[0] ];
	_densityProfileY = new double[ _nNumMidpointPositions[1] ];
	_densityProfileZ = new double[ _nNumMidpointPositions[2] ];

	_densityProfileXAvg = new double[ _nNumMidpointPositions[0] ];
	_densityProfileYAvg = new double[ _nNumMidpointPositions[1] ];
	_densityProfileZAvg = new double[ _nNumMidpointPositions[2] ];

	_nNumAverageTimesteps = nNumAverageTimesteps;
	_nAveragedTimestepsDensity = new unsigned long[3];

	// mean values, standard deviations
	_dDensityMean           = new double[3];
	_dDensityMeanAvg        = new double[3];
	_dDensityStandardDev    = new double[3];
	_dDensityStandardDevAvg = new double[3];

	for(unsigned long d=0; d < 3; d++)
	{
		// will be set to zero again, when max number of average timesteps is reached
		_nAveragedTimestepsDensity[d] = 0;

		// mean values, standard deviations
		_dDensityMean[d]           = 0.0;
		_dDensityMeanAvg[d]        = 0.0;
		_dDensityStandardDev[d]    = 0.0;
		_dDensityStandardDevAvg[d] = 0.0;
	}

	// concentration profiles
	_concentrationProfileX = new double*[nNumComponents];
	_concentrationProfileY = new double*[nNumComponents];
	_concentrationProfileZ = new double*[nNumComponents];

	for(unsigned long c=0; c < nNumComponents; c++)
	{
		_concentrationProfileX[c] = new double[ _nNumMidpointPositions[0] ];
		_concentrationProfileY[c] = new double[ _nNumMidpointPositions[1] ];
		_concentrationProfileZ[c] = new double[ _nNumMidpointPositions[2] ];
	}

	_concentrationProfileXAvg = new double*[nNumComponents];
	_concentrationProfileYAvg = new double*[nNumComponents];
	_concentrationProfileZAvg = new double*[nNumComponents];

	for(unsigned long c=0; c < nNumComponents; c++)
	{
		_concentrationProfileXAvg[c] = new double[ _nNumMidpointPositions[0] ];
		_concentrationProfileYAvg[c] = new double[ _nNumMidpointPositions[1] ];
		_concentrationProfileZAvg[c] = new double[ _nNumMidpointPositions[2] ];
	}

	// will be set to zero again, when max number of average timesteps is reached
	_nAveragedTimestepsConcentration = new unsigned long*[3];

	for(unsigned long d=0; d<3; d++)
	{
		_nAveragedTimestepsConcentration[d] = new unsigned long[nNumComponents];
	}

	for(unsigned long d=0; d < 3; d++)
	{
		for(unsigned long c=0; c < nNumComponents; c++)
		{
			_nAveragedTimestepsConcentration[d][c] = 0;
		}
	}

	// temperature profiles
	_temperatureProfileX = new long double*[nNumComponents];
	_temperatureProfileY = new long double*[nNumComponents];
	_temperatureProfileZ = new long double*[nNumComponents];

	_temperatureProfileAvgX = new long double*[nNumComponents];
	_temperatureProfileAvgY = new long double*[nNumComponents];
	_temperatureProfileAvgZ = new long double*[nNumComponents];

	for(unsigned long c=0; c < nNumComponents; c++)
	{
		_temperatureProfileX[c] = new long double[ _nNumMidpointPositions[0] ];
		_temperatureProfileY[c] = new long double[ _nNumMidpointPositions[1] ];
		_temperatureProfileZ[c] = new long double[ _nNumMidpointPositions[2] ];

		_temperatureProfileAvgX[c] = new long double[ _nNumMidpointPositions[0] ];
		_temperatureProfileAvgY[c] = new long double[ _nNumMidpointPositions[1] ];
		_temperatureProfileAvgZ[c] = new long double[ _nNumMidpointPositions[2] ];
	}

	_nAveragedTimestepsTemperature = 0;  // will be set to zero again, when max number of average timesteps is reached

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
	unsigned int nNumComponents;
	nNumComponents = _ptrDomain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

	// update profiles with respect to update frequency
	if ( (simstep % _updateFrequency == 0) )
	{
		// count only for dimensions of interest
		if(_nNumShellsX > 0)
		{
			CountMoleculesInsideShells( particleContainer, domainDecomp, domain, XYZPFD_DIMENSION_X);
		}

		if(_nNumShellsY > 0)
		{
			CountMoleculesInsideShells( particleContainer, domainDecomp, domain, XYZPFD_DIMENSION_Y);
		}

		if(_nNumShellsZ > 0)
		{
			CountMoleculesInsideShells( particleContainer, domainDecomp, domain, XYZPFD_DIMENSION_Z);
		}

		// calc density, concentration, temperature profiles
		CalcDensityProfiles();
		CalcConcentrationProfiles();
		CalcTemperatureProfiles();
	}


	global_log->info() << "calculation done." << flush << endl;


	// write out profiles with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;


	// writing .prf-files
	std::stringstream outputstream_rhoX, outputstream_rhoY, outputstream_rhoZ;
	std::stringstream outputstreamC, outputstreamT;
	std::stringstream filenamestream_rhoX, filenamestream_rhoY, filenamestream_rhoZ;
	std::stringstream filenamestreamC, filenamestreamT;

	filenamestream_rhoX << "density-profile_TS-" << simstep << ".prfx";
	filenamestream_rhoY << "density-profile_TS-" << simstep << ".prfy";
	filenamestream_rhoZ << "density-profile_TS-" << simstep << ".prfz";

	filenamestreamC << "concentration-profile_TS-" << simstep << ".prf";
	filenamestreamT << "temperature-profile_TS-" << simstep << ".prf";

	char filename_rhoX[filenamestream_rhoX.str().size()+1];
	char filename_rhoY[filenamestream_rhoY.str().size()+1];
	char filename_rhoZ[filenamestream_rhoZ.str().size()+1];

	char filenameC[filenamestreamC.str().size()+1];
	char filenameT[filenamestreamT.str().size()+1];

	strcpy(filename_rhoX, filenamestream_rhoX.str().c_str());
	strcpy(filename_rhoY, filenamestream_rhoY.str().c_str());
	strcpy(filename_rhoZ, filenamestream_rhoZ.str().c_str());

	strcpy(filenameC,filenamestreamC.str().c_str());
    strcpy(filenameT,filenamestreamT.str().c_str());

	#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		// int numprocs = domainDecomp->getNumProcs();
		if (rank== 0)
		{
	#endif
			outputstream_rhoX << std::setw(36) <<  "           x        rhoX        avgX";
			outputstream_rhoX << std::setw(48) <<  "      rhoX_m      avgX_m      stdDev  stdDev_avg" << endl;
			outputstream_rhoY << std::setw(36) <<  "           y        rhoY        avgY";
			outputstream_rhoY << std::setw(48) <<  "      rhoY_m      avgY_m      stdDev  stdDev_avg" << endl;
			outputstream_rhoZ << std::setw(36) <<  "           z        rhoZ        avgZ";
			outputstream_rhoZ << std::setw(48) <<  "      rhoZ_m      avgZ_m      stdDev  stdDev_avg" << endl;

			outputstreamC << "x                   ";
			for(unsigned long c=1; c < nNumComponents; c++)
			{
				outputstreamC << "c" << c << "                  ";
			}
			outputstreamC << "\n";

            outputstreamT <<  "x                   t                   avg                 \n";


			// concentration, temperature profile
			for(unsigned long p=0; p < _nNumMidpointPositions[0]; p++)
			{
				// concentration profile
				outputstreamC << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsX[p];
				for(unsigned long c=1; c < nNumComponents; c++)
				{
					outputstreamC << std::setw(12) << std::setprecision(6) << _concentrationProfileX[c][p];
				}
				outputstreamC << endl;

				// temperature profile
				outputstreamT << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsX[p];
				outputstreamT << std::setw(12) << std::setprecision(6) << _temperatureProfileX[0][p];
				outputstreamT << std::setw(12) << std::setprecision(6) << _temperatureProfileAvgX[0][p] << endl;
			}


			// density profile, x
			if(_nNumShellsX > 0)
			{
				// first line with mean values
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsX[0];
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _densityProfileX[0];
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _densityProfileXAvg[0];
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _dDensityMean[0];
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _dDensityMeanAvg[0];
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _dDensityStandardDev[0];
				outputstream_rhoX << std::setw(12) << std::setprecision(6) << _dDensityStandardDevAvg[0] << endl;

				for(unsigned long p=1; p < _nNumMidpointPositions[0]; p++)
				{
					// density profile
					outputstream_rhoX << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsX[p];
					outputstream_rhoX << std::setw(12) << std::setprecision(6) << _densityProfileX[p];
					outputstream_rhoX << std::setw(12) << std::setprecision(6) << _densityProfileXAvg[p] << endl;
				}
			}

			// density profile, y
			if(_nNumShellsY > 0)
			{
				// first line with mean values
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsY[0];
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _densityProfileY[0];
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _densityProfileYAvg[0];
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _dDensityMean[1];
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _dDensityMeanAvg[1];
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _dDensityStandardDev[1];
				outputstream_rhoY << std::setw(12) << std::setprecision(6) << _dDensityStandardDevAvg[1] << endl;

				for(unsigned long p=1; p < _nNumMidpointPositions[1]; p++)
				{
					// density profile
					outputstream_rhoY << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsY[p];
					outputstream_rhoY << std::setw(12) << std::setprecision(6) << _densityProfileY[p];
					outputstream_rhoY << std::setw(12) << std::setprecision(6) << _densityProfileYAvg[p] << endl;
				}
			}

			// density profile, z
			if(_nNumShellsZ > 0)
			{
				// first line with mean values
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsZ[0];
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _densityProfileZ[0];
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _densityProfileZAvg[0];
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _dDensityMean[2];
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _dDensityMeanAvg[2];
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _dDensityStandardDev[2];
				outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _dDensityStandardDevAvg[2] << endl;

				for(unsigned long p=1; p < _nNumMidpointPositions[2]; p++)
				{
					// density profile
					outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _dShellMidpointPositionsZ[p];
					outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _densityProfileZ[p];
					outputstream_rhoZ << std::setw(12) << std::setprecision(6) << _densityProfileZAvg[p] << endl;
				}
			}


			// typumwandlung
			long outputsize_rhoX = outputstream_rhoX.str().size();
			long outputsize_rhoY = outputstream_rhoY.str().size();
			long outputsize_rhoZ = outputstream_rhoZ.str().size();

			long outputsizeC = outputstreamC.str().size();
            long outputsizeT = outputstreamT.str().size();

			//cout << "rank: " << rank << "; step: " << simstep << "; outputsize: " << outputsize << endl;
			char output_rhoX[outputsize_rhoX+1];
			char output_rhoY[outputsize_rhoY+1];
			char output_rhoZ[outputsize_rhoZ+1];

			char outputC[outputsizeC+1];
            char outputT[outputsizeT+1];

			strcpy(output_rhoX,outputstream_rhoX.str().c_str());
			strcpy(output_rhoY,outputstream_rhoY.str().c_str());
			strcpy(output_rhoZ,outputstream_rhoZ.str().c_str());

			strcpy(outputC,outputstreamC.str().c_str());
            strcpy(outputT,outputstreamT.str().c_str());

			// open files to stream data to
            if(_nNumShellsX > 0)
            {
            	ofstream fileout_rhoX(filename_rhoX, ios::out);
            	fileout_rhoX << output_rhoX;
            	fileout_rhoX.close();
            }
            if(_nNumShellsY > 0)
            {
            	ofstream fileout_rhoY(filename_rhoY, ios::out);
            	fileout_rhoY << output_rhoY;
            	fileout_rhoY.close();
            }
            if(_nNumShellsZ > 0)
            {
            	ofstream fileout_rhoZ(filename_rhoZ, ios::out);
            	fileout_rhoZ << output_rhoZ;
            	fileout_rhoZ.close();
            }

            // open files
			ofstream fileoutC(filenameC, ios::out);
            ofstream fileoutT(filenameT, ios::out);

            // stream data to files
			fileoutC << outputC;
            fileoutT << outputT;

            // close files
			fileoutC.close();
            fileoutT.close();
	#ifdef ENABLE_MPI
		}
	#endif
}

void XyzProfilesWriter::finishOutput( ParticleContainer* particleContainer,
                                      DomainDecompBase* domainDecomp, Domain* domain)
{
}

/*
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

	// Anzahl der Shells
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
*/

void XyzProfilesWriter::CountMoleculesInsideShells( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
		                                            int nDimension )
{
	list<Molecule*> particlePtrs;
	list<Molecule*>::iterator it;
	double dLowerCorner[3];
	double dUpperCorner[3];
	double dDelta[3];
	double dBoxLength[3];
	unsigned int nNumComponents = domain->getNumberOfComponents() + 1;  // component: 0 stands for all molecules --> +1
	unsigned long nNumMoleculesInsideShell;
	unsigned long* nNumMoleculesInsideShellComponentwise;
	long double* dKineticEnergyInsideShellComponentwise;
	long double dKineticEnergyInsideShell;
	nNumMoleculesInsideShellComponentwise = new unsigned long[nNumComponents];
	dKineticEnergyInsideShellComponentwise = new long double[nNumComponents];

    // variables to calculate kinetic energy
	double mv2;
	double Iw2;
	double dKineticEnergy;
	
	dDelta[0] = _dShellWidth[0] / _nMoveStepsX;
	dDelta[1] = _dShellWidth[1] / _nMoveStepsY;
	dDelta[2] = _dShellWidth[2] / _nMoveStepsZ;
	
	dBoxLength[0] = _ptrDomain->getGlobalLength(0);
	dBoxLength[1] = _ptrDomain->getGlobalLength(1);
	dBoxLength[2] = _ptrDomain->getGlobalLength(2);

	// clean number of molecules matrix
	for(unsigned long c=0; c < nNumComponents; c++)
	{
		for(unsigned long p=0; p < _nNumMidpointPositions[nDimension]; p++)
		{
			_nNumMoleculesInsideShells[nDimension][c][p] = 0;
			_dKineticEnergyInsideShells[nDimension][c][p] = 0.0;
		}
	}


	dLowerCorner[0] = 0.0;
	dLowerCorner[1] = 0.0;
	dLowerCorner[2] = 0.0;
	dUpperCorner[0] = dBoxLength[0];
	dUpperCorner[1] = dBoxLength[1];
	dUpperCorner[2] = dBoxLength[2];

	dUpperCorner[nDimension] = _dShellWidth[nDimension];  // only one dimension divided into slices

	//global_log->info() << "counting molecules inside shells..." << endl;

	for(unsigned long p=0; p < _nNumMidpointPositions[nDimension]; p++)
	{
		particlePtrs.clear();
		particleContainer->getRegion(dLowerCorner, dUpperCorner, particlePtrs);

		// global_log->info() << "particlePtrs - size: " << particlePtrs.size() << endl;

		// clean count variables
		nNumMoleculesInsideShell = 0;
		dKineticEnergyInsideShell = 0;
		for(unsigned long c=0; c < nNumComponents; c++)
		{
			nNumMoleculesInsideShellComponentwise[c] = 0;
			dKineticEnergyInsideShellComponentwise[c] = 0.0;
		}

		// global_log->info() << "counting molecules inside shell " << p << " ..." << endl;

		for(it = particlePtrs.begin(); it != particlePtrs.end(); ++it)
		{
			nNumMoleculesInsideShell++;
			unsigned int cid = (*it)->componentid() + 1;  // why +1? --> componentid starts at 0.

			// global_log->info() << "cid: " << cid << endl;
			// global_log->info() << "v2: " << (*it)->v2() << endl;			

			// record _twice_ the total (ordered + unordered) kinetic energy
			mv2 = 0.0;
			Iw2 = 0.0;
			dKineticEnergy = 0.0;
			(*it)->calculate_mv2_Iw2(mv2, Iw2);
			dKineticEnergy += mv2+Iw2;

			// calculate kinetic energy sum
//			dKineticEnergyInsideShellComponentwise[0] = (dKineticEnergyInsideShellComponentwise[0] * nNumMoleculesInsideShellComponentwise[0] + dKineticEnergy ) / (nNumMoleculesInsideShellComponentwise[0] + 1);
//			dKineticEnergyInsideShellComponentwise[cid] = (dKineticEnergyInsideShellComponentwise[cid] * nNumMoleculesInsideShellComponentwise[cid] + dKineticEnergy ) / (nNumMoleculesInsideShellComponentwise[cid] + 1);

			dKineticEnergyInsideShellComponentwise[0]   += dKineticEnergy;
			dKineticEnergyInsideShellComponentwise[cid] += dKineticEnergy;

			// count molecules
			nNumMoleculesInsideShellComponentwise[0]++;
			nNumMoleculesInsideShellComponentwise[cid]++;
		}

		// reduce data from all processes, for every component
		for(unsigned long c=0; c < nNumComponents; c++)
		{
			// reduce number of molecules
			nNumMoleculesInsideShell = nNumMoleculesInsideShellComponentwise[c];
			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendUnsLong(nNumMoleculesInsideShell);
			domainDecomp->collCommAllreduceSum();
			nNumMoleculesInsideShell = domainDecomp->collCommGetUnsLong();
			domainDecomp->collCommFinalize();

			// reduce kinetic energy
			dKineticEnergyInsideShell = dKineticEnergyInsideShellComponentwise[c];
			domainDecomp->collCommInit(1);
			domainDecomp->collCommAppendLongDouble(dKineticEnergyInsideShell);
			domainDecomp->collCommAllreduceSum();
			dKineticEnergyInsideShell = domainDecomp->collCommGetLongDouble();
			domainDecomp->collCommFinalize();

			// store reduced values
			_nNumMoleculesInsideShells[nDimension][c][p] = nNumMoleculesInsideShell;
            _dKineticEnergyInsideShells[nDimension][c][p] = dKineticEnergyInsideShell / nNumMoleculesInsideShell;

        }

		// calculate corners of next region (shell)
		dLowerCorner[nDimension] += dDelta[nDimension];
		dUpperCorner[nDimension] += dDelta[nDimension];
	}

	// free memory
	delete[] nNumMoleculesInsideShellComponentwise;
	delete[] dKineticEnergyInsideShellComponentwise;
}


void XyzProfilesWriter::CalcShellVolumes(Domain* ptrDomain)
{
	double dBoxLength[3];
	dBoxLength[0] = ptrDomain->getGlobalLength(0);
	dBoxLength[1] = ptrDomain->getGlobalLength(1);
	dBoxLength[2] = ptrDomain->getGlobalLength(2);

	_dShellVolumeX = dBoxLength[0] / _nNumShellsX * dBoxLength[1] * dBoxLength[2];
	_dShellVolumeY = dBoxLength[1] / _nNumShellsY * dBoxLength[0] * dBoxLength[2];
	_dShellVolumeZ = dBoxLength[2] / _nNumShellsZ * dBoxLength[0] * dBoxLength[1];
}

void XyzProfilesWriter::CalcShellMidpointPositions(Domain* ptrDomain)  // coordinates
{
	double dDelta[3];

	dDelta[0] = _dShellWidth[0] / _nMoveStepsX;
	dDelta[1] = _dShellWidth[1] / _nMoveStepsY;
	dDelta[2] = _dShellWidth[2] / _nMoveStepsZ;

	// coordinates: X
	_dShellMidpointPositionsX[0] = _dShellWidth[0] * 0.5;

	for(unsigned int i=1; i < _nNumMidpointPositions[0]; i++)
	{
		_dShellMidpointPositionsX[i] = _dShellMidpointPositionsX[i-1] + dDelta[0];
	}

	// coordinates: Y
	_dShellMidpointPositionsY[0] = _dShellWidth[1] * 0.5;

	for(unsigned int i=1; i < _nNumMidpointPositions[1]; i++)
	{
		_dShellMidpointPositionsY[i] = _dShellMidpointPositionsY[i-1] + dDelta[1];
	}

	// coordinates: Z
	_dShellMidpointPositionsZ[0] = _dShellWidth[2] * 0.5;

	for(unsigned int i=1; i < _nNumMidpointPositions[2]; i++)
	{
		_dShellMidpointPositionsZ[i] = _dShellMidpointPositionsZ[i-1] + dDelta[2];
	}
}

void XyzProfilesWriter::CalcDensityProfiles()
{
	double rho_old, rho_new, n, dDeltaRho, dDeltaRhoAvg, dSumStdDev, dSumStdDevAvg, numVals;

	// x dimension
	if(_nNumShellsX > 0)
	{
		// mean values
		_dDensityMean[0]    = 0.0;
		_dDensityMeanAvg[0] = 0.0;

		// number of averaged values
		n = (double) _nAveragedTimestepsDensity[0];

		for(unsigned int p=0; p < _nNumMidpointPositions[0]; p++)
		{
			_densityProfileX[p] = _nNumMoleculesInsideShells[0][0][p] / _dShellVolumeX;
			_dDensityMean[0] += _densityProfileX[p];

			// average
			rho_old = _densityProfileXAvg[p];
			rho_new = _densityProfileX[p];

			_densityProfileXAvg[p] = (rho_old * n + rho_new) / (n + 1.0);
			_dDensityMeanAvg[0] += _densityProfileXAvg[p];
		}
		// mean values
		_dDensityMean[0]    = _dDensityMean[0]    / (double)(_nNumMidpointPositions[0]);
		_dDensityMeanAvg[0] = _dDensityMeanAvg[0] / (double)(_nNumMidpointPositions[0]);

		// calculate standard deviation
		// Papula Formelsamlung, S.294: Standardabweichung des Mittelwerts
		dSumStdDev = 0.0;
		dSumStdDevAvg = 0.0;
		numVals = _nNumMidpointPositions[0];

		for(unsigned int p=0; p < _nNumMidpointPositions[0]; p++)
		{
			dDeltaRho = _densityProfileX[p] - _dDensityMean[0];
			dSumStdDev += dDeltaRho * dDeltaRho;

			dDeltaRhoAvg = _densityProfileXAvg[p] - _dDensityMeanAvg[0];
			dSumStdDevAvg += dDeltaRhoAvg * dDeltaRhoAvg;
		}
		_dDensityStandardDev[0]    = sqrt( dSumStdDev    / (numVals * (numVals - 1.0) ) );
		_dDensityStandardDevAvg[0] = sqrt( dSumStdDevAvg / (numVals * (numVals - 1.0) ) );

		// inkrement averaged timesteps counter
		if(_nAveragedTimestepsDensity[0] == _nNumAverageTimesteps - 1)  // => z.B.: (n + 1) = 1000
		{
			_nAveragedTimestepsDensity[0] = 0;
		}
		else
		{
			_nAveragedTimestepsDensity[0]++;
		}
	}


	// y dimension
	if(_nNumShellsY > 0)
	{
		// mean values
		_dDensityMean[1]    = 0.0;
		_dDensityMeanAvg[1] = 0.0;

		// number of averaged values
		n = (double) _nAveragedTimestepsDensity[1];

		for(unsigned int p=0; p < _nNumMidpointPositions[1]; p++)
		{
			_densityProfileY[p] = _nNumMoleculesInsideShells[1][0][p] / _dShellVolumeY;
			_dDensityMean[1] += _densityProfileY[p];

			// average
			rho_old = _densityProfileYAvg[p];
			rho_new = _densityProfileY[p];

			_densityProfileYAvg[p] = (rho_old * n + rho_new) / (n + 1.0);
			_dDensityMeanAvg[1] += _densityProfileYAvg[p];
		}
		// mean values
		_dDensityMean[1]    = _dDensityMean[1]    / (double)(_nNumMidpointPositions[1]);
		_dDensityMeanAvg[1] = _dDensityMeanAvg[1] / (double)(_nNumMidpointPositions[1]);

		// calculate standard deviation
	    // Papula Formelsamlung, S.294: Standardabweichung des Mittelwerts
		dSumStdDev = 0.0;
		dSumStdDevAvg = 0.0;
		numVals = _nNumMidpointPositions[1];

		for(unsigned int p=0; p < _nNumMidpointPositions[1]; p++)
		{
			dDeltaRho = _densityProfileY[p] - _dDensityMean[1];
			dSumStdDev += dDeltaRho * dDeltaRho;

			dDeltaRhoAvg = _densityProfileYAvg[p] - _dDensityMeanAvg[1];
			dSumStdDevAvg += dDeltaRhoAvg * dDeltaRhoAvg;
		}
		_dDensityStandardDev[1]    = sqrt( dSumStdDev    / (numVals * (numVals - 1.0) ) );
		_dDensityStandardDevAvg[1] = sqrt( dSumStdDevAvg / (numVals * (numVals - 1.0) ) );

		// inkrement averaged timesteps counter
		if(_nAveragedTimestepsDensity[1] == _nNumAverageTimesteps)  // => z.B.: (n + 1) = 1000
		{
			_nAveragedTimestepsDensity[1] = 0;
		}
		else
		{
			_nAveragedTimestepsDensity[1]++;
		}
	}


	// z dimension
	if(_nNumShellsZ > 0)
	{
		// mean values
		_dDensityMean[2]    = 0.0;
		_dDensityMeanAvg[2] = 0.0;

		// number of averaged values
		n = (double) _nAveragedTimestepsDensity[2];

		for(unsigned int p=0; p < _nNumMidpointPositions[2]; p++)
		{
			_densityProfileZ[p] = _nNumMoleculesInsideShells[2][0][p] / _dShellVolumeZ;
			_dDensityMean[2] += _densityProfileZ[p];

			// average
			rho_old = _densityProfileZAvg[p];
			rho_new = _densityProfileZ[p];

			_densityProfileZAvg[p] = (rho_old * n + rho_new) / (n + 1.0);
			_dDensityMeanAvg[2] += _densityProfileZAvg[p];
		}
		// mean values
		_dDensityMean[2]    = _dDensityMean[2]    / (double)(_nNumMidpointPositions[2]);
		_dDensityMeanAvg[2] = _dDensityMeanAvg[2] / (double)(_nNumMidpointPositions[2]);

		// calculate standard deviation
		// Papula Formelsamlung, S.294: Standardabweichung des Mittelwerts
		dSumStdDev = 0.0;
		dSumStdDevAvg = 0.0;
		numVals = _nNumMidpointPositions[2];

		for(unsigned int p=0; p < _nNumMidpointPositions[2]; p++)
		{
			dDeltaRho = _densityProfileZ[p] - _dDensityMean[2];
			dSumStdDev += dDeltaRho * dDeltaRho;

			dDeltaRhoAvg = _densityProfileZAvg[p] - _dDensityMeanAvg[2];
			dSumStdDevAvg += dDeltaRhoAvg * dDeltaRhoAvg;
		}
		_dDensityStandardDev[2]    = sqrt( dSumStdDev    / (numVals * (numVals - 1.0) ) );
		_dDensityStandardDevAvg[2] = sqrt( dSumStdDevAvg / (numVals * (numVals - 1.0) ) );

		// inkrement averaged timesteps counter
		if(_nAveragedTimestepsDensity[2] == _nNumAverageTimesteps - 1)  // => z.B.: (n + 1) = 1000
		{
			_nAveragedTimestepsDensity[2] = 0;
		}
		else
		{
			_nAveragedTimestepsDensity[2]++;
		}
	}
}

void XyzProfilesWriter::CalcConcentrationProfiles()
{
	unsigned int nNumComponents;
	nNumComponents = _ptrDomain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components
	double c_old, c_new, n;

	for(unsigned long c=1; c < nNumComponents; c++)
	{
		// x dimension
		if(_nNumShellsX > 0)
		{
			for(unsigned long p=0; p < _nNumMidpointPositions[0]; p++)
			{
				_concentrationProfileX[c][p] = (double)(_nNumMoleculesInsideShells[0][c][p]) / (double)(_nNumMoleculesInsideShells[0][0][p]);

				// average
				c_old = _concentrationProfileXAvg[c][p];
				c_new = _concentrationProfileX[c][p];
				n = (double) _nAveragedTimestepsConcentration[0][c];

				_concentrationProfileXAvg[c][p] = (c_old * n + c_new) / (n + 1.0);

				if(_nAveragedTimestepsConcentration[0][c] == _nNumAverageTimesteps)
				{
					_nAveragedTimestepsConcentration[0][c] = 0;
				}
				else
				{
					_nAveragedTimestepsConcentration[0][c]++;
				}
			}
		}


		// y dimension
		if(_nNumShellsY > 0)
		{
			for(unsigned long p=0; p < _nNumMidpointPositions[0]; p++)
			{
				_concentrationProfileY[c][p] = (double)(_nNumMoleculesInsideShells[1][c][p]) / (double)(_nNumMoleculesInsideShells[0][0][p]);

				// average
				c_old = _concentrationProfileYAvg[c][p];
				c_new = _concentrationProfileY[c][p];
				n = (double) _nAveragedTimestepsConcentration[1][c];

				_concentrationProfileYAvg[c][p] = (c_old * n + c_new) / (n + 1.0);

				if(_nAveragedTimestepsConcentration[1][c] == _nNumAverageTimesteps)
				{
					_nAveragedTimestepsConcentration[1][c] = 0;
				}
				else
				{
					_nAveragedTimestepsConcentration[1][c]++;
				}
			}
		}

		// z dimension
		if(_nNumShellsZ > 0)
		{
			for(unsigned long p=0; p < _nNumMidpointPositions[0]; p++)
			{
				_concentrationProfileZ[c][p] = (double)(_nNumMoleculesInsideShells[2][c][p]) / (double)(_nNumMoleculesInsideShells[0][0][p]);


				// average
				c_old = _concentrationProfileZAvg[c][p];
				c_new = _concentrationProfileZ[c][p];
				n = (double) _nAveragedTimestepsConcentration[2][c];

				_concentrationProfileZAvg[c][p] = (c_old * n + c_new) / (n + 1.0);

				if(_nAveragedTimestepsConcentration[2][c] == _nNumAverageTimesteps)
				{
					_nAveragedTimestepsConcentration[2][c] = 0;
				}
				else
				{
					_nAveragedTimestepsConcentration[2][c]++;
				}

			}
		}
	}
}

void XyzProfilesWriter::CalcTemperatureProfiles()
{
	unsigned int nNumComponents;
	nNumComponents = _ptrDomain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components

	double dInverseDOF = 1 / 3.0;  // only translatoric movement

	// x dimension
	for(unsigned long c=0; c < nNumComponents; c++)
	{
		for(unsigned int p=0; p < _nNumMidpointPositions[0]; p++)
		{
			_temperatureProfileX[c][p] = _dKineticEnergyInsideShells[0][0][p] * dInverseDOF;

			// average
			_temperatureProfileAvgX[c][p] = (_temperatureProfileAvgX[c][p] * _nAveragedTimestepsTemperature + _dKineticEnergyInsideShells[0][0][p] * dInverseDOF) / (_nAveragedTimestepsTemperature + 1);
		}
	}

	if(_nAveragedTimestepsTemperature == _nNumAverageTimesteps - 1)  // => z.B.: (n + 1) = 1000
	{
		_nAveragedTimestepsTemperature = 0;
	}
	else
	{
		_nAveragedTimestepsTemperature++;
	}
}




