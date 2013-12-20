#ifndef XYZPROFILESWRITER_H_
#define XYZPROFILESWRITER_H_

#include "io/OutputBase.h"
#include "io/xyVal.h"
#include "io/Region.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

enum XyzProfilesWriterDimensions
{
	XYZPFD_DIMENSION_X = 0,
	XYZPFD_DIMENSION_Y = 1,
	XYZPFD_DIMENSION_Z = 2,
};

class XyzProfilesWriter : public OutputBase
{

public:
	//! @brief writes .prf files
	//!
	XyzProfilesWriter() {}
	XyzProfilesWriter( unsigned long writeFrequency, unsigned long updateFrequency, unsigned long nNumAverageTimesteps,
			           unsigned int* nNumShells, unsigned int* nMoveSteps, Domain* ptrDomain);
	~XyzProfilesWriter();
	//! @todo comment
	
	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu
	);
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("XyzProfilesWriter");
	}

/*
	void CalcConcentrProfile( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
			                  unsigned int nCompID, unsigned int nDimension);
	void ClearCxyzVector(unsigned int nDimension);
*/
	void CountMoleculesInsideShells( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
                                         int nDimension );
	void CalcShellVolumes(Domain* ptrDomain);
	void CalcShellMidpointPositions(Domain* ptrDomain);  // coordinates
	void CalcDensityProfiles();
	void CalcConcentrationProfiles();
        void CalcTemperatureProfiles();

private:
	// frequencys
	unsigned long _writeFrequency;
	unsigned long _updateFrequency;

    // number of shells that divide the domain
	unsigned int _nNumShellsX;
	unsigned int _nNumShellsY;
	unsigned int _nNumShellsZ;

    // moving shell, steps to next shell --> increasing datapoints
	unsigned int _nMoveStepsX;
	unsigned int _nMoveStepsY;
	unsigned int _nMoveStepsZ;

	// shell width
	double* _dShellWidth;

	// shell volumes
	double _dShellVolumeX;
	double _dShellVolumeY;
	double _dShellVolumeZ;

	// number of molecules inside shells, componentwise
	unsigned long*** _nNumMoleculesInsideShells;

	// mean squred velocity inside shells, componentwise
	long double*** _dKineticEnergyInsideShells;

	// number of coordinates
	unsigned int* _nNumMidpointPositions;

	// values of coordinates
	double* _dShellMidpointPositionsX;
	double* _dShellMidpointPositionsY;
	double* _dShellMidpointPositionsZ;

	// density profiles
	double* _densityProfileX;
	double* _densityProfileY;
	double* _densityProfileZ;

	double* _densityProfileXAvg;
	double* _densityProfileYAvg;
	double* _densityProfileZAvg;

	unsigned long _nNumAverageTimesteps;
	unsigned long _nAveragedTimesteps;

	// concentration profiles
	double** _concentrationProfileX;
	double** _concentrationProfileY;
	double** _concentrationProfileZ;

	// temperature profiles
	long double** _temperatureProfileX;
	long double** _temperatureProfileY;
	long double** _temperatureProfileZ;

	long double** _temperatureProfileAvgX;
	long double** _temperatureProfileAvgY;
	long double** _temperatureProfileAvgZ;

	unsigned long _nAveragedTimestepsTemperature;

	// pointer to domain
	Domain* _ptrDomain;

};  // class XyzProfilesWriter : public OutputBase

#endif /*XYZPROFILESWRITER_H_*/
