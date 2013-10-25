#ifndef RADIALPROFILESWRITER_H_
#define RADIALPROFILESWRITER_H_

#include "io/OutputBase.h"
#include "io/xyVal.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

class RadialProfilesWriter : public OutputBase
{

public:
	//! @brief writes .rpf files
	//!
	RadialProfilesWriter() {}
	RadialProfilesWriter( unsigned long writeFrequency, unsigned long updateFrequency, unsigned long writeFreqAverage,
			              unsigned int nNumShells, unsigned int nNumDiscretSteps, unsigned long nNumTimestepsMaxValDetection,
                          bool incremental );
	~RadialProfilesWriter();
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
	
	std::string getPluginName() { return std::string("RadialProfilesWriter"); }

	void CalcRadialProfiles( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain );
	double CalcMaxVelocity( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain );
	void DoDiscretisation(Domain* domain);

private:
	unsigned long _writeFrequency;
	unsigned long _updateFrequency;
	unsigned long _writeFreqAverage;
	unsigned int  _nNumShells;
	unsigned int  _nNumDiscretSteps;  // resolution of velocity value

	// first determine max value for discretisation
	unsigned long _nNumDetectionTimesteps;
	unsigned long _nNumTimestepsDetected;
	double _dVeloMax;
	double* _dDiscreteRadiusValues;
	double* _dDiscreteVelocityValues;
	double* _dDensityProfile;
	double* _dShellVolumesReduced;

	double _dAverageFactor; // == 0.5 or 1.0

	bool _incremental;
	bool _bMaxVeloDetermined;
	bool _bDiscretisationDone;
	bool _bWroteHeaderVmax;

	// radial density profile
	unsigned long* _nNumMoleculesInsideShellLocal;
#ifdef ENABLE_MPI
	unsigned long* _nNumMoleculesInsideShellGlobal;
#endif

	unsigned long** _veloDistrMatrix;
	unsigned long** _veloDistrMatrixAve;

	bool _appendTimestamp;

};  // class RadialProfilesWriter : public OutputBase

#endif /*RADIALPROFILESWRITER_H_*/
