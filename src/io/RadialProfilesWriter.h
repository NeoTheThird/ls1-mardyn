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

	void CalcVelocityDistributaion( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain );
	double CalcMaxVelocity( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain );
	void DoDiscretisation();

private:
	unsigned long _writeFrequency;
	unsigned long _updateFrequency;
	unsigned long _writeFreqAverage;
	unsigned int  _nNumShells;
	unsigned int  _nNumDiscretSteps;

	// first determine max value for discretisation
	unsigned long _nNumDetectionTimesteps;
	unsigned long _nNumTimestepsDetected;
	double _dVeloMax;

	bool _incremental;
	bool _bMaxVeloDetermined;
	bool _bDiscretisationDone;
	bool _bWroteHeaderVmax;

	unsigned int** _veloDistrMatrix;
	unsigned int** _veloDistrMatrixAve;

	bool _appendTimestamp;

};  // class RadialProfilesWriter : public OutputBase

#endif /*RADIALPROFILESWRITER_H_*/
