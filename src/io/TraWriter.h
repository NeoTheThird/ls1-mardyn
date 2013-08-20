#ifndef TRAWRITER_H_
#define TRAWRITER_H_

#include "io/OutputBase.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

class TraWriter : public OutputBase {
public:
	//! @brief writes a .tra file
	//!
    TraWriter() {}
	TraWriter( unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~TraWriter();
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
		return std::string("TraWriter");
	}

	void SetMidPoint(double* dMidpoint)
	{
		_dMidpoint[0] = dMidpoint[0];
		_dMidpoint[1] = dMidpoint[1];
		_dMidpoint[2] = dMidpoint[2];
	}

	void SetBoxLength(double* dBoxLength)
	{
		_dBoxLength[0] = dBoxLength[0];
		_dBoxLength[1] = dBoxLength[1];
		_dBoxLength[2] = dBoxLength[2];
	}

private:
	std::string _filename;
	unsigned long _writeFrequency;
	unsigned long _numberOfTimesteps;
	bool	_incremental;
	bool	_filenameisdate;
	bool  _wroteTra;

	double _dMidpoint[3];
	double _dBoxLength[3];

	std::string _outputPrefix;
	bool _appendTimestamp;
};

#endif /*TRAWRITER_H_*/
