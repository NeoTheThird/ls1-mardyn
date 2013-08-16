#ifndef VISWRITER_H_
#define VISWRITER_H_

#include "io/OutputBase.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

class VISWriter : public OutputBase {
public:
	//! @brief writes a .vis_ file
	//!
	// VISWriter(unsigned long writeFrequency, std::string filename, unsigned long numberOfTimesteps,
// bool incremental);
        VISWriter() {}
	VISWriter( unsigned long writeFrequency, std::string outputPrefix, bool incremental);
	~VISWriter();
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
		return std::string("VISWriter");
	}
private:
	std::string _filename;
	unsigned long _writeFrequency;
	unsigned long _numberOfTimesteps;
	bool	_incremental;
	bool	_filenameisdate;
	bool  _wroteVIS;

	std::string _outputPrefix;
	bool _appendTimestamp;
};

#endif /*VISWRITER_H_*/
