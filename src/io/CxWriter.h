#ifndef CXWRITER_H_
#define CXWRITER_H_

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

class TraWriter;
class CxWriter : public OutputBase
{

friend class TraWriter;  // to access C(x) profile.

public:
	//! @brief writes a .cox file
	//!
    CxWriter() {}
	CxWriter( unsigned long writeFrequency, double DeltaX, bool incremental);
	~CxWriter();
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
		return std::string("CxWriter");
	}

	void CalcConcentrAtX( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned int nCompID);
	void ClearCxVector(void);

private:
	unsigned long _writeFrequency;
	bool	_incremental;
	bool	_filenameisdate;

	// zur Bestimmung von Konzentrationsprofil
	unsigned int _nNumMolsSlice;
	unsigned int _nNumMolsSliceComp;
	double _dDeltaX;
	std::vector<xyVal*> _Cx;

	bool _appendTimestamp;

};  // class CxWriter : public OutputBase

#endif /*CXWRITER_H_*/
