#ifndef CXWRITER_H_
#define CXWRITER_H_

#include "io/OutputBase.h"
#include "io/xyVal.h"
#include "io/Region.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

enum UpdateRegionListStates
{
	URLS_INITIAL_STEP = 0,
	URLS_LEFT_BOUNDARY_SET = 1,
	URLS_RIGHT_BOUNDARY_SET = 2,
	URLS_NUM_VALUES_MIN_REACHED = 3,
	URLS_LAST_ELEMENT = 4,
	URLS_END_OF_DOMAIN_REACHED = 5,
};

enum PresentPhasesStates
{
	PPS_UNKNOWN = 0,
	PPS_COMPONENT_ONE_ONLY = 1,
	PPS_COMPONENT_TWO_ONLY = 2,
	PPS_COMPONENT_ONE_AND_TWO_MIXED = 3,
};

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
	CxWriter( unsigned long writeFrequency, unsigned long updateFrequency, double dDeltaX, unsigned int nDeltaXSteps,
			  bool incremental);
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
	int UpdateRegionList(Domain* domain);
	void ClearRegionList(void);

private:
	unsigned long _writeFrequency;
	unsigned long _updateFrequency;
	bool	_incremental;
	bool	_filenameisdate;

	// zur Bestimmung von Konzentrationsprofil
	unsigned int _nNumMolsSlice;
	unsigned int _nNumMolsSliceComp;
	double _dDeltaX;
	unsigned int _nDeltaXSteps;

	std::vector<xyVal*> _Cx;
	std::vector<Region*> _regionList;

	bool _appendTimestamp;

};  // class CxWriter : public OutputBase

#endif /*CXWRITER_H_*/
