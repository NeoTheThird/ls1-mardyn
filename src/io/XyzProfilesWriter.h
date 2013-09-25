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
	XyzProfilesWriter( unsigned long writeFrequency, unsigned long updateFrequency, double dSliceWidth, unsigned int nMoveSteps,
                       unsigned int* nRecSwitches, bool incremental );
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

	void CalcConcentrProfile( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
			                  unsigned int nCompID, unsigned int nDimension);
	void ClearCxyzVector(unsigned int nDimension);

private:
	unsigned long _writeFrequency;
	unsigned long _updateFrequency;
	bool	_incremental;
	bool	_filenameisdate;

	// zur Bestimmung von Konzentrationsprofil
	unsigned int _nNumMolsSlice;
	unsigned int _nNumMolsSliceComp;
	double _dSliceWidth;
	unsigned int _nMoveSteps;

	std::vector<xyVal*> _Cxyz[3];
	std::vector<xyVal*> _Txyz[3];
	unsigned int _nRecSwitches[6];

	bool _appendTimestamp;

};  // class XyzProfilesWriter : public OutputBase

#endif /*XYZPROFILESWRITER_H_*/
