#ifndef TRAWRITER_H_
#define TRAWRITER_H_

#include "io/OutputBase.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>
#include "io/TraTrack.h"

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

class CxWriter;
class TraWriter : public OutputBase
{
	friend class CxWriter;  // to access box length

public:
	//! @brief writes a .tra file
	//!
    TraWriter() {}
	TraWriter( unsigned long writeFrequency, unsigned long nPhaseBoundaryTrackFreq, unsigned long nNumTrackTimesteps,
			   bool incremental);
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

	void SetBoxLength(double* dBoxLength)
	{
		_dBoxLength[0] = dBoxLength[0];
		_dBoxLength[1] = dBoxLength[1];
		_dBoxLength[2] = dBoxLength[2];
	}

	void UpdateTrackList(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, const unsigned long &simstep);
	int FindPhaseBoundaryMidpointsX( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
                                     std::list<double*> &phaseBoundMidpointList);

private:
	unsigned long _writeFrequency;
	unsigned long _nPhaseBoundaryTrackFreq;
	unsigned long _nNumTrackTimesteps;
	bool	_incremental;
	bool	_filenameisdate;

	std::list<TraTrack*> _TraTrackList;
	double _dBoxLength[3];
	unsigned int _nPhaseBoundID;
	unsigned int _nTrackID;        // number to identify Track
	unsigned int _nNumCompInBox;

	// zur Bestimmung von Konzentrationsprofil
	unsigned int _nNumMolsSlice;
	unsigned int _nNumMolsSliceComp;

	std::vector<double[2]> _concentrAtX;

	bool _appendTimestamp;

};  // class TraWriter : public OutputBase

#endif /*TRAWRITER_H_*/
