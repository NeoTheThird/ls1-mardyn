#ifndef MMPLDWRITER_H_
#define MMPLDWRITER_H_

#include <string>
#include <vector>
#include <array>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"

enum InitSphereData : uint8_t
{
	ISD_USE_DEFAULT = 1,
	ISD_READ_FROM_FILE = 2
};

class Simulation;
class MmpldWriter : public OutputBase
{
protected:
    MmpldWriter(){};
	//! @brief: writes a mmspd file used by MegaMol
	//!
	//! Depending on write frequency (for example: every timestep, or every 10th, 100th, 1000th ...) number of frames
	//! can be controlled. The *.mmspd-file can be visualized by visualization software like MegaMol.
	//! (for detail information visit: https://svn.vis.uni-stuttgart.de/trac/megamol/)
	//!
	//! @param filename Name of the *.mmspd-file (including path)
	//! @param particleContainer The molecules that have to be written to the file are stored here
	//! @param domainDecomp In the parallel version, the file has to be written by more than one process.
	//!                     Methods to achieve this are available in domainDecomp
	//! @param writeFrequency Controls the frequency of writing out the data (every timestep, every 10th, 100th, ... timestep)
    MmpldWriter(unsigned long writeFrequency, std::string outputPrefix);
	virtual ~MmpldWriter();
	virtual void SetNumSphereTypes() = 0;
	virtual void CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol) = 0;
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex) = 0;

	void InitSphereData();

public:
	void readXML(XMLfileUnits& xmlconfig);

	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu,
			std::map<unsigned, CavityEnsemble>* mcav
	);
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("MmpldWriter");
	}

	void SetInitSphereDataParameters(const uint8_t &bInitSphereData, const std::string &strSphereDataFilename) {
		_bInitSphereData = bInitSphereData; _strSphereDataFilename = strSphereDataFilename;
	}

protected:
	std::string _outputPrefix;
	unsigned long _writeFrequency;
	bool _appendTimestamp;
	std::string _timestampString;
	uint32_t _numSeekEntries;
	uint32_t _frameCount;
#ifdef ENABLE_MPI
	uint64_t *_seekTable;

#endif
	uint8_t  _numComponents;
	uint8_t  _numSitesTotal;
	uint8_t  _numSphereTypes;
	uint8_t* _numSitesPerComp;
	uint8_t* _nCompSitesOffset;
	std::vector<float> _vfSphereRadius;
	std::vector< std::array<uint8_t, 4> > _vaSphereColors;
	std::string _strSphereDataFilename;
	uint8_t _bInitSphereData;
};

class MmpldWriterSimpleSphere : public MmpldWriter
{
public:
	MmpldWriterSimpleSphere(unsigned long writeFrequency, std::string outputPrefix)
			: MmpldWriter(writeFrequency, outputPrefix)
	{
//		MmpldWriter::_numSphereTypes = &(MmpldWriter::_numComponents);
	}
	virtual ~MmpldWriterSimpleSphere() {}

	virtual void SetNumSphereTypes() {_numSphereTypes = _numComponents;}
	virtual void CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol);
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex);
};

class MmpldWriterMultiSphere : public MmpldWriter
{
public:
	MmpldWriterMultiSphere(unsigned long writeFrequency, std::string outputPrefix)
			: MmpldWriter(writeFrequency, outputPrefix)
	{
//		MmpldWriter::_numSphereTypes = &(MmpldWriter::_numSitesTotal);
	}
	virtual ~MmpldWriterMultiSphere() {}

	virtual void SetNumSphereTypes() {_numSphereTypes = _numSitesTotal;}
	virtual void CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol);
	virtual bool GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex);
};

#endif /* MMPLDWRITER_H_ */




