#include "io/MmpldWriter.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "/usr/include/endian.h"

// set mmpld file version. possible values: 100 or 102
#define MMPLD_FILE_VERSION 100

using Log::global_log;
using namespace std;

MmpldWriter::MmpldWriter(unsigned long writeFrequency, string outputPrefix)
{
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

MmpldWriter::~MmpldWriter(){}

void MmpldWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

//Header Information
void MmpldWriter::initOutput(ParticleContainer* /*particleContainer*/,
                           DomainDecompBase* domainDecomp, Domain* domain)
{
	// number of components / sites
	_numComponents = (uint8_t)domain->getNumberOfComponents();
	_numSitesPerComp  = new uint8_t[_numComponents];
	_nCompSitesOffset = new uint8_t[_numComponents];
	_numSitesTotal = 0;

	vector<Component>& vComponents = *(_simulation.getEnsemble()->getComponents() );
	vector<Component>::iterator cit;
	cit=vComponents.begin();

	_nCompSitesOffset[0] = 0;
	for(uint8_t cid=0; cid<_numComponents; cid++)
	{
		uint8_t numSites = (*cit).numLJcenters();
		_numSitesPerComp[cid] = numSites;
		if(cid>0)
			_nCompSitesOffset[cid] = _nCompSitesOffset[cid-1] + _numSitesPerComp[cid-1];
		// total number of sites (all components)
		_numSitesTotal += numSites;
		cit++;
	}

#ifndef NDEBUG
	for(uint8_t cid=0; cid<_numComponents; cid++)
		cout << "_nCompSitesOffset[" << (uint32_t)cid << "] = " << (uint32_t)_nCompSitesOffset[cid] << endl;
#endif

	// init radius and color of spheres
	this->InitSphereData();
	this->SetNumSphereTypes();

	stringstream filenamestream;
	filenamestream << _outputPrefix;

	if(_appendTimestamp) {
		_timestampString = gettimestring();
		filenamestream << "-" << _timestampString;
	}
	filenamestream << ".mmpld";

	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
#endif
	ofstream mmpldfstream(filename, ios::binary|ios::out);

  //format marker
  uint8_t magicIdentifier[6] = {0x4D, 0x4D, 0x50, 0x4C, 0x44, 0x00};
  mmpldfstream.write((char*)magicIdentifier, sizeof(magicIdentifier));
  
  //version number
  uint16_t versionNumber;
  switch (MMPLD_FILE_VERSION){
	  case 100:
		versionNumber = htole16(100);
		break;

	  case 102:
		versionNumber = htole16(102);
		break;
	  
	  default:
		cout << "Error mmpld-writer: file version " << MMPLD_FILE_VERSION << " not supported." << endl;
		return;
		break;
  }
  
  mmpldfstream.write((char*)&versionNumber, sizeof(versionNumber));
  
  //calculate the number of frames
  unsigned long numTimesteps = _simulation.getNumTimesteps();
  uint32_t numframes;
  uint32_t numframes_le;
  if (_writeFrequency == 0){
	  numframes = numTimesteps;
  }else{
	  numframes = (numTimesteps/_writeFrequency)+1;
  }
  _frameCount = 0;
  _numSeekEntries = numframes+1;
  numframes_le = htole32(numframes);
  mmpldfstream.write((char*)&numframes_le,sizeof(numframes_le));
  #ifdef ENABLE_MPI
	_seekTable = new uint64_t[_numSeekEntries];
  #endif
  //boundary box
  float minbox[3] = {0, 0, 0};
  float maxbox[3];
  for (unsigned short d = 0; d < 3; ++d) maxbox[d] = domain->getGlobalLength(d);
  mmpldfstream.write((char*)&minbox,sizeof(minbox));
  mmpldfstream.write((char*)&maxbox,sizeof(maxbox));
    
  //clipping box
  float inflateRadius = *(_vfSphereRadius.begin() );
  std::vector<float>::iterator it;
  for(it=_vfSphereRadius.begin(); it!=_vfSphereRadius.end(); ++it)
	  if(inflateRadius < (*it) ) inflateRadius = (*it);
  
  for (unsigned short d = 0; d < 3; ++d){
	  maxbox[d] = maxbox[d] + inflateRadius;
	  minbox[d] = minbox[d] - inflateRadius;
  }
  mmpldfstream.write((char*)&minbox,sizeof(minbox));
  mmpldfstream.write((char*)&maxbox,sizeof(maxbox));
  
  //preallocate seektable
  uint64_t seekNum = 0;
  for (uint32_t i = 0; i <= numframes; ++i){
	  mmpldfstream.write((char*)&seekNum,sizeof(seekNum));
  }
  mmpldfstream.close();
#ifdef ENABLE_MPI
	}
#endif

}

void MmpldWriter::doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* /*lmu*/,
		   map<unsigned, CavityEnsemble>* /*mcav*/){
	if (simstep % _writeFrequency == 0) {
		stringstream filenamestream, outputstream;
		filenamestream << _outputPrefix;

		if(_appendTimestamp) {
			filenamestream << "-" << _timestampString;
		}
		filenamestream << ".mmpld";
		
		char filename[filenamestream.str().size()+1];
		strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
		int rank = domainDecomp->getRank();
		int numprocs = domainDecomp->getNumProcs();
		unsigned long numberParticles = particleContainer->getNumberOfParticles();
		long outputsize = 0;

        uint64_t numSpheresPerType[_numSphereTypes];
        for (uint8_t ti = 0; ti < _numSphereTypes; ++ti){
			numSpheresPerType[ti] = 0;
		}

		//calculate number of spheres per component|siteType
		uint32_t molcid = 0;
		for (Molecule* mol = particleContainer->begin(); mol != particleContainer->end(); mol = particleContainer->next())
			this->CalcNumSpheresPerType(numSpheresPerType, mol);
		
		//distribute global component particle count
		uint64_t globalNumCompSpheres[_numSphereTypes];
		if (rank == 0){
			for (uint32_t i = 0; i < _numSphereTypes; ++i){
				globalNumCompSpheres[i] = numSpheresPerType[i];
			}
			MPI_Status status;
			uint64_t numSpheresPerTypeTmp[_numSphereTypes];
			for (int source = rank+1; source < numprocs; ++source){
				int recvcount = sizeof(numSpheresPerTypeTmp);
				int recvtag = 1;
				MPI_Recv(numSpheresPerTypeTmp, recvcount, MPI_BYTE, source, recvtag, MPI_COMM_WORLD, &status);
				for (uint8_t ti = 0; ti < _numSphereTypes; ++ti)
					globalNumCompSpheres[ti] = globalNumCompSpheres[ti] + numSpheresPerTypeTmp[ti];
			}
		}else{
				int dest = 0;
				int sendcount = sizeof(numSpheresPerType);
				int sendtag = 1;
				MPI_Request request;
				MPI_Isend(numSpheresPerType, sendcount, MPI_BYTE, dest, sendtag, MPI_COMM_WORLD, &request);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY|MPI_MODE_APPEND|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
		
		
		//write particle list for each component|site (sphere type)
		for (uint8_t nSphereTypeIndex=0; nSphereTypeIndex<_numSphereTypes; ++nSphereTypeIndex){
			//add space for particle data
			outputsize = (long)numSpheresPerType[nSphereTypeIndex]*12;
			
			//add space for particle list header
			if (rank == 0){
				//add particle list header
				outputsize += 18;
				if (nSphereTypeIndex == 0){
					
					switch (MMPLD_FILE_VERSION){
						case 100:
							//add space for number of particle lists
							outputsize += 4;
							break;
							
						case 102:
							//add space for timestamp and number of particle lists
							outputsize += 8;
							break;
							
						default:
							cout << "Error mmpld-writer: file version " << MMPLD_FILE_VERSION << " not supported." << endl;
							return;
							break;
					}
				}
			}
			
			//send outputsize of current rank to next rank
			for (int dest = rank+1; dest < numprocs; ++dest){
				int sendcount = 1;
				int sendtag = 0;
				MPI_Request request;
				MPI_Isend(&outputsize, sendcount, MPI_LONG, dest, sendtag, MPI_COMM_WORLD, &request);
			}
			//accumulate outputsizes of previous ranks and use it as offset for output file
			MPI_Status status;
			long offset = 0;
			long outputsize_get;
			for (int source = 0; source < rank; ++source){
				int recvcount = 1;
				int recvtag = 0;
				MPI_Recv(&outputsize_get, recvcount, MPI_LONG, source, recvtag, MPI_COMM_WORLD, &status);
				offset += outputsize_get;
			}

			global_log->debug() << "MmpldWriter rank: " << rank << "; step: " << simstep << "; sphereTypeIndex: " << nSphereTypeIndex << "; offset: " << offset << endl;

			MPI_File_seek(fh, offset, MPI_SEEK_END);

			MPI_Barrier(MPI_COMM_WORLD);
			
			//write particle list header
			if (rank == 0){
				
				//write frame header if we are before the first particle list
				if (nSphereTypeIndex == 0){
					//store file position for seek table
					if (_frameCount < _numSeekEntries){
						MPI_Offset entry;
						MPI_File_get_position(fh, &entry);
						_seekTable[_frameCount] = (uint64_t)entry;
					}
					_frameCount = _frameCount + 1;
					
					float frameHeader_timestamp = simstep;
					
					switch (MMPLD_FILE_VERSION){
						case 100:
							//do not write timestamp to frame header
							break;
							
						case 102:
							//write timestamp to frame header
							MPI_File_write(fh, &frameHeader_timestamp, 1, MPI_FLOAT, &status);
							break;
							
						default:
							cout << "Error mmpld-writer: file version " << MMPLD_FILE_VERSION << " not supported." << endl;
							return;
							break;
					}
					
					uint32_t frameHeader_numPLists = htole32(_numSphereTypes);
					MPI_File_write(fh, &frameHeader_numPLists, 1, MPI_UNSIGNED, &status);
					
				}
				uint8_t pListHeader_vortexType;
				uint8_t pListHeader_colorType;
				float pListHeader_globalRadius;
				uint8_t pListHeader_red;
				uint8_t pListHeader_green;
				uint8_t pListHeader_blue;
				uint8_t pListHeader_alpha;
				uint64_t pListHeader_particleCount;
				
				//set vortex data type to FLOAT_XYZ
				pListHeader_vortexType = 1;
				
				//set color data type to NONE (only global color used)
				pListHeader_colorType = 0;

				//select different colors depending on component|site id
				pListHeader_globalRadius = _vfSphereRadius[nSphereTypeIndex];
				pListHeader_red   = _vaSphereColors[nSphereTypeIndex][0];
				pListHeader_green = _vaSphereColors[nSphereTypeIndex][1];
				pListHeader_blue  = _vaSphereColors[nSphereTypeIndex][2];
				pListHeader_alpha = _vaSphereColors[nSphereTypeIndex][3];
				
				//store componentParticleCount
				pListHeader_particleCount = htole64(globalNumCompSpheres[nSphereTypeIndex]);
				
				MPI_File_write(fh, &pListHeader_vortexType, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_colorType, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_globalRadius, 1, MPI_FLOAT, &status);
				MPI_File_write(fh, &pListHeader_red, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_green, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_blue, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_alpha, 1, MPI_BYTE, &status);
				MPI_File_write(fh, &pListHeader_particleCount, 1, MPI_LONG_LONG_INT, &status);
			}
			
			float spherePos[3];
			for (Molecule* mol = particleContainer->begin(); mol != particleContainer->end(); mol = particleContainer->next())
			{
				if(true == GetSpherePos(spherePos, mol, nSphereTypeIndex) )
					MPI_File_write(fh, spherePos, 3, MPI_FLOAT, &status);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_File_close(&fh);
#endif
	}
}

void MmpldWriter::finishOutput(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* /*domain*/) {
	//fill seektable
	stringstream filenamestream;
	filenamestream << _outputPrefix;

	if(_appendTimestamp) {
		filenamestream << "-" << _timestampString;
	}
	filenamestream << ".mmpld";

	char filename[filenamestream.str().size()+1];
	strcpy(filename,filenamestream.str().c_str());

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	if (rank == 0){
		
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_Offset endPosition;
		MPI_File_get_position(fh, &endPosition);
		_seekTable[_numSeekEntries-1] = (uint64_t)endPosition;
		MPI_File_seek(fh, 0x3C, MPI_SEEK_SET);
		uint64_t seekPosition;
		MPI_Status status;
		for (uint32_t i = 0; i < _numSeekEntries; ++i){
			seekPosition = htole64(_seekTable[i]);
			MPI_File_write(fh, &seekPosition, 1, MPI_LONG_LONG_INT, &status);
		}
		delete[] _seekTable;
		MPI_File_close(&fh);
	}else{
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_close(&fh);
	}
		
#endif
}

void MmpldWriter::InitSphereData()
{
	if(_bInitSphereData == ISD_USE_DEFAULT)
	{
		for(uint8_t i=0; i<6; i++)
			_vfSphereRadius.push_back(0.5);

		//                            R    G    B  alpha
		_vaSphereColors.push_back( {255,   0,   0, 255} );  // red
		_vaSphereColors.push_back( {  0, 205, 255, 255} );  // lightblue
		_vaSphereColors.push_back( {255,   0, 255, 255} );  // blue
		_vaSphereColors.push_back( {  0, 155,   0, 255} );  // green
		_vaSphereColors.push_back( {105,   0, 205, 255} );  // purple
		_vaSphereColors.push_back( {255, 125,   0, 255} );  // orange

		return;
	}

	std::ifstream filein(_strSphereDataFilename.c_str(), ios::in);
	std::string strLine, strToken;
	std::string strTokens[6];
	std::array<uint8_t, 4> arrColors;

	while (getline (filein, strLine))
	{
		stringstream sstr;
		sstr << strLine;
//		cout << sstr.str() << endl;

		uint8_t ti=0;
		while (sstr >> strToken)
			if(ti<6) strTokens[ti++] = strToken;

		if(ti==6 && strTokens[0][0] != '#')
		{
			_vfSphereRadius.push_back( (float)(atof( strTokens[1].c_str() ) ) );
			for(uint8_t ci=0; ci<4; ci++)
				arrColors[ci] = (uint8_t)(atoi(strTokens[ci+2].c_str() ) );
			_vaSphereColors.push_back(arrColors);
		}
	}

#ifndef NDEBUG
	std::vector<float>::iterator it; int i=0; cout << "radii" << endl;
	for(it=_vfSphereRadius.begin(); it!=_vfSphereRadius.end(); ++it)
		cout << i++ << ": " << (*it) << endl;

	std::vector< std::array<uint8_t, 4> >::iterator cit; i=0; cout << "colors" << endl;
	for(cit=_vaSphereColors.begin(); cit!=_vaSphereColors.end(); ++cit)
	{
		cout << i++ << ":";
		for(int j=0; j<4; j++)
			 cout << setw(4) << (int)(*cit).data()[j];
		cout << endl;
	}
#endif
}


// derived classes
void MmpldWriterSimpleSphere::CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol)
{
	uint8_t cid = mol->componentid();
	numSpheresPerType[cid]++;
}

bool MmpldWriterSimpleSphere::GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex)
{
	uint8_t cid = mol->componentid();
	for (unsigned short d = 0; d < 3; ++d) spherePos[d] = (float)mol->r(d);
	return (cid == nSphereTypeIndex);
}


void MmpldWriterMultiSphere::CalcNumSpheresPerType(uint64_t* numSpheresPerType, Molecule* mol)
{
	uint8_t cid = mol->componentid();
	uint8_t offset = _nCompSitesOffset[cid];
	for (uint8_t si = 0; si < _numSitesPerComp[cid]; ++si)
		numSpheresPerType[offset+si]++;
}

bool MmpldWriterMultiSphere::GetSpherePos(float (&spherePos)[3], Molecule* mol, uint8_t& nSphereTypeIndex)
{
	bool ret = false;
	uint8_t cid = mol->componentid();
	uint8_t numSites =  _numSitesPerComp[cid];
	uint8_t offset  = _nCompSitesOffset[cid];
	for (uint8_t si = 0; si < numSites; ++si)
	{
		if(offset+si == nSphereTypeIndex)
		{
			const std::array<double,3> arrSite = mol->ljcenter_d_abs(si);
			const double* posSite = arrSite.data();
			for (unsigned short d = 0; d < 3; ++d) spherePos[d] = (float)posSite[d];
			ret = true;
		}
	}
	return ret;
}



















