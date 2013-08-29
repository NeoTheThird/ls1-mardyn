/*
 * TraTrack.cpp
 *
 *  Created on: 21.08.2013
 *      Author: mheinen
 */

#include "TraTrack.h"
#include <list>
#include <string>
#include <sstream>
#include <sys/types.h> // types for mkdir()
#include <sys/stat.h>  // to use mkdir()
#include <cstdlib>  // to use atoi()

TraTrack::TraTrack()
{
	// control params
	_nNumTrackTimesteps = 0;
	_nNumTimestepsTracked = 0;

	// directory names
	_strTimestepDirName = "";
	_strPhaseBoundDirName = "";

	// status variables
	_bTrackFinished = false;
	_bWroteHeader = false;
	_bTimestepDirCreated = false;
	_bPhaseBoundDirCreated = false;

	// Track ID
	_nTrackID = 0;
}

TraTrack::TraTrack( unsigned long nNumTrackTimesteps, std::list<unsigned int> lisTrackMoleculeIDs,
		            unsigned long nInitSimstep, unsigned int nPhaseBoundID)
{
	// control params
	_nNumTrackTimesteps = nNumTrackTimesteps;
	_nNumTimestepsTracked = 0;

	// status variables
	_bTrackFinished = false;
	_bWroteHeader = false;
	_bTimestepDirCreated = false;
	_bPhaseBoundDirCreated = false;

	// directory names
	std::ostringstream oss1, oss2;

	oss1 << "TS_" << nInitSimstep;
	oss2 << "PB_" << nPhaseBoundID;
	_strTimestepDirName   = oss1.str();
	_strPhaseBoundDirName = oss2.str();

	std::list<unsigned int>::iterator it;
	for(it = lisTrackMoleculeIDs.begin(); it != lisTrackMoleculeIDs.end(); ++it)
	{
		_lisTrackMoleculeIDs.push_back(*it);
	}

	// Track ID
	std::stringstream sst;
	sst << nInitSimstep << nPhaseBoundID;
	_nTrackID = atoi( sst.str().c_str() );
}

TraTrack::~TraTrack()
{
}

void TraTrack::SetStatus(int nStatusType, bool bStatus)
{
	switch(nStatusType)
	{
	case TTST_WROTE_HEADER:
		_bWroteHeader = bStatus;
		break;
	case TTST_TIMESTEP_DIRECTORY_CREATED:
		_bTimestepDirCreated = bStatus;
		break;
	case TTST_PHASE_BOUNDARY_DIRECTORY_CREATED:
		_bPhaseBoundDirCreated = bStatus;
		break;
	default:
		break;
	}
}

void TraTrack::InformTrack(int nMessage)
{
	switch(nMessage)
	{
	case TTM_TIMESTEP_RECORDED:
		_nNumTimestepsTracked++;
		_bTrackFinished = _nNumTimestepsTracked == _nNumTrackTimesteps+1;
		break;
	default:
		break;
	}
}

void TraTrack::CreateTimestepDir()
{
	// check if allready exist
	if( this->TimestepDirCreated() )
		return;

	// int status;
	// status = mkdir("/home/cnd/mod1", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir( this->GetTimestepDirName().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	this->SetStatus(TTST_TIMESTEP_DIRECTORY_CREATED, true);
}

void TraTrack::CreatePhaseBoundDir(void)
{
	// check if allready exist
	if( this->PhaseBoundDirCreated() )
		return;

	// check if timestep directory exist
	if( !this->TimestepDirCreated() )
	{
		this->CreateTimestepDir();
	}

	std::stringstream tmpstream;
	tmpstream << "./" << this->GetTimestepDirName() << "/" << this->GetPhaseBoundDirName();
	// int status;
	// status = mkdir("/home/cnd/mod1", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//	global_log->info() << "GetTimestepDirName: " << (*itTrack)->GetTimestepDirName() << endl;
//	global_log->info() << "GetPhaseBoundDirName: " << (*itTrack)->GetPhaseBoundDirName() << endl;
//	global_log->info() << "tmpstream: " << tmpstream.str().c_str() << endl;
	mkdir( tmpstream.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	this->SetStatus(TTST_TIMESTEP_DIRECTORY_CREATED, true);
}


