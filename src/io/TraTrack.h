/*
 * TraTrack.h
 *
 *  Created on: 21.08.2013
 *      Author: mheinen
 */

#ifndef TRATRACK_H_
#define TRATRACK_H_

#include <list>
#include <string>

enum TraTrackStatusTypes
{
	TTST_WROTE_HEADER = 1,
	TTST_TIMESTEP_DIRECTORY_CREATED = 2,
	TTST_PHASE_BOUNDARY_DIRECTORY_CREATED = 3,
};

enum TraTrackMessages
{
	TTM_TIMESTEP_RECORDED = 1,
};

class TraTrack
{
public:
	TraTrack();
	TraTrack( unsigned long nNumTrackTimesteps, std::list<unsigned int> lisTrackMoleculeIDs,
			  unsigned long nInitSimstep, unsigned int nPhaseBoundID);
	~TraTrack();

	std::list<unsigned int> GetTrackMoleculeIDs(){return _lisTrackMoleculeIDs;}
	std::string GetTimestepDirName(void)   {return _strTimestepDirName;}
	std::string GetPhaseBoundDirName(void) {return _strPhaseBoundDirName;}

	void InformTrack(int nMessage);
	bool TrackFinished(void) {return _bTrackFinished;}
	bool WroteHeader(void) {return _bWroteHeader;}
	bool TimestepDirCreated(void) {return _bTimestepDirCreated;}
	bool PhaseBoundDirCreated(void) {return _bPhaseBoundDirCreated;}
	void SetStatus(int nStatusType, bool bStatus);
	unsigned int GetTrackID(void) {return _nTrackID;}

	void CreateTimestepDir(void);
	void CreatePhaseBoundDir(void);

private:
	unsigned long _nNumTrackTimesteps;
	unsigned long _nNumTimestepsTracked;
	unsigned int  _nTrackID;
	std::string _strTimestepDirName;
	std::string _strPhaseBoundDirName;
	std::list<unsigned int> _lisTrackMoleculeIDs;
	bool _bTrackFinished;
	bool _bWroteHeader;
	bool _bTimestepDirCreated;
	bool _bPhaseBoundDirCreated;

};  // class TraTrack


#endif /* TRATRACK_H_ */
