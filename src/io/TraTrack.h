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
	TTST_TIMESTEP_DIRECTORY_CREATED = 1,
	TTST_PHASE_BOUNDARY_DIRECTORY_CREATED = 2,
};

enum TraTrackMessages
{
	TTM_TIMESTEP_RECORDED = 1,
};

class TrackMoleculeData
{
public:
	TrackMoleculeData(unsigned long nID)
	{
		_nID = nID;
		_bWroteHeader = false;
	}
	~TrackMoleculeData() {}

	unsigned long GetID(void) {return _nID;}
	bool WroteHeader(void) {return _bWroteHeader;}
	void SetWroteHeader(bool bWroteHeader) {_bWroteHeader = bWroteHeader;}

private:
	unsigned long _nID;
	bool _bWroteHeader;

};  // class TrackMoleculeData

class TraTrack
{
public:
	TraTrack();
	TraTrack( unsigned long nNumTrackTimesteps, std::list<TrackMoleculeData*> trackMolecules,
			  unsigned long nInitSimstep, unsigned int nPhaseBoundID);
	~TraTrack();

	std::list<TrackMoleculeData*> GetTrackMolecules(){return _trackMolecules;}
	std::string GetTimestepDirName(void)   {return _strTimestepDirName;}
	std::string GetPhaseBoundDirName(void) {return _strPhaseBoundDirName;}

	void InformTrack(int nMessage);
	bool TrackFinished(void) {return _bTrackFinished;}
	bool WroteHeader(unsigned long nMolID);
	void SetWroteHeader(unsigned long nMolID, bool bWroteHeader);
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
	std::list<TrackMoleculeData*> _trackMolecules;
	bool _bTrackFinished;
	bool _bTimestepDirCreated;
	bool _bPhaseBoundDirCreated;

};  // class TraTrack


#endif /* TRATRACK_H_ */
