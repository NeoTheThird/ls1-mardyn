// TraWriter.cpp

#include "io/TraWriter.h"
#include "io/TraTrack.h"
// #include "io/CxWriter.h"  // to access C(x) profile
#include "io/Region.h"
#include "Common.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include <list>

#ifdef ENABLE_MPI
#include <mpi.h>
#include <sys/resource.h>
#endif

using namespace std;

TraWriter::TraWriter( unsigned long writeFrequency, unsigned long nPhaseBoundaryTrackFreq, unsigned long nNumTrackTimesteps,
		              bool incremental)
{
	_writeFrequency = writeFrequency;
	_nPhaseBoundaryTrackFreq = nPhaseBoundaryTrackFreq;
	_nNumTrackTimesteps = nNumTrackTimesteps;

	_incremental = incremental;
	_filenameisdate = false;

	// TODO: IDs für Phasengrenze anders, an anderer Stelle vergeben
	_nPhaseBoundID = 0;
}

TraWriter::~TraWriter(){}

void TraWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain)
{
	/*
	string filename = _filename + ".tra";
	ofstream fileout(filename.c_str(), ios::out);
	fileout.close();
	_wroteTra = false;
	*/
}

void TraWriter::readXML(XMLfileUnits& xmlconfig)
{
	/*
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
*/
}

void TraWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu)
{
	// write output with respect to write frequency
	if ( !(simstep % _writeFrequency == 0) )
		return;

	// TraTrackliste aktualisieren
	std::list<TraTrack*>::iterator itTrack;
	std::vector<Region*> regionList = domain->GetRegionList();
	std::vector<Region*>::iterator itReg;
	unsigned short nRegionType;
	int nNumUnfinishedTracks;
	unsigned int nCompID;

	UpdateTrackList(particleContainer, domainDecomp, domain, simstep);
	nNumUnfinishedTracks = 0;

	// global_log->info() << "Counting unfinished tracks..." << endl;

	for(itTrack = _TraTrackList.begin(); itTrack != _TraTrackList.end(); itTrack++)
	{
		if( !(*itTrack)->TrackFinished() )
			nNumUnfinishedTracks++;

		// create simstep directory
		(*itTrack)->CreateTimestepDir();

		// create phase boundary directory
		(*itTrack)->CreatePhaseBoundDir();
	}
	// global_log->info() << "Working on: " << nNumUnfinishedTracks << " Tracks." << endl;


	// TraTrackliste durchgehen
	for(itTrack = _TraTrackList.begin(); itTrack != _TraTrackList.end(); ++itTrack)
	{
		// für abgeschlossene Tracks muss nichts gemacht werden
		if( (*itTrack)->TrackFinished() )
			continue;

		struct rlimit limit;
		/* Set the stack limit in bytes. */
		//
		limit.rlim_cur = 1099511627776;
		limit.rlim_max = 1099511627776;
		if (setrlimit(RLIMIT_STACK, &limit) != 0)
		{
			cout << "setrlimit() failed.";
			exit(1);
		}

		std::list<TrackMoleculeData*>::iterator itMol;
		std::list<TrackMoleculeData*> trackMolecules = (*itTrack)->GetTrackMolecules();

		for(itMol = trackMolecules.begin(); itMol != trackMolecules.end(); itMol++)
		{
			std::stringstream outputstream;
			std::stringstream filepathstream;

			// write header
			if ( (*itTrack)->WroteHeader( (*itMol)->GetID() ) == false )
			{
				outputstream << "ZS         x          y          z          R \n";
				(*itTrack)->SetWroteHeader( (*itMol)->GetID(), true);
			}
			else
			{
				// outputstream << "#" << endl;
			}

			outputstream << simstep;

			for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next() )
			{
				bool halo = false;
				for (unsigned short d = 0; d < 3; d++) {
					if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
						halo = true;
						break;
					}
				}
				if (!halo && pos->id() == (*itMol)->GetID() )
				{
					for (unsigned short d = 0; d < 3; d++)
					{
						outputstream << std::setw(11) << std::setprecision(6) << pos->r(d);
					}

					// get region type
					nRegionType = -1;  // invalid type

					for(itReg = regionList.begin(); itReg != regionList.end(); itReg++)
					{
						if( (*itReg)->MoleculeIsInside(pos) )
							nRegionType = (*itReg)->GetType();
					}

					outputstream << std::setw(11) << nRegionType << endl;

					// get component id
					nCompID = pos->component()->ID() + 1;  // ls1 component ids starts with 0, input files with 1 --> transform with: +1
				}

			}  // for each molecule in particleContainer

			filepathstream << "./" << (*itTrack)->GetTimestepDirName() << "/" << (*itTrack)->GetPhaseBoundDirName() << "/"
						   << "mol-id_" << (*itMol)->GetID() << "_Comp-id_" << nCompID << ".tra";

			char filepath[filepathstream.str().size()+1];
			strcpy(filepath,filepathstream.str().c_str());

			// convert string
			long outputsize = outputstream.str().size();
			//cout << "rank: " << rank << "; step: " << simstep << "; outputsize: " << outputsize << endl;
			char output[outputsize+1];
			strcpy(output,outputstream.str().c_str());

			// write to file
			ofstream fileout(filepath, ios::out|ios::app);
			fileout << output;
			fileout.close();

		}  // for each track molecule of track

		// global_log->info() << "TrackID: " << (*itTrack)->GetTrackID() << "Timestep " << simstep << " recorded." << endl;
		(*itTrack)->InformTrack(TTM_TIMESTEP_RECORDED);

	}  // work on whole tracklist

}

void TraWriter::finishOutput(ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp, Domain* domain) {
}

void TraWriter::UpdateTrackList(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, const unsigned long &simstep)
{
	// Nur entsprechend update-frequenz updaten
	if ( !(simstep % _nPhaseBoundaryTrackFreq == 0) )
	{
		return;
	}

	stringstream strStream;
	std::vector<Region*> regionList = domain->GetRegionList();
	std::vector<Region*>::iterator itReg;
	int nRet;

	TraTrack* track;
	double dLowerCorner[3];
	double dUpperCorner[3];
	double dBoundaryMidpoint[3];

	global_log->info() << "durchlaufe Regionsliste mit " << regionList.size() << " Regionen..." << endl;

	for(itReg = regionList.begin(); itReg != regionList.end(); itReg++)
	{
		// skip region with right boundary equals domain right boundary
		if( (*itReg)->GetRightBoundary() == domain->getGlobalLength(0) )
			continue;

		dBoundaryMidpoint[0] = (*itReg)->GetRightBoundary();
		dBoundaryMidpoint[1] = (*itReg)->GetMidpoint(1);
		dBoundaryMidpoint[2] = (*itReg)->GetMidpoint(2);

		dLowerCorner[0] = dBoundaryMidpoint[0] - 50.0;
		dLowerCorner[1] = dBoundaryMidpoint[1] - 50.0;
		dLowerCorner[2] = dBoundaryMidpoint[2] - 50.0;

		dUpperCorner[0] = dBoundaryMidpoint[0] + 50.0;
		dUpperCorner[1] = dBoundaryMidpoint[1] + 50.0;
		dUpperCorner[2] = dBoundaryMidpoint[2] + 50.0;

		global_log->info() << "LowerCorner x: " << dLowerCorner[0] <<  ", y: " << dLowerCorner[1] <<  ", z: " << dLowerCorner[2] << endl;
		global_log->info() << "UpperCorner x: " << dUpperCorner[0] <<  ", y: " << dUpperCorner[1] <<  ", z: " << dUpperCorner[2] << endl;

		// get molecules from region
		std::list<Molecule*> particlePtrs;
		std::list<Molecule*>::iterator itMol;
		particleContainer->getRegion(dLowerCorner, dUpperCorner, particlePtrs);

		global_log->info() << "particlePtrs - size: " << particlePtrs.size() << endl;

		// create list of Molecule IDs
		std::list<TrackMoleculeData*> trackMolecules;

		for(itMol = particlePtrs.begin(); itMol != particlePtrs.end(); itMol++)
		{
			trackMolecules.push_back( new TrackMoleculeData( (*itMol)->id() ) );
		}

		_nPhaseBoundID++;
		track = new TraTrack(_nNumTrackTimesteps, trackMolecules, simstep, _nPhaseBoundID);
		_TraTrackList.push_back(track);
	}
}


/*
 *
void TraWriter::UpdateTrackList(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, const unsigned long &simstep)
{
	// Nur entsprechend update-frequenz updaten
	if ( !(simstep % _nPhaseBoundaryTrackFreq == 0) )
	{
		return;
	}

	stringstream strStream;
	std::list<unsigned int> nTrackMoleculeIDs;
	std::list<double*> phaseBoundMidpointList;
	std::list<double*>::iterator itMid;
	double* dBoundaryMidpoint;
	int nRet;

	// nRet = domain->FindPhaseBoundarys(phaseBoundMidpointList);
	nRet = this->FindPhaseBoundaryMidpointsX(particleContainer, domainDecomp, domain, phaseBoundMidpointList);

	// check if midpoints found
	if(nRet < 1)
	{
		global_log->info() << "No midpoints found, returning..." << endl;
		return;
	}

	TraTrack* track;
	std::list<Molecule*> particlePtrs;
	std::list<Molecule*>::iterator it;
	double dLowerCorner[3];
	double dUpperCorner[3];

	for(itMid = phaseBoundMidpointList.begin(); itMid != phaseBoundMidpointList.end(); itMid++)
	{
		dBoundaryMidpoint = *itMid;

		dLowerCorner[0] = dBoundaryMidpoint[0] - 5.0;
		dLowerCorner[1] = dBoundaryMidpoint[1] - 5.0;
		dLowerCorner[2] = dBoundaryMidpoint[2] - 5.0;

		dUpperCorner[0] = dBoundaryMidpoint[0] + 5.0;
		dUpperCorner[1] = dBoundaryMidpoint[1] + 5.0;
		dUpperCorner[2] = dBoundaryMidpoint[2] + 5.0;

		global_log->info() << "LowerCorner x: " << dLowerCorner[0] <<  ", y: " << dLowerCorner[1] <<  ", z: " << dLowerCorner[2] << endl;
		global_log->info() << "UpperCorner x: " << dUpperCorner[0] <<  ", y: " << dUpperCorner[1] <<  ", z: " << dUpperCorner[2] << endl;

		particlePtrs.clear();
		particleContainer->getRegion(dLowerCorner, dUpperCorner, particlePtrs);

		global_log->info() << "particlePtrs - size: " << particlePtrs.size() << endl;

		// create list of Molecule IDs
		nTrackMoleculeIDs.empty();

		for(it = particlePtrs.begin(); it != particlePtrs.end(); ++it)
		{
			nTrackMoleculeIDs.push_back( (*it)->id() );
		}

		_nPhaseBoundID++;
		track = new TraTrack(_nNumTrackTimesteps, nTrackMoleculeIDs, simstep, _nPhaseBoundID);
		_TraTrackList.push_back(track);
	}

}

*/


int TraWriter::FindPhaseBoundaryMidpointsX( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
		                                    std::list<double*> &phaseBoundMidpointList)
{
	// CxWriter::_Cx.begin();

	double dBoxMidpoint[3];
	double dBoxLength[3];
	double dSliceWidthX;
	int nNumSlices = 100;
	int prevID;

	std::list<Molecule*> particlePtrs;
	std::list<Molecule*>::iterator it;
	double dLowerCorner[3];
	double dUpperCorner[3];
	double dPhaseBoundMidpoint[3];

	dBoxLength[0] = domain->getGlobalLength(0);
	dBoxLength[1] = domain->getGlobalLength(1);
	dBoxLength[2] = domain->getGlobalLength(2);

	dSliceWidthX = dBoxLength[0] / nNumSlices;

	dBoxMidpoint[0] = dBoxLength[0] / 2;
	dBoxMidpoint[1] = dBoxLength[1] / 2;
	dBoxMidpoint[2] = dBoxLength[2] / 2;

	// phase midpoint stays on x-axis through mid of y,z-area
	dPhaseBoundMidpoint[1] = dBoxMidpoint[1];
	dPhaseBoundMidpoint[2] = dBoxMidpoint[2];

	global_log->info() << "searching for phase boundary..." << endl;

	dLowerCorner[1] = 0;
	dLowerCorner[2] = 0;
	dUpperCorner[1] = dBoxLength[1];
	dUpperCorner[2] = dBoxLength[2];

	for(int i=0; i < nNumSlices / 2; i++)  // sucht in beide Richtungen
	{
		global_log->info() << "slice nr: " << i << endl;

		dLowerCorner[0] = dBoxMidpoint[0] - (i+1) * dSliceWidthX;
		dUpperCorner[0] = dBoxMidpoint[0] + (i+1) * dSliceWidthX;

		global_log->info() << "LowerCorner x: " << dLowerCorner[0] <<  ", y: " << dLowerCorner[1] <<  ", z: " << dLowerCorner[2] << endl;
		global_log->info() << "UpperCorner x: " << dUpperCorner[0] <<  ", y: " << dUpperCorner[1] <<  ", z: " << dUpperCorner[2] << endl;

		particlePtrs.clear();
		particleContainer->getRegion(dLowerCorner, dUpperCorner, particlePtrs);
		prevID = -1;

		global_log->info() << "particlePtrs - size: " << particlePtrs.size() << endl;

		_nNumCompInBox = 1;

		for(it = particlePtrs.begin(); it != particlePtrs.end(); ++it)
		{
			int cid = (*it)->componentid();
			if(prevID > 0 && prevID != cid)
				_nNumCompInBox++;
			prevID = cid;
		}

		global_log->info() << "Number of Components found: " << _nNumCompInBox << endl;


		// Wert akkumulieren mit anderen Prozessen
		domainDecomp->collCommInit(1);
		domainDecomp->collCommAppendInt(_nNumCompInBox);
		domainDecomp->collCommAllreduceSum();
		_nNumCompInBox = domainDecomp->collCommGetInt();
		domainDecomp->collCommFinalize();

		global_log->info() << "_nNumCompInBox: " << _nNumCompInBox << endl;
		global_log->info() << "Number of Processes: " << domainDecomp->getNumProcs() << endl;

		_nNumCompInBox = _nNumCompInBox / domainDecomp->getNumProcs();
		if(_nNumCompInBox > 1)  // mehr als eine Phase anwesend
		{
			dPhaseBoundMidpoint[0] = dLowerCorner[0];
			global_log->info() << "dPhaseBoundMidpoint: x= " << dPhaseBoundMidpoint[0] << ", y= " << dPhaseBoundMidpoint[1] << ", z= " << dPhaseBoundMidpoint[2] << endl;
			phaseBoundMidpointList.push_back(dPhaseBoundMidpoint);

			dPhaseBoundMidpoint[0] = dUpperCorner[0];
			global_log->info() << "dPhaseBoundMidpoint: x= " << dPhaseBoundMidpoint[0] << ", y= " << dPhaseBoundMidpoint[1] << ", z= " << dPhaseBoundMidpoint[2] << endl;
			phaseBoundMidpointList.push_back(dPhaseBoundMidpoint);
			break;
		}
	}

	std::list<double*>::iterator itM;
	double* point;

	for(itM = phaseBoundMidpointList.begin(); itM != phaseBoundMidpointList.end(); itM++)
	{
		point = (*itM);
		global_log->info() << "new phase boundary midpoint: x= " << point[0] << ", y= " << point[1] << ", z= " << point[2] << endl;
	}

	return phaseBoundMidpointList.size();
}























