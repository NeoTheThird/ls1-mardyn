/*
 * DensityControl.cpp
 *
 *  Created on: 29.05.2015
 *      Author: mheinen
 */

#include "NEMD/DensityControl.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
//#include <cmath>
//#include <list>
//#include <map>
#include "Domain.h"
//#include "NEMD/ParticleInsertion.h"
//#include "utils/Random.h"
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
//#include <iterator>  // std::advance

#include <cstdlib>

using namespace std;

// init static ID --> instance counting
unsigned short dec::ControlRegion::_nStaticID = 0;

// class dec::ControlRegion

dec::ControlRegion::ControlRegion(ControlInstance* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity)
: CuboidRegionObs(parent, dLowerCorner, dUpperCorner)
{
	// ID
	_nID = ++_nStaticID;

    // target density
    _dTargetDensity = dTargetDensity;

    // target component ID (inert gas)
    _nTargetComponentID = nTargetComponentID;

    // init process relevance
    _bProcessIsRelevant = true;

    // init rank array
    _ranks = NULL;

    // reset local values
    _nDeletedNumMoleculesLocal = 0;

    for(unsigned int d=0; d<3; ++d)
    {
        _dDeletedVelocityLocal[d] = 0.;
        _dDeletedVelocitySquaredLocal[d] = 0.;
        _dDeletedEkinLocal[d] = 0.;
    }

    this->WriteHeaderDeletedMolecules();
}


dec::ControlRegion::~ControlRegion()
{
}

void dec::ControlRegion::CheckBounds()
{
	// check if control region is outside of simulation volume
	double dBoxLengthY = _parent->GetDomain()->getGlobalLength(1);

	if ( (_dLowerCorner[1] > dBoxLengthY && _dUpperCorner[1] > dBoxLengthY) || (_dLowerCorner[1] < 0. && _dUpperCorner[1] < 0.) )
	{
		cout << "dec::ControlRegion::ControlRegion: Control region dimensions are outside of simulation volume! Programm exit..." << endl;
		exit(-1);
	}
}

void dec::ControlRegion::Init()
{
	// update region volume
	this->UpdateVolume();
}

void dec::ControlRegion::InitMPI()
{
	// domain decomposition
    DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

//    int numprocs = domainDecomp->getNumProcs();
    int nRank = domainDecomp->getRank();

    // get number of relevant processes
    int nRelevant;
    int nNumRelevantGlobal;

    double bbMin[3];
    double bbMax[3];

    domainDecomp->getBoundingBoxMinMax(_parent->GetDomain(), bbMin, bbMax);

    if( (bbMin[1] < _dLowerCorner[1] && bbMax[1] < _dLowerCorner[1]) ||  (bbMin[1] > _dUpperCorner[1] && bbMax[1] > _dUpperCorner[1]) )
        nRelevant = 0;
    else
        nRelevant = 1;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendInt(nRelevant);
    domainDecomp->collCommAllreduceSum();
    nNumRelevantGlobal = domainDecomp->collCommGetInt();
    domainDecomp->collCommFinalize();

    // cout << "nNumRelevantGlobal = " << nNumRelevantGlobal << endl;

    // allocate rank array
    if(_ranks != NULL)
    	delete _ranks;

    _ranks = new int[nNumRelevantGlobal];

    int nUnregistered = 1;

    for(int r=nNumRelevantGlobal-1; r>=0; --r)
    {
        int nRankLocal = (nRank+1) * nRelevant * nUnregistered;
        int nRankGlobal;

       // cout << "nRankLocal = " << nRankLocal << endl;

        MPI_Allreduce( &nRankLocal, &nRankGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        _ranks[r] = nRankGlobal-1;

        if(nRankGlobal-1 == nRank)
            nUnregistered = 0;
    }

    // check / set process relevance
    bool bProcessIsRelevant = false;

    for(int r=0; r<nNumRelevantGlobal; ++r)
    {
      //  cout << "_ranks[" << r << "] = " << _ranks[r] << endl;

        if(nRank == _ranks[r])
            bProcessIsRelevant = true;
    }

    if(true == bProcessIsRelevant)
        _bProcessIsRelevant = true;
    else
        _bProcessIsRelevant = false;

    MPI_Group group, newgroup;

    MPI_Comm_group(MPI_COMM_WORLD, &group);

    MPI_Group_incl(group, nNumRelevantGlobal, _ranks, &newgroup);
    MPI_Comm_create(MPI_COMM_WORLD, newgroup, &_newcomm);

    int groupRank;

    MPI_Group_rank(group, &groupRank);

  //  cout << "[" << nRank << "]: " << "groupRank = " << groupRank << endl;

    MPI_Group_rank(newgroup, &groupRank);

  //  cout << "[" << nRank << "]: " << "groupRank = " << groupRank << endl;

    int nRankLocal = nRank;
    int nRankGlobal;

    if (_newcomm != MPI_COMM_NULL)
    {
        MPI_Allreduce( &nRankLocal, &nRankGlobal, 1, MPI_INT, MPI_MAX, _newcomm);
  //      cout << "[" << nRank << "]: " << "NEW COMM! nRankGlobal = " << nRankGlobal << endl;
    }
    else
    {
 //       cout << "new comm is NULL!" << endl;
    }
}

void dec::ControlRegion::CalcGlobalValues()
{
	// In vacuum evaporation simulation global density of control region has not to be calculated
	if( 0. == _dTargetDensity)
		return;

#ifdef ENABLE_MPI
//    if (_newcomm == MPI_COMM_NULL)
//        return;

    MPI_Allreduce( &_nNumMoleculesLocal, &_nNumMoleculesGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); // _newcomm);

    _dDensityGlobal = _nNumMoleculesGlobal * _dInvertVolume;
#else
    _dDensityGlobal = _nNumMoleculesLocal * _dInvertVolume;
#endif

}

void dec::ControlRegion::UpdateGlobalDensity(bool bDeleteMolecule)
{
	// domain decomposition
	DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

    double dDensityLocal = _dDensityGlobal;
    int nDeletionsLocal = 0 + bDeleteMolecule;
    int nDeletionsGlobal = 0;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendInt(nDeletionsLocal);
    domainDecomp->collCommAllreduceSum();
    nDeletionsGlobal = domainDecomp->collCommGetInt();
    domainDecomp->collCommFinalize();

    // update density
    dDensityLocal -= nDeletionsGlobal * _dInvertVolume;

    // collective Communication
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendDouble(dDensityLocal);
    domainDecomp->collCommAllreduceSum();
    _dDensityGlobal = domainDecomp->collCommGetDouble();
    domainDecomp->collCommFinalize();
}

void dec::ControlRegion::MeasureDensity(Molecule* mol)
{
	// In vacuum evaporation simulation global density of control region has not to be calculated
	if( 0. == _dTargetDensity)
		return;

    // check if molecule inside control region
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    // count num molecules
    _nNumMoleculesLocal++;
}


void dec::ControlRegion::ControlDensity(Molecule* mol, Simulation* simulation, bool& bDeleteMolecule)
{
    /*
    // check component ID: if inert gas --> do nothing (return)
    if(mol->componentid()+1 == 3)  // program intern component ID starts with 0
        return;
*/

//    int nRank = domainDecomp->getRank();
//    double dPosY = mol->r(1);

    // check if molecule is inside
    for(unsigned short d = 0; d<3; ++d)
    {
        double dPos = mol->r(d);

        if(dPos <= _dLowerCorner[d] || dPos >= _dUpperCorner[d] )
            return;
    }

    if( 0. == _dTargetDensity)
    {
        bDeleteMolecule = true;

        // sample deleted molecules data
        _nDeletedNumMoleculesLocal++;
        _dDeletedEkinLocal[0] += mol->U_kin();
        _dDeletedEkinLocal[1] += mol->U_trans();
        _dDeletedEkinLocal[2] += mol->U_rot();

        for(unsigned short d = 0; d<3; ++d)
        {
        	double v = mol->v(d);
            _dDeletedVelocityLocal[d] += v;
            _dDeletedVelocitySquaredLocal[d] += v*v;
        }
    }
    else
    {
        /* initialize random seed: */
        srand (time(NULL) + mol->id() );

        /* generate secret number between 0 and 99999: */
        int nRand = rand() % 100000;

        double dRand = (double) (nRand / 100000.);

    //    cout << "dRand = " << dRand << endl;
    //
    //    cout << "_dDensityGlobal = " << _dDensityGlobal << endl;
    //    cout << "_dTargetDensity = " << _dTargetDensity << endl;

        double dPercentToTakeOut = (_dDensityGlobal - _dTargetDensity) / _dDensityGlobal;

    //    cout << "dPercentToTakeOut = " << dPercentToTakeOut << endl;

        if(dPercentToTakeOut > 0. && dRand < abs(dPercentToTakeOut) )
            bDeleteMolecule = true;
        else
            bDeleteMolecule = false;
    }
}

void dec::ControlRegion::ResetLocalValues()
{
    _nNumMoleculesLocal = 0;
}

void dec::ControlRegion::WriteHeaderDeletedMolecules()
{
	// domain decomposition
	DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

#ifdef ENABLE_MPI
int rank = domainDecomp->getRank();
// int numprocs = domainDecomp->getNumProcs();
if (rank != 0)
    return;
#endif

    // write header
    stringstream outputstream;
    stringstream sstrFilename;    
    sstrFilename << "DensityControl_del-mol-data_region" << this->GetID() << ".dat";

    outputstream << "             simstep";
    outputstream << "         numMols";
    outputstream << "           U_kin";
    outputstream << "         U_trans";
    outputstream << "           U_rot";
    outputstream << "              vx";
    outputstream << "              vy";
    outputstream << "              vz";
    outputstream << "             vx2";
    outputstream << "             vy2";
    outputstream << "             vz2";
    outputstream << endl;

    ofstream fileout(sstrFilename.str().c_str(), ios::out);
    fileout << outputstream.str();
    fileout.close();
}

void dec::ControlRegion::WriteDataDeletedMolecules(unsigned long simstep)
{
	// domain decomposition
	DomainDecompBase* domainDecomp = _parent->GetDomainDecomposition();

    // calc global values 
#ifdef ENABLE_MPI

    MPI_Reduce( &_nDeletedNumMoleculesLocal,    &_nDeletedNumMoleculesGlobal,    1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  _dDeletedEkinLocal,             _dDeletedEkinGlobal,            3, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  _dDeletedVelocityLocal,         _dDeletedVelocityGlobal,        3, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(  _dDeletedVelocitySquaredLocal,  _dDeletedVelocitySquaredGlobal, 3, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);

#else
    _nDeletedNumMoleculesGlobal = _nDeletedNumMoleculesLocal;

    for(unsigned int d=0; d<3; ++d)
    {
    	_dDeletedEkinGlobal[d] = _dDeletedEkinLocal[d];
        _dDeletedVelocityGlobal[d] = _dDeletedVelocityLocal[d];
        _dDeletedVelocitySquaredGlobal[d] = _dDeletedVelocitySquaredLocal[d];
    }
#endif

    // reset local values
    _nDeletedNumMoleculesLocal = 0;

    for(unsigned int d=0; d<3; ++d)
    {
        _dDeletedEkinLocal[d] = 0.;
        _dDeletedVelocityLocal[d] = 0.;
        _dDeletedVelocitySquaredLocal[d] = 0.;
    }

    // write out data
    #ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank != 0)
        return;
    #endif

    stringstream outputstream;
    stringstream sstrFilename;
    sstrFilename << "DensityControl_del-mol-data_region" << this->GetID() << ".dat";

    outputstream << std::setw(20) << simstep;
    outputstream << std::setw(16) << _nDeletedNumMoleculesGlobal;

    for(unsigned int d=0; d<3; ++d)
        outputstream << std::setw(16) << fixed << std::setprecision(3) << _dDeletedEkinGlobal[d];

    for(unsigned int d=0; d<3; ++d)
        outputstream << std::setw(16) << fixed << std::setprecision(3) << _dDeletedVelocityGlobal[d];

    for(unsigned int d=0; d<3; ++d)
        outputstream << std::setw(16) << fixed << std::setprecision(3) << _dDeletedVelocitySquaredGlobal[d];

    outputstream << endl;

    ofstream fileout(sstrFilename.str().c_str(), ios::app);
    fileout << outputstream.str();
    fileout.close();
}
        
// class DensityControl

DensityControl::DensityControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop)
: ControlInstance(domain, domainDecomp)
{
    // control frequency
    _nControlFreq = nControlFreq;

    // start/stop timestep
    _nStart = nStart;
    _nStop  = nStop;

    // deleted molecules data
    _nWriteFreqDeleted = 1000;
}

DensityControl::~DensityControl()
{
}

void DensityControl::AddRegion(dec::ControlRegion* region)
{
    _vecControlRegions.push_back(region);

    // check / set process relevance
    std::vector<dec::ControlRegion*>::iterator it;
    bool bProcessIsRelevant = false;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        if ( true == (*it)->ProcessIsRelevant() );
            bProcessIsRelevant = true;
    }

    if(true == bProcessIsRelevant)
        _bProcessIsRelevant = true;
    else
        _bProcessIsRelevant = false;
}

void DensityControl::MeasureDensity(Molecule* mol, unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // measure drift in each control region
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->MeasureDensity(mol);
    }
}

void DensityControl::CalcGlobalValues(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CalcGlobalValues();
    }
}

void DensityControl::UpdateGlobalDensities(unsigned long simstep, bool bDeleteMolecule)
{
    if(simstep % _nControlFreq != 0)
        return;

    // calc global values for control region
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->UpdateGlobalDensity(bDeleteMolecule);
    }
}

void DensityControl::ControlDensity(Molecule* mol, Simulation* simulation, unsigned long simstep, bool& bDeleteMolecule)
{
    if(simstep % _nControlFreq != 0)
        return;

    // control drift of all regions
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ControlDensity(mol, simulation, bDeleteMolecule);
    }
}

void DensityControl::CheckRegionBounds()
{
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->CheckBounds();
    }
}

void DensityControl::Init(unsigned long simstep)
{
    if(simstep % _nControlFreq != 0)
        return;

    // reset local values
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->ResetLocalValues();

        // Init
        (*it)->Init();
    }
}

void DensityControl::InitMPI()
{
    // reset local values
    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->InitMPI();
    }
}

void DensityControl::WriteDataDeletedMolecules(unsigned long simstep)
{
    if(simstep % _nWriteFreqDeleted != 0)
        return;

    std::vector<dec::ControlRegion*>::iterator it;

    for(it=_vecControlRegions.begin(); it!=_vecControlRegions.end(); ++it)
    {
        (*it)->WriteDataDeletedMolecules(simstep);
    }
}

