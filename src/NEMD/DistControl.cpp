/*
 * DistControl.cpp
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#include "DistControl.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


DistControl::DistControl(Domain* domain, unsigned int nUpdateFreq, unsigned int nNumShells, int nSimtype, double dVaporDensity, unsigned int nMethod)
{
    _nUpdateFreq = nUpdateFreq;
    _nNumShells = nNumShells;

    // number of molecules
    _nNumMoleculesLocal  = new unsigned long[_nNumShells];
    _nNumMoleculesGlobal = new unsigned long[_nNumShells];

    // density profile
    _dDensityProfile = new double[_nNumShells];
    _dDensityProfileSmoothed = new double[_nNumShells];
    _dDensityProfileSmoothedDerivation = new double[_nNumShells];

    // force profile
    _dForceSumLocal  = new double [_nNumShells];
    _dForceSumGlobal = new double [_nNumShells];
    _dForceProfile   = new double [_nNumShells];
    _dForceProfileSmoothed   = new double [_nNumShells];

    // reset local values
    this->ResetLocalValues();

    // calc shell width and shell volume
    double dNumShells = (double) (_nNumShells);

    _dShellWidth = domain->getGlobalLength(1) / dNumShells;
    _dInvertShellWidth = 1. / _dShellWidth;
    _dShellVolume = _dShellWidth * domain->getGlobalLength(0) * domain->getGlobalLength(2);

    // profile midpoint positions
    _dMidpointPositions = new double[_nNumShells];

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _dMidpointPositions[s] = (0.5 + s)*_dShellWidth;
    }

    // write data
    _nWriteFreq = _nUpdateFreq / 10;
    _nWriteFreqProfiles = 100000;

    _strFilename = "DistControl.dat";
    _strFilenameProfilesPrefix = "DistControlProfiles";

    // simtype
    _nSimType = nSimtype;

    // vapor density
    _dVaporDensity = dVaporDensity;

    // method of interface midpoint determination
    _nMethod = nMethod;  // 1: old one, 2: new one (density derivation)
}

DistControl::~DistControl()
{

}


void DistControl::SampleProfiles(Molecule* mol)
{
    unsigned int nPosIndex;
    unsigned int nIndexMax = _nNumShells - 1;

    // calc position index
    double dPos = mol->r(1);

    nPosIndex = (unsigned int) floor(dPos * _dInvertShellWidth);

    // ignore outer (halo) molecules
    if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
        return;

    _nNumMoleculesLocal[nPosIndex]++;
    _dForceSumLocal[nPosIndex] += mol->F(1);
}

void DistControl::EstimateInterfaceMidpointsByForce()
{
#ifdef ENABLE_MPI

    MPI_Allreduce( _nNumMoleculesLocal, _nNumMoleculesGlobal, _nNumShells, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _dForceSumLocal, _dForceSumGlobal, _nNumShells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _nNumMoleculesGlobal[s] = _nNumMoleculesLocal[s];
        _dForceSumGlobal[s] = _dForceSumLocal[s];
    }
#endif

    // calc profiles
    double dInvertShellVolume = 1. / _dShellVolume;
    double dInvertSampleTimesteps = 1. / ( (double)(_nUpdateFreq) );

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _dDensityProfile[s] = _nNumMoleculesGlobal[s] * dInvertSampleTimesteps * dInvertShellVolume;

        if(_nNumMoleculesGlobal[s] > 0) 
            _dForceProfile[s] = _dForceSumGlobal[s] / (double)(_nNumMoleculesGlobal[s]);
        else
            _dForceProfile[s] = 0.;
    }

    {
		// smooth profiles
		unsigned int nNumNeighbourVals = 20;
		double dNumSmoothedVals = (double)(2.*nNumNeighbourVals+1);

		unsigned int l = 0;
		unsigned int r = nNumNeighbourVals;

		for(unsigned int s = 0; s < _nNumShells; ++s)
		{
			double dDensitySum = 0.;
			double dForceSum   = 0.;

			for(unsigned int t = s-l; t <= s+r; ++t)
			{
				dDensitySum += _dDensityProfile[t];
				dForceSum   += _dForceProfile[t];
			}
			_dDensityProfileSmoothed[s] = dDensitySum / dNumSmoothedVals;
			_dForceProfileSmoothed[s]   = dForceSum   / dNumSmoothedVals;

			if(l < nNumNeighbourVals)
				l++;

			if(s >= (_nNumShells-1 - nNumNeighbourVals) )
				r--;
		}
    }

    // density derivation
    unsigned int nNumNeighbourVals = 3;
    
    unsigned int l = 0;
    unsigned int r = nNumNeighbourVals;

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        double dy = _dDensityProfileSmoothed[s+r] - _dDensityProfileSmoothed[s-l];
        double dx = (r+l)*_dShellWidth;

        _dDensityProfileSmoothedDerivation[s] = dy / dx;

        if(l < nNumNeighbourVals)
            l++;

        if(s >= (_nNumShells-1 - nNumNeighbourVals) )
            r--;
    }

    // find min/max in drho/dy profile
    double dMin = 0;
    double dMax = 0;
    unsigned int nIndexMin = 0;
    unsigned int nIndexMax = 0;

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        if(_dDensityProfileSmoothed[s] > _dVaporDensity)
        {
            double drhody = _dDensityProfileSmoothedDerivation[s];

            if(drhody > dMax)
            {
                dMax = drhody;
                nIndexMax = s;
            }

            if(drhody < dMin)
            {
                dMin = drhody;
                nIndexMin = s;
            }
        }
    }

/*
    // find min/max in force profile
    double dFmin = 0;
    double dFmax = 0;
    unsigned int nIndexMin = 0;
    unsigned int nIndexMax = 0;

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        if(_dDensityProfileSmoothed[s] > _dVaporDensity)
        {
            double dF = _dForceProfileSmoothed[s];

            if(dF > dFmax)
            {
                dFmax = dF;
                nIndexMax = s;
            }

            if(dF < dFmin)
            {
                dFmin = dF;
                nIndexMin = s;
            }
        }
    }
*/

    _dInterfaceMidLeft  = _dMidpointPositions[nIndexMax];
    _dInterfaceMidRight = _dMidpointPositions[nIndexMin];
}

void DistControl::EstimateInterfaceMidpoint(Domain* domain)
{
#ifdef ENABLE_MPI

    MPI_Allreduce( _nNumMoleculesLocal, _nNumMoleculesGlobal, _nNumShells, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

#else

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _nNumMoleculesGlobal[s] = _nNumMoleculesLocal[s];
    }

#endif

    // calc density profile
    double dInvertShellVolume = 1. / _dShellVolume;
    double dInvertSampleTimesteps = 1. / ( (double)(_nUpdateFreq) );

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _dDensityProfile[s] = _nNumMoleculesGlobal[s] * dInvertSampleTimesteps * dInvertShellVolume;
    }

/*
    // mheinen_2015-03-17 --> DEBUG_TEST_CASE
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _dDensityProfile[s] = 0.7;
    }

    _dDensityProfile[0] = 0.1;
    _dDensityProfile[1] = 0.1;
    _dDensityProfile[2] = 0.1;
    _dDensityProfile[3] = 0.1;
    _dDensityProfile[4] = 0.2;
    _dDensityProfile[5] = 0.3;
    _dDensityProfile[6] = 0.4;
    _dDensityProfile[7] = 0.5;
    _dDensityProfile[8] = 0.6;

    _dDensityProfile[211] = 0.6;
    _dDensityProfile[212] = 0.5;
    _dDensityProfile[213] = 0.4;
    _dDensityProfile[214] = 0.3;
    _dDensityProfile[215] = 0.2;
    _dDensityProfile[216] = 0.1;
    _dDensityProfile[217] = 0.1;
    _dDensityProfile[218] = 0.1;
    _dDensityProfile[219] = 0.1;


    // result:
    // ----------
    // nIndexSlabGreater = 7
    // x1 = 6.5
    // dx1m = 0
    // _dInterfaceMidLeft = 6.5
    // _dInterfaceMidRight = 213.5

    // <-- DEBUG_TEST_CASE
*/


    // determine left midpoint
    {
        // determine min max
         double dMin = _dDensityProfile[0];
         double dMax = _dDensityProfile[0];

//         int nIndexMin = 0;
//         int nIndexMax = 0;

         for(unsigned int s = 0; s < _nNumShells/2; ++s)
         {
             double dDensity = _dDensityProfile[s];

             // find min
             if(dDensity < dMin)
             {
                 dMin = dDensity;
//                 nIndexMin = s;
             }

             // find max
             if(dDensity > dMax)
             {
                 dMax = dDensity;
//                 nIndexMax = s;
             }
         }

         // dN/dShellWidth = (Nsoll - Nlower) / dInterpol
         //
         // => dInterpol = (Nsoll - Nlower) / dN * dShellWidth
         //
         // dMidPos = dLeftBound + nNumSlabs_lower * dShellWidth + dInterpol

         // 1. ym bestimmen
         double ym = dMin + (dMax - dMin) / 2;

//         cout << "dMin = " << dMin << endl;
//         cout << "dMax = " << dMax << endl;
//         cout << "ym = " << ym << endl;

         // 2. Slab mit gerade niedrigerer Anzahl finden --> von links angefangen ersten Slab größer numMax finden --> Index -1 rechnen
         int nIndexSlabGreater = 0;

         for(unsigned int s = 0; s < _nNumShells; ++s)
         {
             if(_dDensityProfile[s] > ym)
             {
                 nIndexSlabGreater = s;
                 break;
             }
         }

         if(nIndexSlabGreater < 1)
         {
             cout << "ERROR in MEXRegion::EstimateInterfaceMidpoint(): could not find valid index in density profile" << endl;
             return;
         }

//         cout << "nIndexSlabGreater = " << nIndexSlabGreater << endl;

         //3. Interpolierte Distanz dInterpol berechnen, die von dem Slab mit N < Nsoll ausgehend in Richtung N > Nsoll gegangen werden muss
         double dx12, dx1m, y1, y2, dy12, dy1m;

         y1 = _dDensityProfile[nIndexSlabGreater-1];
         y2 = _dDensityProfile[nIndexSlabGreater];

         dy1m = ym - y1;

         dy12 = y2 - y1;
         dx12 = _dShellWidth;

         // dy12 / dx12 = dy1m / dx1m => dx1m = dy1m / dy12 * dx12
         dx1m = dy1m / dy12 * dx12;

         // 4. Berechnung der Mittelpunktsposition
         // nNumSlabslower == nIndexSlabLower + 1 == nIndexSlabGreater
         double xm, x1;
         double dOuterBoundarySampleZone = 0.;  //this->GetLowerCorner()[1];

         x1 = dOuterBoundarySampleZone + nIndexSlabGreater * _dShellWidth - _dShellWidth*0.5;

//         cout << "x1 = " << x1 << endl;
//         cout << "dx1m = " << dx1m << endl;

         xm = x1 + dx1m;
         _dInterfaceMidLeft = xm;

//         cout << "_dInterfaceMidLeft = " << _dInterfaceMidLeft << endl;
     }


    // determine right midpoint
    {
        // determine min max
        double dMin = _dDensityProfile[0];
        double dMax = _dDensityProfile[0];

//         int nIndexMin = 0;
//         int nIndexMax = 0;

         for(unsigned int s = _nNumShells/2; s < _nNumShells; ++s)
         {
             double dDensity = _dDensityProfile[s];

             // find min
             if(dDensity < dMin)
             {
                 dMin = dDensity;
//                 nIndexMin = s;
             }

             // find max
             if(dDensity > dMax)
             {
                 dMax = dDensity;
//                 nIndexMax = s;
             }
         }

         // dN/dShellWidth = (Nsoll - Nlower) / dInterpol
         //
         // => dInterpol = (Nsoll - Nlower) / dN * dShellWidth
         //
         // dMidPos = dLeftBound + nNumSlabs_lower * dShellWidth + dInterpol

         // 1. ym bestimmen
         double ym = dMin + (dMax - dMin) / 2;

//         cout << "ym = " << ym << endl;

         // 2. Slab mit gerade niedrigerer Anzahl finden --> von links angefangen ersten Slab größer numMax finden --> Index -1 rechnen
         unsigned int nIndexSlabGreater = 0;

         for(int s = _nNumShells-1; s >= 0; --s)
         {
             if(_dDensityProfile[s] > ym)
             {
                 nIndexSlabGreater = s;
                 break;
             }
         }

         if(nIndexSlabGreater < 1)
         {
             cout << "ERROR in MEXRegion::EstimateInterfaceMidpoint(): could not find valid index in density profile" << endl;
             return;
         }

         //3. Interpolierte Distanz dInterpol berechnen, die von dem Slab mit N < Nsoll ausgehend in Richtung N > Nsoll gegangen werden muss
         double dx12, dx1m, y1, y2, dy12, dy1m;

         y1 = _dDensityProfile[nIndexSlabGreater+1];
         y2 = _dDensityProfile[nIndexSlabGreater];

         dy1m = ym - y1;

         dy12 = y2 - y1;
         dx12 = _dShellWidth;

         // dy12 / dx12 = dy1m / dx1m => dx1m = dy1m / dy12 * dx12
         dx1m = dy1m / dy12 * dx12;

         // 4. Berechnung der Mittelpunktsposition
         // nNumSlabslower == nIndexSlabLower + 1 == nIndexSlabGreater
         double xm, x1;
         double dOuterBoundarySampleZone = domain->getGlobalLength(1);  //this->GetUpperCorner()[1];
         unsigned int numShellsLower = _nNumShells-1 - nIndexSlabGreater;

//         cout << "numShellsLower = " << numShellsLower << endl;

         x1 = dOuterBoundarySampleZone - numShellsLower * _dShellWidth + _dShellWidth*0.5;

//         cout << "x1 = " << x1 << endl;
//         cout << "dx1m = " << dx1m << endl;

         xm = x1 - dx1m;
         _dInterfaceMidRight = xm;

//         cout << "_dInterfaceMidRight = " << _dInterfaceMidRight << endl;
    }
}


// init
void DistControl::InitPositions(double dInterfaceMidLeft, double dInterfaceMidRight)
{
    // set interface midpoints manually
    _dInterfaceMidLeft  = dInterfaceMidLeft;
    _dInterfaceMidRight = dInterfaceMidRight;

    // control volume (CV)
    _dControlVolumeLeft  = _dInterfaceMidLeft  - _dCVFactor * _d1090Thickness;
    _dControlVolumeRight = _dInterfaceMidRight + _dCVFactor * _d1090Thickness;

    // thermostating zone
    _dTZoneLeft  = _dInterfaceMidLeft  + _dTZoneFactor * _d1090Thickness;
    _dTZoneRight = _dInterfaceMidRight - _dTZoneFactor * _d1090Thickness;

    // sampling zone
    _dSZoneLeft_uc  = _dInterfaceMidLeft  + _dSZoneFactor * _d1090Thickness;
    _dSZoneRight_lc = _dInterfaceMidRight - _dSZoneFactor * _d1090Thickness;

    _dSZoneLeft_lc  = _dControlVolumeLeft  - _dCVWidth;
    _dSZoneRight_uc = _dControlVolumeRight + _dCVWidth;

    // drift control region
    _dDriftControlLeft  = _dInterfaceMidLeft  + _dTZoneFactor * _d1090Thickness;
    _dDriftControlRight = _dInterfaceMidRight - _dTZoneFactor * _d1090Thickness;
}

void DistControl::UpdatePositions(unsigned long simstep, Domain* domain)
{
    // update with respect to update frequency
    if(simstep % _nUpdateFreq != 0 || simstep == 0 || simstep == 1)  // TODO init timestep
        return;

    // TODO: check for initial timestep???

    if(_nMethod == 1)
    	this->EstimateInterfaceMidpoint(domain);
    else if (_nMethod == 2)
    	this->EstimateInterfaceMidpointsByForce();
    else
    {
    	cout << "DistControl::UpdatePositions: Corrupted code!!! Programm exit..." << endl;
    	exit(-1);
    }

    // update positions

    if(_nSimType == SIMTYPE_EVAPORATION)
    {
        // control volume (CV)
        _dControlVolumeLeft  = _dInterfaceMidLeft  - _dCVFactor * _d1090Thickness;
        _dControlVolumeRight = _dInterfaceMidRight + _dCVFactor * _d1090Thickness;

        // thermostating zone
        _dTZoneLeft  = _dInterfaceMidLeft  + _dTZoneFactor * _d1090Thickness;
        _dTZoneRight = _dInterfaceMidRight - _dTZoneFactor * _d1090Thickness;

        // sampling zone
        _dSZoneLeft_uc  = _dInterfaceMidLeft  + _dSZoneFactor * _d1090Thickness;
        _dSZoneRight_lc = _dInterfaceMidRight - _dSZoneFactor * _d1090Thickness;

        _dSZoneLeft_lc  = _dControlVolumeLeft  - _dCVWidth;
        _dSZoneRight_uc = _dControlVolumeRight + _dCVWidth;
    }

    // drift control region
    _dDriftControlLeft  = _dInterfaceMidLeft  + _dTZoneFactor * _d1090Thickness;
    _dDriftControlRight = _dInterfaceMidRight - _dTZoneFactor * _d1090Thickness;

    // DEBUG
#ifdef DEBUG
    cout << "_dInterfaceMidLeft  = " << _dInterfaceMidLeft << endl;
    cout << "_dInterfaceMidRight = " << _dInterfaceMidRight << endl;


    cout << "_dVacuumLeft      = " << _dVacuumLeft << endl;
    cout << "_dVacuumRight     = " << _dVacuumRight << endl;
    cout << "_dTzoneLeft       = " << _dTzoneLeft << endl;
    cout << "_dTzoneRight      = " << _dTzoneRight << endl;
    cout << "_dLeftSzoneLeft   = " << _dLeftSzoneLeft << endl;
    cout << "_dLeftSzoneRight  = " << _dLeftSzoneRight << endl;
    cout << "_dRightSzoneLeft  = " << _dRightSzoneLeft << endl;
    cout << "_dRightSzoneRight = " << _dRightSzoneRight << endl;
    cout << "_dFluxAreaLeft    = " << _dFluxAreaLeft << endl;
    cout << "_dFluxAreaRight   = " << _dFluxAreaRight << endl;
#endif

    // reset local values
    this->ResetLocalValues();
}


void DistControl::ResetLocalValues()
{
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _nNumMoleculesLocal[s] = 0;
        _dForceSumLocal[s] = 0.;
    }
}

void DistControl::WriteData(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{
    // write out data
    stringstream outputstream;

    // write data
    if(simstep % _nWriteFreq == 0)
    {
        // calc global values
//        this->CalcGlobalValues(domainDecomp);

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif
        outputstream << std::setw(20) << simstep;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dInterfaceMidLeft;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dInterfaceMidRight;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dControlVolumeLeft;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dControlVolumeRight;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dTZoneLeft;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dTZoneRight;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dSZoneLeft_lc;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dSZoneLeft_uc;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dSZoneRight_lc;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dSZoneRight_uc;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dDriftControlLeft;
        outputstream << std::setw(16) << fixed << std::setprecision(7) << _dDriftControlRight;
        outputstream << endl;

        ofstream fileout(_strFilename.c_str(), ios::out|ios::app);
        fileout << outputstream.str();
        fileout.close();
#ifdef ENABLE_MPI
    }
#endif
    }
}

void DistControl::WriteHeader(DomainDecompBase* domainDecomp, Domain* domain)
{
    // write header
    stringstream outputstream;

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif

    outputstream << "             simstep";
    outputstream << "         midLeft" << "        midRight";
    outputstream << "          cvLeft" << "         cvRight";
    outputstream << "          tzLeft" << "         tzRight";
    outputstream << "       szLeft_lc" << "       szLeft_uc";
    outputstream << "      szRight_lc" << "      szRight_uc";
    outputstream << "          dcLeft" << "         dcRight";
    outputstream << endl;

    ofstream fileout(_strFilename.c_str(), ios::out);
    fileout << outputstream.str();
    fileout.close();

#ifdef ENABLE_MPI
    }
#endif
}


void DistControl::WriteDataProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep)
{
    // write out data
    stringstream outputstream;
    stringstream filenamestream;

    filenamestream << _strFilenameProfilesPrefix << "_TS";
    filenamestream.fill('0');
    filenamestream.width(9);
    filenamestream << right << simstep;
    filenamestream << ".dat";

    // write data
    if(simstep % _nWriteFreqProfiles == 0)
    {
        // calc global values
//        this->CalcGlobalValues(domainDecomp);

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif

        // write header
        outputstream << "           y" << "         rho" << "  rho_smooth" << "     drho/dy" << "          Fy" << "   Fy_smooth" << endl;

        for(unsigned int s=0; s<_nNumShells; s++)
        {
            outputstream << std::setw(12) << fixed << std::setprecision(3) << _dMidpointPositions[s];
            outputstream << std::setw(12) << fixed << std::setprecision(3) << _dDensityProfile[s];
            outputstream << std::setw(12) << fixed << std::setprecision(3) << _dDensityProfileSmoothed[s];
            outputstream << std::setw(12) << fixed << std::setprecision(3) << _dDensityProfileSmoothedDerivation[s];
            outputstream << std::setw(12) << fixed << std::setprecision(3) << _dForceProfile[s];
            outputstream << std::setw(12) << fixed << std::setprecision(3) << _dForceProfileSmoothed[s];
            outputstream << endl;
        }

        ofstream fileout(filenamestream.str().c_str(), ios::out|ios::out);
        fileout << outputstream.str();
        fileout.close();
#ifdef ENABLE_MPI
    }
#endif
    }
}


void DistControl::Init(DomainDecompBase* domainDecomp, Domain* domain, ParticleContainer* particleContainer)
{
    // write output file with header
    this->WriteHeader(domainDecomp, domain);
}


void DistControl::AlignSystemCenterOfMass(Domain* domain, Molecule* mol, unsigned long simstep)
{
    // update with respect to update frequency
    if(simstep % _nUpdateFreq != 0 || simstep == 0 || simstep == 1)  // TODO init timestep
        return;

    double dDeltaLeftY  = _dInterfaceMidLeft;  // (minus 0)
    double dDeltaRightY = domain->getGlobalLength(1) - _dInterfaceMidRight;
    double dDeltaY = 0.5 * (dDeltaRightY - dDeltaLeftY);

    double dNewPosition = mol->r(1) + dDeltaY;

/*
 * seems not to work: simulation crashes!!!
 *
    double dBoxLengthY = domain->getGlobalLength(1);

    if (dNewPosition > dBoxLengthY )
        dNewPosition -= dBoxLengthY;
    else if (dNewPosition < 0. )
        dNewPosition += dBoxLengthY;
*/

    mol->setr(1, dNewPosition);

#ifdef DEBUG
    cout << "_dInterfaceMidLeft = " << _dInterfaceMidLeft << endl;
    cout << "_dInterfaceMidRight = " << _dInterfaceMidRight << endl;
    cout << "dDeltaLeftY = " << dDeltaLeftY << endl;
    cout << "dDeltaRightY = " << dDeltaRightY << endl;
    cout << "dDeltaY = " << dDeltaY << endl;
#endif
}





