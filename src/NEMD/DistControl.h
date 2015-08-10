/*
 * DistControl.h
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#ifndef DISTCONTROL_H_
#define DISTCONTROL_H_

#include <string>

using namespace std;

class Domain;
class ParticleContainer;
class DomainDecompBase;
class Molecule;

class DistControl
{
public:
    DistControl(Domain* domain, unsigned int nUpdateFreq, unsigned int nNumShells);
    ~DistControl();

    // init
    void InitPositions(double dInterfaceMidLeft, double dInterfaceMidRight);

    double GetInterfaceMidLeft() {return _dInterfaceMidLeft;}
    double GetInterfaceMidRight() {return _dInterfaceMidRight;}
    void SetDistanceParameters(double d1090Thickness, double dCVFactor, double dTZoneFactor, double dSZoneFactor, double dCVWidth)
    {
        _d1090Thickness = d1090Thickness;
        _dCVFactor      = dCVFactor;
        _dTZoneFactor   = dTZoneFactor;
        _dSZoneFactor   = dSZoneFactor;
        _dCVWidth   = dCVWidth;
    }

    // get positions
    // positions
    double GetCVLeft() {return _dControlVolumeLeft;}
    double GetCVRight() {return _dControlVolumeRight;}
    double GetTZoneLeft()  {return _dTZoneLeft;}
    double GetTZoneRight() {return _dTZoneRight;}
    double GetSZoneLeft_lc()  {return _dSZoneLeft_lc;}
    double GetSZoneLeft_uc()  {return _dSZoneLeft_uc;}
    double GetSZoneRight_lc() {return _dSZoneRight_lc;}
    double GetSZoneRight_uc() {return _dSZoneRight_uc;}

    void Init(DomainDecompBase* domainDecomp, Domain* domain, ParticleContainer* particleContainer);
    void WriteHeader(DomainDecompBase* domainDecomp, Domain* domain);
    void WriteData(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
    void WriteDataDensity(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);


    // place method inside loop over molecule container
    void SampleDensityProfile(Molecule* mol);

    // place methods after the loop
private:
    void EstimateInterfaceMidpoint(Domain* domain);  // called by UpdatePositions
public:
    void UpdatePositions(unsigned long simstep, Domain* domain);
    void AlignSystemCenterOfMass(Domain* domain, Molecule* mol, unsigned long simstep);

private:
    void ResetLocalValues();

private:
    double _dInterfaceMidLeft;
    double _dInterfaceMidRight;

    unsigned long* _nNumMoleculesLocal;
    unsigned long* _nNumMoleculesGlobal;
    double* _dDensityProfile;
    unsigned int _nNumShells;
    double _dShellWidth;
    double _dInvertShellWidth;
    double _dShellVolume;
    unsigned int _nUpdateFreq;
    unsigned int _nWriteFreq;
    unsigned int _nWriteFreqDensity;

    // distance parameters
    // all distances are related to the 10-90 thickness of the interface
    double _d1090Thickness;
    double _dCVFactor;
    double _dTZoneFactor;
    double _dSZoneFactor;
    double _dCVWidth;

    // positions
    double _dControlVolumeLeft;
    double _dControlVolumeRight;
    double _dTZoneLeft;
    double _dTZoneRight;
    double _dSZoneLeft_lc;
    double _dSZoneLeft_uc;
    double _dSZoneRight_lc;
    double _dSZoneRight_uc;

    // write data
    string _strFilename;
    string _strFilenameDensityPrefix;

};  // class DistControl


#endif /* DISTCONTROL_H_ */
