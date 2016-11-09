/*
 * DistControl.h
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#ifndef DISTCONTROL_H_
#define DISTCONTROL_H_

#include <string>
#include <vector>
#include "utils/ObserverBase.h"

using namespace std;

class Domain;
class ParticleContainer;
class DomainDecompBase;
class Molecule;

class DistControl : public SubjectBase
{
public:
    DistControl(Domain* domain, unsigned int nUpdateFreq, unsigned int nNumShells, double dVaporDensity, unsigned int nMethod);
    ~DistControl();

    // init
    void InitPositions(double dInterfaceMidLeft, double dInterfaceMidRight);

    double GetInterfaceMidLeft() {return _dInterfaceMidLeft;}
    double GetInterfaceMidRight() {return _dInterfaceMidRight;}
    void SetWriteFreqProfiles(unsigned int nVal) {_nWriteFreqProfiles = nVal;}
    unsigned int GetUpdateFreq() {return _nUpdateFreq;}

    void Init(DomainDecompBase* domainDecomp, Domain* domain, ParticleContainer* particleContainer);
    void WriteHeader(DomainDecompBase* domainDecomp, Domain* domain);
    void WriteData(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
    void WriteDataProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);


    // place method inside loop over molecule container
    void SampleProfiles(Molecule* mol);

    void UpdatePositions(unsigned long simstep, Domain* domain);
    void AlignSystemCenterOfMass(Domain* domain, Molecule* mol, unsigned long simstep);

    // SubjectBase methods
	virtual void registerObserver(ObserverBase* observer);
	virtual void deregisterObserver(ObserverBase* observer);
	virtual void informObserver();

private:
    // place methods after the loop
    void EstimateInterfaceMidpoint(Domain* domain);  // called by UpdatePositions
    void EstimateInterfaceMidpointsByForce();
    void ResetLocalValues();

private:
    double _dInterfaceMidLeft;
    double _dInterfaceMidRight;

    unsigned long* _nNumMoleculesLocal;
    unsigned long* _nNumMoleculesGlobal;
    double* _dMidpointPositions;
	double* _dForceSumLocal;
	double* _dForceSumGlobal;
    double* _dDensityProfile;
    double* _dDensityProfileSmoothed;
    double* _dDensityProfileSmoothedDerivation;
	double* _dForceProfile;
    double* _dForceProfileSmoothed;
    unsigned int _nNumShells;
    double _dShellWidth;
    double _dInvertShellWidth;
    double _dShellVolume;
    unsigned int _nUpdateFreq;
    unsigned int _nWriteFreq;
    unsigned int _nWriteFreqProfiles;
    double _dVaporDensity;
    unsigned int _nMethod;

    // write data
    string _strFilename;
    string _strFilenameProfilesPrefix;

    // simtype
    int _nSimType;

    // observer
	std::vector<ObserverBase*> _observer;

};  // class DistControl


#endif /* DISTCONTROL_H_ */
