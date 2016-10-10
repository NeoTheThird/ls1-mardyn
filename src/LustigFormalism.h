/*
 * LustigFormalism.h
 *
 *  Created on: 18.05.2016
 *      Author: mheinen
 */

#ifndef LUSTIGFORMALISM_H_
#define LUSTIGFORMALISM_H_

#include "molecules/Comp2Param.h"
#include <string>

class Domain;
class DomainDecompBase;

class LustigFormalism
{
public:
	LustigFormalism();
	~LustigFormalism();

	void SetWriteFreq(unsigned long nWriteFreq, unsigned long nStart, unsigned long nStop, unsigned long nWriteFreqSums)
	{_nWriteFreq = nWriteFreq; _nStart = nStart; _nStop = nStop; _nWriteFreqSums = nWriteFreqSums;}
	void InitRestart(std::string strRestartFilenameSums, unsigned long nRestartTimestep);
	unsigned long GetStart() {return _nStart;}
	unsigned long GetStop()  {return _nStop;}
	void InitNVT(Domain* domain, unsigned long N, double V, double T, double cutoffRadiusLJ);
	void Init(const double& U6, const double& dUdV, const double& d2UdV2);
	void InitWidom(const double& DU, const double& T);
	void InitMSD(const unsigned long& nID, double const * dr);
	void CalcGlobalValues(DomainDecompBase* domainDecomp);
    void CalcDerivatives();
	void WriteHeader(DomainDecompBase* domainDecomp, Domain* domain);
	void WriteData(DomainDecompBase* domainDecomp, unsigned long simstep);
    void TriggerNewSimstep(unsigned long simstep){_bSimstepTrigger = true; _simstep = simstep;}
    void ResetTriggerNewSimstep(){_bSimstepTrigger = false;}
    bool IsNewSimstep(){return _bSimstepTrigger;}
    void SetNumWidomTests(unsigned int nNumWidomTests){_nNumWidomTestsGlobal = nNumWidomTests;}
    void EventConfigurationSampled();
    void StoreDomainDecompPointer(DomainDecompBase* domainDecomp) {_domainDecomp = domainDecomp;}
    void StoreTimestepLength(double dTimestepLength) { _Dt_ts = dTimestepLength;}
    void StartSamplingMSD(){_bSampleMSD = true;}
    bool SamplingStartedMSD(){return _bSampleMSD;}

private:
    // reset local values
    void ResetSums();
    void ResetLocalValues();
    void InitDatastructures();
    void InitSums();


private:
	unsigned long _N;
	double _InvN;
	double _V;
    double _InvV;
	double _mInvV;
	double _InvV2;
	double _T;
	unsigned long _tsBufferIndex;
	unsigned long _nWriteFreq;
	unsigned long _nWriteFreqSums;
	unsigned long _nStart;
	unsigned long _nStop;
	unsigned long _nNumConfigs;
    double _rho;
	double _v;
	double _v2;
	double _beta;
	double _beta2;
	double _beta3;

	// local
	double* _ULocal;
	double* _dUdVLocal;
	double* _d2UdV2Local;

    // global
	double* _UGlobal;
	double* _dUdVGlobal;
	double* _d2UdV2Global;

	double* _U2Global;
    double* _U3Global;
    double* _dUdV2Global;
    double* _UdUdVGlobal;
    double* _U2dUdVGlobal;
    double* _UdUdV2Global;
    double* _Ud2UdV2Global;

    // sums
    double _UGlobalSum;
    double _U2GlobalSum;
    double _U3GlobalSum;
    double _dUdVGlobalSum;
    double _d2UdV2GlobalSum;
    double _dUdV2GlobalSum;
    double _UdUdVGlobalSum;
    double _U2dUdVGlobalSum;
    double _UdUdV2GlobalSum;
    double _Ud2UdV2GlobalSum;

    // LRC
    double _U_LRC;
    double _dUdV_LRC;
    double _d2UdV2_LRC;
    double _dU_LRC;

	// ideal
	double _A00i;
	double _A10i;
	double _A01i;
	double _A20i;
	double _A11i;
	double _A02i;
	double _A30i;
	double _A21i;
	double _A12i;

	// residual
	double _A00r;
	double _A10r;
	double _A01r;
	double _A20r;
	double _A11r;
	double _A02r;
	double _A30r;
	double _A21r;
	double _A12r;
	// thermodynamic properties
	double _Cv;

	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;

	// mu
    double* _WidomEnergyLocal;
    double* _WidomEnergyGlobal;
    double _WidomEnergyGlobalSum;
    unsigned int _nNumWidomTestsLocal;
    unsigned int _nNumWidomTestsGlobal;
    double _mu_res;

    bool _bSimstepTrigger;
    unsigned long _simstep;

    DomainDecompBase* _domainDecomp;

    // restart
    std::string _strRestartFilenameSums;
    unsigned long _nRestartTimestep;

    // MSD
    //double _D;           // self-diffusion coefficient
    double _MSD;          // Mean Squared Displacement
    double _Dt_ts;       // time step length
    double _dInv6NDt_ts;  // 1 / (6*N*Dt_ts)
    double** _dDisplacementVecLocal;
    double** _dDisplacementVecGlobal;
    double _dDisplacementGlobalSum;
    bool _bSampleMSD;

};  // class LustigFormalism



#endif /* LUSTIGFORMALISM_H_ */
