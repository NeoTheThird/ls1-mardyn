/*
 * LustigFormalism.cpp
 *
 *  Created on: 18.05.2016
 *      Author: mheinen
 */

#include "LustigFormalism.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cmath>

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/ParaStrm.h"

using namespace std;

//typedef std::numeric_limits<double> dbl;

LustigFormalism::LustigFormalism()
{
	_tsBufferIndex = 0;
	_nWriteFreq = 1000;
	_nWriteFreqSums = 100000;

	// reset sums
	this->ResetSums();

	// init num widom tests
	_nNumWidomTestsGlobal = 0;

	// allocate memory
	_ULocal       = new double[_nWriteFreq];
	_dUdVLocal    = new double[_nWriteFreq];
	_d2UdV2Local  = new double[_nWriteFreq];
	_UGlobal      = new double[_nWriteFreq];
	_dUdVGlobal   = new double[_nWriteFreq];
	_d2UdV2Global = new double[_nWriteFreq];

	_U2Global      = new double[_nWriteFreq];
    _U3Global      = new double[_nWriteFreq];
    _dUdV2Global   = new double[_nWriteFreq];
    _UdUdVGlobal   = new double[_nWriteFreq];
    _U2dUdVGlobal  = new double[_nWriteFreq];
    _UdUdV2Global  = new double[_nWriteFreq];
    _Ud2UdV2Global = new double[_nWriteFreq];

    _WidomEnergyLocal  = new double[_nWriteFreq];
    _WidomEnergyGlobal = new double[_nWriteFreq];

	// reset (init) local values
	this->ResetLocalValues();

	_bSimstepTrigger = false;
}

LustigFormalism::~LustigFormalism()
{
	// free memory
	delete[] _ULocal;
	delete[] _dUdVLocal;
	delete[] _d2UdV2Local;
	delete[] _UGlobal;
	delete[] _dUdVGlobal;
	delete[] _d2UdV2Global;

	delete[] _U2Global;
	delete[] _U3Global;
	delete[] _dUdV2Global;
	delete[] _UdUdVGlobal;
	delete[] _U2dUdVGlobal;
	delete[] _UdUdV2Global;
	delete[] _Ud2UdV2Global;

	delete[] _WidomEnergyLocal;
	delete[] _WidomEnergyGlobal;
}

void LustigFormalism::InitSums(string strFilename)
{
	ifstream filein(strFilename.c_str(), ios::in);

	string zeile1, zeile2;
	while (getline (filein, zeile1)) {
	  zeile2.swap (zeile1);
	}

	if (filein.bad () || !filein.eof ()) {
	  // Datei konnte nicht gelesen werden o. sonst. Fehler
		cout << "Couldnt read file!" << endl;
	} else {
	  // Variable zeile2 enthält hier die letzte Zeile der Datei.
		char * cstr = new char [zeile2.length()+1];
		std::strcpy (cstr, zeile2.c_str());
		char* pch;
		string str[20];

		pch = strtok(cstr, " ");
		str[0] = pch;

		int i=0;
		while (pch != NULL)
		{
			str[i] = string(pch);
#ifdef DEBUG
			cout << "str[" << i << "] = " << str[i] << endl;
			cout << "i = " << i << endl;
#endif
			i++;
			pch = strtok (NULL, " ");
		}

		delete[] cstr;

		_nNumConfigs      = atoi(str[0].c_str());
		_UGlobalSum       = atof(str[1].c_str());
		_U2GlobalSum      = atof(str[2].c_str());
		_U3GlobalSum      = atof(str[3].c_str());
		_dUdVGlobalSum    = atof(str[4].c_str());
		_d2UdV2GlobalSum  = atof(str[5].c_str());
		_dUdV2GlobalSum   = atof(str[6].c_str());
		_UdUdVGlobalSum   = atof(str[7].c_str());
		_U2dUdVGlobalSum  = atof(str[8].c_str());
		_UdUdV2GlobalSum  = atof(str[9].c_str());
		_Ud2UdV2GlobalSum = atof(str[10].c_str());
		_WidomEnergyGlobalSum = atof(str[11].c_str()) * _nNumWidomTestsGlobal;

#ifdef DEBUG
		cout << "_nNumConfigs      = " << _nNumConfigs << endl;
		cout << "_UGlobalSum       = " << _UGlobalSum << endl;
		cout << "_U2GlobalSum      = " << _U2GlobalSum << endl;
		cout << "_U3GlobalSum      = " << _U3GlobalSum << endl;
		cout << "_dUdVGlobalSum    = " << _dUdVGlobalSum << endl;
		cout << "_d2UdV2GlobalSum  = " << _d2UdV2GlobalSum << endl;
		cout << "_dUdV2GlobalSum   = " << _dUdV2GlobalSum << endl;
		cout << "_UdUdVGlobalSum   = " << _UdUdVGlobalSum << endl;
		cout << "_U2dUdVGlobalSum  = " << _U2dUdVGlobalSum << endl;
		cout << "_UdUdV2GlobalSum  = " << _UdUdV2GlobalSum << endl;
		cout << "_Ud2UdV2GlobalSum = " << _Ud2UdV2GlobalSum << endl;
		cout << "_WidomEnergyGlobalSum = " << _WidomEnergyGlobalSum << endl;
#endif
	}

	filein.close();
}

void LustigFormalism::InitNVT(Domain* domain, unsigned long N, double V, double T, double cutoffRadiusLJ)
{
	_N = N;
	_V = V;
	_T = T;

#ifdef DEBUG
	cout << ">>> Lustig formalism <<<" << endl;
	cout << "N = " << N << endl;
	cout << "V = " << V << endl;
	cout << "T = " << T << endl;
#endif

    // temperature
	_beta = 1./_T;
	_beta2 = _beta*_beta;
	_beta3 = _beta*_beta2;

    // density
    _rho = _N/_V;
	_v   = 1./_rho;
	_v2  = _v*_v;

#ifdef DEBUG
	cout << "_rho = " << _rho << endl;
	cout << "_v = " << _v << endl;
	cout << "_beta = " << _beta << endl;
#endif

	_InvN = 1./(double)(_N);

	// invert volume
    _InvV  = 1./_V;
	_mInvV = -1.*_InvV;
	_InvV2 = _InvV*_InvV;
    
    // LRC
    // TODO: get these values from elsewhere 
	_comp2params = domain->getComp2Params();
	ParaStrm& params = _comp2params(0, 0);
	params.reset_read();
	double eps24;
	params >> eps24;
	double sig2;
	params >> sig2;
	double uLJshift6;
	params >> uLJshift6;  // 0 unless TRUNCATED_SHIFTED

    double rc = cutoffRadiusLJ;
    double eps = eps24 / 24.;

#ifdef DEBUG
    cout << "cutoffRadiusLJ = " << cutoffRadiusLJ << endl;
    cout << "eps24 = " << eps24 << endl;
    cout << "sig2 = " << sig2 << endl;
#endif

    double rc2 = rc*rc;
    double rc3 = rc*rc2;
    double invrc2 = 1. / rc2;
    double lj6 = sig2 * invrc2; lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    
    const double PI = 3.14159265358979323846;
	if(uLJshift6 == 0.)
		_U_LRC = PI*_rho*eps24*rc3*(1./3.*lj12 - lj6) * _N/9.;
	else
		_U_LRC = 0.;
    _dUdV_LRC   = -8.*PI/9.*_InvV *(_N-1.)*_rho*rc3*eps*( 4.*lj12 -  6.*lj6);
    _d2UdV2_LRC =  8.*PI/9.*_InvV2*(_N-1.)*_rho*rc3*eps*(20.*lj12 - 18.*lj6);

#ifdef DEBUG
    cout << "_U_LRC = " << _U_LRC << endl;
    cout << "_dUdV_LRC = " << _dUdV_LRC << endl;
    cout << "_d2UdV2_LRC = " << _d2UdV2_LRC << endl;
#endif
}

void LustigFormalism::Init(const double& U6, const double& dUdVm3, const double& d2UdV2m3)
{
	_ULocal[_tsBufferIndex]      = U6 / 6.;
	_dUdVLocal[_tsBufferIndex]   = _mInvV * dUdVm3 / 3.;
	_d2UdV2Local[_tsBufferIndex] = _InvV2 * d2UdV2m3 / 3.;
}

void LustigFormalism::EventConfigurationSampled()
{
	_tsBufferIndex++;
	_nNumConfigs++;
}

void LustigFormalism::InitWidom(const double& dU)
{
	_WidomEnergyLocal[_tsBufferIndex] += dU;
}

void LustigFormalism::CalcGlobalValues(DomainDecompBase* domainDecomp)
{

#ifdef DEBUG

#ifdef ENABLE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	unsigned long tsBufferIndexGlobal;
    domainDecomp->collCommInit(1);
    domainDecomp->collCommAppendUnsLong(_tsBufferIndex);
    domainDecomp->collCommAllreduceSum();
    tsBufferIndexGlobal = domainDecomp->collCommGetUnsLong();
    domainDecomp->collCommFinalize();
    cout << "tsBufferIndexGlobal = " << tsBufferIndexGlobal << endl;
#endif

    // Collective Communiation
    domainDecomp->collCommInit(4*_nWriteFreq);

    for(unsigned int i=0; i<_nWriteFreq; ++i)
    {
        domainDecomp->collCommAppendDouble(_ULocal[i]);
        domainDecomp->collCommAppendDouble(_dUdVLocal[i]);
        domainDecomp->collCommAppendDouble(_d2UdV2Local[i]);
        domainDecomp->collCommAppendDouble(_WidomEnergyLocal[i]);
    }
    domainDecomp->collCommAllreduceSum();

    for(unsigned int i=0; i<_nWriteFreq; ++i)
    {
        _UGlobal[i]      = domainDecomp->collCommGetDouble();
        _dUdVGlobal[i]   = domainDecomp->collCommGetDouble();
        _d2UdV2Global[i] = domainDecomp->collCommGetDouble();
        _WidomEnergyGlobal[i] = domainDecomp->collCommGetDouble();
    }
    domainDecomp->collCommFinalize();

	// reset local values
	this->ResetLocalValues();

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank != 0)
		return;
#endif

	for(unsigned int i=0; i<_nWriteFreq; ++i)
	{
		// LRC
		_UGlobal[i]      += _U_LRC;
		_dUdVGlobal[i]   += _dUdV_LRC;
		_d2UdV2Global[i] += _d2UdV2_LRC;

		// squared values, 3rd potenz
		_U2Global[i] = _UGlobal[i]*_UGlobal[i];
	    _U3Global[i] = _UGlobal[i]*_U2Global[i];
	    _dUdV2Global[i] = _dUdVGlobal[i]*_dUdVGlobal[i];

	    // mixed values
	    _UdUdVGlobal[i]     = _UGlobal[i]  * _dUdVGlobal[i];
	    _U2dUdVGlobal[i]    = _U2Global[i] * _dUdVGlobal[i];
	    _Ud2UdV2Global[i]   = _UGlobal[i]  * _d2UdV2Global[i];
	    _UdUdV2Global[i]    = _UGlobal[i]  * _dUdV2Global[i];

	    // accumulate
		_UGlobalSum       += _UGlobal[i];
		_U2GlobalSum      += _U2Global[i];
	    _U3GlobalSum      += _U3Global[i];
		_dUdVGlobalSum    += _dUdVGlobal[i];
		_d2UdV2GlobalSum  += _d2UdV2Global[i];
	    _dUdV2GlobalSum   += _dUdV2Global[i];
	    _UdUdVGlobalSum   += _UdUdVGlobal[i];
	    _U2dUdVGlobalSum  += _U2dUdVGlobal[i];
	    _UdUdV2GlobalSum  += _UdUdV2Global[i];
	    _Ud2UdV2GlobalSum += _Ud2UdV2Global[i];

	    // chem. potential mu
	    _WidomEnergyGlobalSum += _WidomEnergyGlobal[i];
	}
}

void LustigFormalism::CalcDerivatives()
{
	// divide by number of sampled configurations
	double InvNumConfigs = 1. / (double)(_nNumConfigs);

	double U       = _UGlobalSum       * InvNumConfigs;
	double U2      = _U2GlobalSum      * InvNumConfigs;
    double U3      = _U3GlobalSum      * InvNumConfigs;
	double dUdV    = _dUdVGlobalSum    * InvNumConfigs;
	double d2UdV2  = _d2UdV2GlobalSum  * InvNumConfigs;
    double dUdV2   = _dUdV2GlobalSum   * InvNumConfigs;
    double UdUdV   = _UdUdVGlobalSum   * InvNumConfigs;
    double U2dUdV  = _U2dUdVGlobalSum  * InvNumConfigs;
    double UdUdV2  = _UdUdV2GlobalSum  * InvNumConfigs;
    double Ud2UdV2 = _Ud2UdV2GlobalSum * InvNumConfigs;

//    cout << "U = " <<  U << endl;
//    cout << "U2 = " <<  U2 << endl;
//    cout << "U3 = " <<  U3 << endl;
//    cout << "dUdV = " <<  dUdV << endl;
//    cout << "d2UdV2 = " <<  d2UdV2 << endl;
//    cout << "dUdV2 = " <<  dUdV2 << endl;
//    cout << "UdUdV = " <<  UdUdV << endl;
//    cout << "U2dUdV = " <<  U2dUdV << endl;
//    cout << "UdUdV2 = " <<  UdUdV2 << endl;
//    cout << "Ud2UdV2 = " <<  Ud2UdV2 << endl;

	// derivatives
	_A10r = _beta*U*_InvN;
	_A01r = -1.*_beta*_v*dUdV;
	_A20r = _beta2*_InvN*(U*U - U2);
	_A11r = -1.*_v*_beta*dUdV + _v*_beta2*UdUdV - _v*_beta2*U*dUdV;
	_A02r = _v2*_N*_beta*d2UdV2 -_v2*_N*_beta2*dUdV2 + _v2*_N*_beta2*dUdV*dUdV + 2.*_v*_beta*dUdV;
	_A30r = _beta3*_InvN*(U3 - 3.*U*U2 + 2.*U*U*U);
	_A21r = 2.*_v*_beta2*UdUdV - 2.*_v*_beta2*U*dUdV + _v*_beta3*U2*dUdV -_v*_beta3*U2dUdV + 2.*_v*_beta3*U*UdUdV - 2.*_v*_beta3*U*U*dUdV;

	_A12r  =    _v2*_N*_beta3*UdUdV2    + 2.*_v2*_N*_beta3*U*dUdV*dUdV -    _v2*_N*_beta3*U*dUdV2 - 2.*_v2*_N*_beta3*UdUdV*dUdV;
	_A12r += 2.*_v2*_N*_beta2*dUdV*dUdV +    _v2*_N*_beta2*U*d2UdV2    - 2.*_v2*_N*_beta2*dUdV2   -    _v2*_N*_beta2*Ud2UdV2;
	_A12r +=    _v2*_N*_beta *d2UdV2;
	_A12r += 2.*_v *   _beta2*U*dUdV    - 2.*_v    *_beta2*UdUdV;
	_A12r += 2.*_v *   _beta *dUdV;

	// chem. potential
	if(_nNumWidomTestsGlobal < 1)
		_mu_res = NAN;
	else
	{
		double dInvNumWidomTests = 1. / (double)(_nNumWidomTestsGlobal * _nNumConfigs);
		_mu_res = -log(_WidomEnergyGlobalSum * dInvNumWidomTests);
	}

	// thermodynamic properties
}

void LustigFormalism::ResetSums()
{
	_nNumConfigs = 0;

	_UGlobalSum       = 0.;
	_U2GlobalSum      = 0.;
	_U3GlobalSum      = 0.;
	_dUdVGlobalSum    = 0.;
	_d2UdV2GlobalSum  = 0.;
    _dUdV2GlobalSum   = 0.;
    _UdUdVGlobalSum   = 0.;
    _U2dUdVGlobalSum  = 0.;
    _UdUdV2GlobalSum  = 0.;
    _Ud2UdV2GlobalSum = 0.;
    _WidomEnergyGlobalSum = 0.;
}

void LustigFormalism::ResetLocalValues()
{
	for(unsigned int i=0; i<_nWriteFreq; ++i)
		_WidomEnergyLocal[i] = 0.;

	_tsBufferIndex = 0;
}

void LustigFormalism::WriteHeader(DomainDecompBase* domainDecomp, Domain* domain)
{
#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank != 0)
    	return;
#endif

    {
		// write header
		stringstream outputstream;
		std::stringstream filenamestream;

		filenamestream << "LustigFormalism" << ".dat";
		string strFilename = filenamestream.str();

		outputstream << "      numConfigs";
		outputstream << "                  mu_res";
		outputstream << "                   _A10r";
		outputstream << "                   _A01r";
		outputstream << "                   _A20r";
		outputstream << "                   _A11r";
		outputstream << "                   _A02r";
		outputstream << "                   _A30r";
		outputstream << "                   _A21r";
		outputstream << "                   _A12r";
		outputstream << endl;

		ofstream fileout(strFilename.c_str(), ios::out);
		fileout << outputstream.str();
		fileout.close();
    }

    {
		// write header
		stringstream outputstream;
		std::stringstream filenamestream;

		filenamestream << "LustigFormalism_dUdV" << ".dat";
		string strFilename = filenamestream.str();

		outputstream << "      Config no.";
		outputstream << "                       U";
		outputstream << "                    dUdV";
		outputstream << "                  d2UdV2";
		outputstream << "             WidomEnergy";
		outputstream << endl;

		ofstream fileout(strFilename.c_str(), ios::out);
		fileout << outputstream.str();
		fileout.close();
    }

    {
		// write header
		stringstream outputstream;
		std::stringstream filenamestream;

		filenamestream << "LustigFormalism_sums" << ".dat";
		string strFilename = filenamestream.str();

		outputstream << "      numConfigs";
		outputstream << "                       U";
		outputstream << "                      U2";
		outputstream << "                      U3";
		outputstream << "                    dUdV";
		outputstream << "                  d2UdV2";
		outputstream << "                   dUdV2";
		outputstream << "                   UdUdV";
		outputstream << "                  U2dUdV";
		outputstream << "                  UdUdV2";
		outputstream << "                 Ud2UdV2";
		outputstream << "             WidomEnergy";
		outputstream << endl;

		ofstream fileout(strFilename.c_str(), ios::out);
		fileout << outputstream.str();
		fileout.close();
    }
}

void LustigFormalism::WriteData(DomainDecompBase* domainDecomp, unsigned long simstep)
{
	if( _tsBufferIndex == _nWriteFreq)
	{
		// calc global values
		this->CalcGlobalValues(domainDecomp);

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank != 0)
        return;
#endif

		this->CalcDerivatives();

		{
			// writing .dat-files
			std::stringstream outputstream;
			std::stringstream filenamestream;

			filenamestream << "LustigFormalism" << ".dat";
			string strFilename = filenamestream.str();

			// number of sampled configurations
			outputstream << std::setw(16) << _nNumConfigs;

			// data
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _mu_res;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A10r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A01r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A20r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A11r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A02r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A30r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A21r;
			outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _A12r;
			outputstream << endl;

			ofstream fileout(strFilename.c_str(), ios::app);
			fileout << outputstream.str();
			fileout.close();
		}

		{
			// writing .dat-files
			std::stringstream outputstream;
			std::stringstream filenamestream;

			filenamestream << "LustigFormalism_dUdV" << ".dat";
			string strFilename = filenamestream.str();

			// data
			unsigned long nConfigIndexBegin = _nNumConfigs - _nWriteFreq + 1;

			for(unsigned int i=0; i<_nWriteFreq; ++i)
			{
				// number of sampled configurations
				outputstream << std::setw(16) << nConfigIndexBegin + i;

				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _UGlobal[i];
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dUdVGlobal[i];
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2UdV2Global[i];
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _WidomEnergyGlobal[i];
				outputstream << endl;
			}

			ofstream fileout(strFilename.c_str(), ios::app);
			fileout << outputstream.str();
			fileout.close();
		}

	} // if(simstep % _nWriteFreq == 0)

	if( _nNumConfigs % _nWriteFreqSums != 0)
		return;

    {
		// writing .dat-files
		std::stringstream outputstream;
		std::stringstream filenamestream;

		filenamestream << "LustigFormalism_sums" << ".dat";
		string strFilename = filenamestream.str();

		// number of sampled configurations
		outputstream << std::setw(16) << _nNumConfigs;

		// data
		double dInvNumWidomTestsPerConfig = 1. / (double)(_nNumWidomTestsGlobal);

		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _UGlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _U2GlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _U3GlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dUdVGlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2UdV2GlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dUdV2GlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _UdUdVGlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _U2dUdVGlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _UdUdV2GlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _Ud2UdV2GlobalSum;
		outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _WidomEnergyGlobalSum * dInvNumWidomTestsPerConfig;
		outputstream << endl;

		ofstream fileout(strFilename.c_str(), ios::app);
		fileout << outputstream.str();
		fileout.close();
    }
}
