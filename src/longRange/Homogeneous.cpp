/*************************************************************************
 * Copyright (C) 2012 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or (at *
 * your option) any later version.                                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            * 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the Free Software           *
 * Foundation, 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.   *
 *************************************************************************/

#include <cmath>

#include "Domain.h"
#include "longRange/Homogeneous.h"
#include "Simulation.h"
//#include "LongRangeCorrection.h"

#include "utils/Logger.h"
using Log::global_log;

using namespace std;

//Homogeneous::Homogeneous(){
//  
//};

Homogeneous::Homogeneous(double cutoffRadius, double cutoffRadiusLJ, Domain* domain, Simulation _simulation) {
	_domain = domain;
	_components = *_simulation.getEnsemble()->components();
	_comp2params = _domain->getComp2Params();
	global_log->info() << "Long range correction for homogeneous systems is used " << endl;
	double UpotCorrLJ=0.;
	double VirialCorrLJ=0.;
	double MySelbstTerm=0.;
	double _globalRho = domain->getglobalRho();
	double _globalNumMolecules = domain->getglobalNumMolecules();
	double _epsilonRF = _domain->getepsilonRF();
	unsigned int numcomp=_components.size();
	unsigned long nummolecules=0;
	for(unsigned int i=0;i<numcomp;++i) {
		Component& ci=_components[i];
		nummolecules+=ci.getNumMolecules();
		unsigned int numljcentersi=ci.numLJcenters();
		unsigned int numchargesi = ci.numCharges();
		unsigned int numdipolesi=ci.numDipoles();
		unsigned int numtersoffi = ci.numTersoff();

		// effective dipoles computed from point charge distributions
		double chargeBalance[3];
		for(unsigned d = 0; d < 3; d++) chargeBalance[d] = 0;
		for(unsigned int si = 0; si < numchargesi; si++)
		{
			double tq = ci.charge(si).q();
			for(unsigned d = 0; d < 3; d++) chargeBalance[d] += tq * ci.charge(si).r()[d];
		}
		// point dipoles
		for(unsigned int si=0;si<numdipolesi;++si)
		{
			double tmy = ci.dipole(si).absMy();
			double evect = 0;
			for(unsigned d = 0; d < 3; d++) evect += ci.dipole(si).e()[d] * ci.dipole(si).e()[d];
			double norm = 1.0 / sqrt(evect);
			for(unsigned d = 0; d < 3; d++) chargeBalance[d] += tmy * ci.dipole(si).e()[d] * norm;
		}
		double my2 = 0.0;
		for(unsigned d = 0; d < 3; d++) my2 += chargeBalance[d] * chargeBalance[d];
		MySelbstTerm += my2 * ci.getNumMolecules();

		for(unsigned int j=0;j<numcomp;++j) {
			Component& cj=_components[j];
			unsigned numtersoffj = cj.numTersoff();
			// no LJ interaction between Tersoff components
			if(numtersoffi && numtersoffj) continue;
			unsigned int numljcentersj=cj.numLJcenters();
			ParaStrm& params=_comp2params(i,j);
			params.reset_read();
			// LJ centers
			for(unsigned int si=0;si<numljcentersi;++si) {
				double xi=ci.ljcenter(si).rx();
				double yi=ci.ljcenter(si).ry();
				double zi=ci.ljcenter(si).rz();
				double tau1=sqrt(xi*xi+yi*yi+zi*zi);
				for(unsigned int sj=0;sj<numljcentersj;++sj) {
					double xj=cj.ljcenter(sj).rx();
					double yj=cj.ljcenter(sj).ry();
					double zj=cj.ljcenter(sj).rz();
					double tau2=sqrt(xj*xj+yj*yj+zj*zj);
					if(tau1+tau2>=cutoffRadiusLJ){
						global_log->error() << "Error calculating cutoff corrections, rc too small" << endl;
						exit(1);
					}
					double eps24;
					params >> eps24;
					double sig2;
					params >> sig2;
					double uLJshift6;
					params >> uLJshift6;  // 0 unless TRUNCATED_SHIFTED
					if(uLJshift6 == 0.0)
					{
						double fac=double(ci.getNumMolecules())*double(cj.getNumMolecules())*eps24;
						if(tau1==0. && tau2==0.)
						{
							UpotCorrLJ+=fac*(this->_TICCu(-6,cutoffRadiusLJ,sig2)-this->_TICCu(-3,cutoffRadiusLJ,sig2));
							VirialCorrLJ+=fac*(this->_TICCv(-6,cutoffRadiusLJ,sig2)- this->_TICCv(-3,cutoffRadiusLJ,sig2));
						}
						else if(tau1!=0. && tau2!=0.)
						{
							UpotCorrLJ += fac*( this->_TISSu(-6,cutoffRadiusLJ,sig2,tau1,tau2)
									- this->_TISSu(-3,cutoffRadiusLJ,sig2,tau1,tau2) );
							VirialCorrLJ += fac*( this->_TISSv(-6,cutoffRadiusLJ,sig2,tau1,tau2)
									- this->_TISSv(-3,cutoffRadiusLJ,sig2,tau1,tau2) );
						}
						else {
							if(tau2==0.) 
								tau2=tau1;
							UpotCorrLJ+=fac*(this->_TICSu(-6,cutoffRadiusLJ,sig2,tau2)-this->_TICSu(-3,cutoffRadiusLJ,sig2,tau2));
							VirialCorrLJ+=fac*(this->_TICSv(-6,cutoffRadiusLJ,sig2,tau2)-this->_TICSv(-3,cutoffRadiusLJ,sig2,tau2));
						}
					}
				}
			}
		}
	}

	double fac=M_PI*_globalRho/(3.*_globalNumMolecules);
	UpotCorrLJ*=fac;
	VirialCorrLJ*=-fac;

	double epsRFInvrc3=2.*(_epsilonRF-1.)/((cutoffRadius*cutoffRadius*cutoffRadius)*(2.*_epsilonRF+1.));
	MySelbstTerm*=-0.5*epsRFInvrc3;

	_UpotCorr=UpotCorrLJ+MySelbstTerm;
	_VirialCorr=VirialCorrLJ+3.*MySelbstTerm;

	global_log->info() << "Far field terms: U_pot_correction  = " << _UpotCorr << " virial_correction = " << _VirialCorr << endl;
}

void Homogeneous::calculateLongRange(){
      _domain->setUpotCorr(_UpotCorr);
      _domain->setVirialCorr(_VirialCorr);  
}



double Homogeneous::_TICCu(int n,double rc,double sigma2)
{
	return -pow(rc,2*n+3) / (pow(sigma2,n)*(2*n+3));
}

double Homogeneous::_TICSu(int n,double rc,double sigma2,double tau)
{
	return -( pow(rc+tau,2*n+3) - pow(rc-tau,2*n+3) ) * rc / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3) ) 
		+  ( pow(rc+tau,2*n+4) - pow(rc-tau,2*n+4) ) / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3)*(2*n+4) );
}

double Homogeneous::_TISSu(int n,double rc,double sigma2,double tau1,double tau2)
{
	double tauMinus,tauPlus;
	tauPlus = tau1+tau2;
	tauMinus = tau1-tau2;
	return -( pow(rc+tauPlus,2*n+4) - pow(rc+tauMinus,2*n+4) - pow(rc-tauMinus,2*n+4) + pow(rc-tauPlus,2*n+4) ) * rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) ) +  ( pow(rc+tauPlus,2*n+5) - pow(rc+tauMinus,2*n+5) - pow(rc-tauMinus,2*n+5) + pow(rc-tauPlus,2*n+5) ) / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) );
}

double Homogeneous::_TICCv(int n,double rc,double sigma2)
{
	return 2*n * _TICCu(n,rc,sigma2);
}

double Homogeneous::_TICSv(int n,double rc,double sigma2,double tau)
{
	return -( pow(rc+tau,2*n+2) - pow(rc-tau,2*n+2) ) * rc*rc / ( 4*pow(sigma2,n)*tau*(n+1) ) - 3*_TICSu(n,rc,sigma2,tau);
}

double Homogeneous::_TISSv(int n,double rc,double sigma2,double tau1,double tau2){
	double tauMinus,tauPlus;
	tauPlus = tau1+tau2;
	tauMinus = tau1-tau2;
	return -(   pow(rc+tauPlus,2*n+3) - pow(rc+tauMinus,2*n+3) - pow(rc-tauMinus,2*n+3) + pow(rc-tauPlus,2*n+3) ) * rc*rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3) ) - 3*_TISSu(n,rc,sigma2,tau1,tau2);
}
