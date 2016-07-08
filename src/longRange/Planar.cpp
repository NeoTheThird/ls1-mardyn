
#include "Domain.h"
#include "longRange/Planar.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Timer.h"

#include <vector>
#include <cmath>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;

Planar::Planar(double /*cutoffT*/, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned slabs, Simulation _simulation){
	cutoff=cutoffLJ; 
	_domain = domain;
	_domainDecomposition = domainDecomposition;
	_particleContainer = particleContainer;
	_slabs = slabs;
	
	_smooth=true; // Deactivate this for transient simulations!
	global_log->info() << "Long Range Correction for planar interfaces is used" << endl;
	
	vector<Component>&  components = *_simulation.getEnsemble()->components();
	numComp=components.size();
	numLJ = (unsigned *) malloc (sizeof(unsigned)*numComp);
	numDipole = (unsigned *) malloc (sizeof(unsigned)*numComp);
	numLJSum=0;
	numDipoleSum = 0;
	numLJSum2 = (unsigned *) malloc (sizeof(unsigned)*numComp);
	numDipoleSum2 = (unsigned *) malloc (sizeof(unsigned)*numComp);
	for (unsigned i =0; i< numComp; i++){
		numLJSum2[i]=0;
		numDipoleSum2[i]=0;
	}
	for (unsigned i =0; i< numComp; i++){
		numLJ[i]=components[i].numLJcenters();
		if (components[i].numDipoles() > 0 || components[i].numCharges() != 0){
			numDipole[i]=1;
		} else{
			numDipole[i]=0;
		}
		numLJSum+=numLJ[i];
		numDipoleSum+=numDipole[i];
		for (unsigned j=i+1; j< numComp; j++){
			numLJSum2[j]+=numLJ[i];
			numDipoleSum2[j]+=numDipole[i];
		}
	}
	muSquare =  (double*)malloc(sizeof(double)*numDipoleSum);
	uLJ = (double*)malloc(sizeof(double)*_slabs*numLJSum);
	vNLJ = (double*)malloc(sizeof(double)*_slabs*numLJSum);
	vTLJ = (double*)malloc(sizeof(double)*_slabs*numLJSum);
	fLJ = (double*)malloc(sizeof(double)*_slabs*numLJSum);
	rho_g = (double*)malloc(sizeof(double)*_slabs*numLJSum);
	rho_l = (double*)malloc(sizeof(double)*_slabs*numLJSum);
	fDipole = (double*)malloc(sizeof(double)*_slabs*numDipoleSum);
	uDipole = (double*)malloc(sizeof(double)*_slabs*numDipoleSum);
	vNDipole = (double*)malloc(sizeof(double)*_slabs*numDipoleSum);
	vTDipole = (double*)malloc(sizeof(double)*_slabs*numDipoleSum);
	rhoDipole = (double*)malloc(sizeof(double)*_slabs*numDipoleSum);
	rhoDipoleL = (double*)malloc(sizeof(double)*_slabs*numDipoleSum);
	eLong =  (double*)malloc(sizeof(double)*numLJSum);
	
	unsigned counter=0;
	for (unsigned i =0; i< numComp; i++){		// Determination of the elongation of the Lennard-Jones sites
		for (unsigned j=0; j< components[i].numLJcenters(); j++){
			const LJcenter& ljcenteri = static_cast<const LJcenter&>(components[i].ljcenter(j));
			double dX[3];
			dX[0]=ljcenteri.rx();
			dX[1]=ljcenteri.ry();
			dX[2]=ljcenteri.rz();
			for (unsigned d=0; d<3; d++){
				dX[d]*=dX[d];
			}
			eLong[counter]=sqrt(dX[0]+dX[1]+dX[2]);
			counter++;
		}
	}
	for (unsigned i=0; i< _slabs*numLJSum; i++){
		rho_g[i]=0;
	}
	for (unsigned i=0; i < _slabs*numDipoleSum; i++){
		rhoDipole[i]=0;
	}
	
	unsigned dpcount=0;
	for (unsigned i=0; i<numComp; i++){
	//	for (unsigned j=0; j<numDipole[i]; j++){
		if (numDipole[i] != 0){
			double chargeBalance[3];
			for (unsigned d=0; d < 3; d++) chargeBalance[d] = 0;
			for (unsigned si=0; si < components[i].numCharges(); si++){
				double tq = components[i].charge(si).q();
				for (unsigned d=0; d < 3; d++) chargeBalance[d] += tq * components[i].charge(si).r()[d];
			}
			for (unsigned si=0; si < components[i].numDipoles(); si++){
				double tmy = components[i].dipole(si).absMy();
				for (unsigned d=0; d < 3; d++) chargeBalance[d] += tmy * components[i].dipole(si).e()[d];
			}

			// const Dipole& dipolej = static_cast<const Dipole&> (components[i].dipole(j));
			double my2 = 0;
			for (unsigned d=0; d < 3; d++)	my2 += chargeBalance[d]*chargeBalance[d];
			muSquare[dpcount]=my2;//dipolej.absMy()*dipolej.absMy();
			dpcount++;
			cout << dpcount << "\tmu: " << muSquare[dpcount-1] << endl;
		}
	} 	
	
	ymax=_domain->getGlobalLength(1);
	boxlength[0]=_domain->getGlobalLength(0);
	boxlength[2]=_domain->getGlobalLength(2);
	cutoff_slabs=cutoff*_slabs/ymax; // Number of slabs to cutoff
	delta=ymax/_slabs;	
	frequency=10;
	V=ymax*_domain->getGlobalLength(0)*_domain->getGlobalLength(2); // Volume
	
	sint=_slabs;
	simstep = 0;
	
	temp=_domain->getTargetTemperature(0);
}


void Planar::calculateLongRange(){
	Molecule* tempMol;
	if (_smooth){
		for( tempMol = _particleContainer->begin(); tempMol != _particleContainer->end(); tempMol = _particleContainer->next()){
			unsigned cid=tempMol->componentid();
			for (unsigned i=0; i<numLJ[cid]; i++){
				int loc=(tempMol->r(1)+tempMol->ljcenter_d(i)[1])/delta;
				if (loc < 0){
					loc=loc+_slabs;
				}
				else if (loc > sint-1){
					loc=loc-_slabs;
				}
				rho_g[loc+_slabs*(i+numLJSum2[cid])]+=_slabs/V;
			}
			if (numDipole[cid] != 0){
				int loc=tempMol->r(1)/delta;
				rhoDipole[loc+_slabs*(numDipoleSum2[cid])]+=_slabs/V;
			}
		}
	} 
	if (simstep % frequency == 0){	// The Density Profile is only calculated once in 10 simulation steps

		for (unsigned i=0; i<_slabs*numLJSum; i++){
			rho_l[i]=0;
			uLJ[i]=0;
			vNLJ[i]=0;
			vTLJ[i]=0;
			fLJ[i]=0;
		}
		
		for (unsigned i=0; i<_slabs*numDipoleSum; i++){
			fDipole[i]=0;
			uDipole[i]=0;
			vNDipole[i]=0;
			vTDipole[i]=0;
			rhoDipoleL[i]=0;
		}
		// Calculation of the density profile for s slabs
		if (!_smooth){
			for( tempMol = _particleContainer->begin(); tempMol != _particleContainer->end(); tempMol = _particleContainer->next()){
				unsigned cid=tempMol->componentid();
				for (unsigned i=0; i<numLJ[cid]; i++){
					int loc=(tempMol->r(1)+tempMol->ljcenter_d(i)[1])/delta;
					if (loc < 0){
						loc=loc+_slabs;
					}
					else if (loc > sint-1){
						loc=loc-_slabs;
					}
					rho_l[loc+_slabs*(i+numLJSum2[cid])]+=_slabs/V;
				}
				if (numDipole[cid] != 0){
					int loc=tempMol->r(1)/delta;
					rhoDipoleL[loc+_slabs*numDipoleSum2[cid]]+=_slabs/V;
				}
			}
		}
		else{
			for (unsigned i=0; i<_slabs*numLJSum; i++){
				rho_l[i]=rho_g[i]/(simstep+1);
			}
			for (unsigned i=0; i<_slabs*numDipoleSum; i++){
				rhoDipoleL[i]=rhoDipole[i]/(simstep+1);
			}
		}
		
		// Distribution of the Density Profile to every node
		_domainDecomposition->collCommInit(_slabs*(numLJSum+numDipoleSum));
		for (unsigned i=0; i < _slabs*numLJSum; i++) {
			_domainDecomposition->collCommAppendDouble(rho_l[i]);
		}
		for (unsigned i=0; i< _slabs*numDipoleSum; i++){
			_domainDecomposition->collCommAppendDouble(rhoDipoleL[i]);
		}
	
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i=0; i < _slabs*numLJSum; i++) {
			rho_l[i] = _domainDecomposition->collCommGetDouble();
		}
		for (unsigned i=0; i< _slabs*numDipoleSum; i++){
			rhoDipoleL[i] = _domainDecomposition->collCommGetDouble();
		}
		
		_domainDecomposition->collCommFinalize();

		for (unsigned ci = 0; ci < numComp; ++ci){		
			for (unsigned cj = 0; cj < numComp; ++cj){	
				ParaStrm& params = _domain->getComp2Params()(ci,cj);
				params.reset_read();
				for (unsigned si = 0; si < numLJ[ci]; ++si) { // Long Range Correction for Lennard-Jones sites
					for (unsigned sj = 0; sj < numLJ[cj]; ++sj) {
						double eps24;
						double sig2;
						double shift6;
						double eps;
						params >> eps24;
						params >> sig2;
						params >> shift6;
						sig2=sqrt(sig2);
						eps=eps24/24;
						if (eLong[numLJSum2[ci]+si] ==0 && eLong[numLJSum2[cj]+sj] == 0){
							centerCenter(sig2,eps,ci,cj,si,sj);
						}
						else if (eLong[numLJSum2[ci]+si] ==0 || eLong[numLJSum2[cj]+sj] == 0){
							centerSite(sig2,eps,ci,cj,si,sj);
						}
						else{
							siteSite(sig2,eps,ci,cj,si,sj);
						}
					}
				}
		
				for (unsigned si=0; si< numDipole[ci]; si++){	//Long Range Correction for Dipoles
					for (unsigned sj=0; sj< numDipole[cj]; sj++){
						 dipoleDipole(ci,cj,si,sj);	
					}
				}
			}
	 	} 
	     
		// Distribution of the Force, Energy and Virial to every Node
		_domainDecomposition->collCommInit(_slabs*(4*numLJSum+4*numDipoleSum));
		for (unsigned i=0; i<_slabs*numLJSum; i++){
			_domainDecomposition->collCommAppendDouble(uLJ[i]);
			_domainDecomposition->collCommAppendDouble(vNLJ[i]);
			_domainDecomposition->collCommAppendDouble(vTLJ[i]);
			_domainDecomposition->collCommAppendDouble(fLJ[i]);
		}
		for (unsigned i=0; i<_slabs*numDipoleSum; i++){
			_domainDecomposition->collCommAppendDouble(uDipole[i]);
			_domainDecomposition->collCommAppendDouble(fDipole[i]);
			_domainDecomposition->collCommAppendDouble(vNDipole[i]);
			_domainDecomposition->collCommAppendDouble(vTDipole[i]);
		}	
		_domainDecomposition->collCommAllreduceSum();
		for (unsigned i=0; i<_slabs*numLJSum; i++){
			uLJ[i]=_domainDecomposition->collCommGetDouble();
			vNLJ[i]=_domainDecomposition->collCommGetDouble();
			vTLJ[i]=_domainDecomposition->collCommGetDouble();
			fLJ[i]=_domainDecomposition->collCommGetDouble();
		}
		for (unsigned i=0; i<_slabs*numDipoleSum; i++){
		      uDipole[i]=_domainDecomposition->collCommGetDouble();
		      fDipole[i]=_domainDecomposition->collCommGetDouble();
		      vNDipole[i]=_domainDecomposition->collCommGetDouble();
		      vTDipole[i]=_domainDecomposition->collCommGetDouble();
		}
		_domainDecomposition->collCommFinalize();
	}

	// Adding the Force to the Molecules; this is done in every timestep
	double Fa[3]={0};
	double Via[3]={0};
	double Upot_c=0;
	double Virial_c=0;	// Correction used for the Pressure Calculation
	for( tempMol = _particleContainer->begin(); tempMol != _particleContainer->end(); tempMol = _particleContainer->next()){
	        unsigned cid=tempMol->componentid();
		for (unsigned i=0; i<numLJ[cid]; i++){
			int loc=(tempMol->r(1)+tempMol->ljcenter_d(i)[1])/delta;
			if (loc < 0){
				loc=loc+_slabs;
			}
			else if (loc > sint-1){
				loc=loc-_slabs;
			}
			Fa[1]=fLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
			Upot_c+=uLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
			Virial_c+=2*vTLJ[loc+i*_slabs+_slabs*numLJSum2[cid]]+vNLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
			Via[0]=vTLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
			Via[1]=vNLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
			Via[2]=vTLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
			tempMol->Fljcenteradd(i,Fa);
			Virial_c=Via[1];
			tempMol->Viadd(Via);
//			tempMol->Uadd(uLJ[loc+i*s+_slabs*numLJSum2[cid]]);	// Storing potential energy onto the molecules is currently not implemented!
		}
		if (numDipole[cid] != 0){
			int loc = tempMol->r(1)/delta;
			Fa[1]=fDipole[loc+_slabs*numDipoleSum2[cid]];
			Upot_c+=uDipole[loc+_slabs*numDipoleSum2[cid]];
			Virial_c+=2*vTDipole[loc+_slabs*numDipoleSum2[cid]]+vNDipole[loc+_slabs*numDipoleSum2[cid]];
			Via[0]=vTDipole[loc+_slabs*numDipoleSum2[cid]];
			Via[1]=vNDipole[loc+_slabs*numDipoleSum2[cid]];
			Via[2]=vTDipole[loc+_slabs*numDipoleSum2[cid]];
			tempMol->Fadd(Fa); // Force is stored on the center of mass of the molecule!
			tempMol->Viadd(Via);
//			tempMol->Uadd(uDipole[loc+i*_slabs+_slabs*numDipoleSum2[cid]]);	// Storing potential energy onto the molecules is currently not implemented!
		//	}
		}				
	}

	// Summation of the correction terms
	_domainDecomposition->collCommInit(2);
	_domainDecomposition->collCommAppendDouble(Upot_c);
	_domainDecomposition->collCommAppendDouble(Virial_c);
	_domainDecomposition->collCommAllreduceSum();
	Upot_c = _domainDecomposition->collCommGetDouble();
	Virial_c = _domainDecomposition->collCommGetDouble();
	_domainDecomposition->collCommFinalize();

	// Setting the Energy and Virial correction
	_domain->setUpotCorr(Upot_c);
	_domain->setVirialCorr(Virial_c);	
	
	simstep++;
}


void Planar::centerCenter(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	double rc=sig/cutoff;
	double rc2=rc*rc;
	double rc6=rc2*rc2*rc2;
	double rc12=rc6*rc6;
	double r,r2,r6,r12;
	double rhoI,rhoJ;
	double termU = 4*3.1416*delta*eps*sig*sig;
	double termF = 8*3.1416*delta*eps*sig;
	double termVN = 4*3.1416*delta*eps*sig*sig;
	double termVT = 2*3.1416*delta*eps*sig*sig;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()){
		rhoI=rho_l[i+si*_slabs+_slabs*numLJSum2[ci]];
		vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVT*rhoI*(6*rc12/5-3*rc6/2)/rc2;
		uLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termU*rhoI*(rc12/5-rc6/2)/rc2;
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			r=sig/((j-i)*delta);
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j> i+cutoff_slabs){
				r2=r*r;
				r6=r2*r2*r2;
				r12=r6*r6;
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=termVT*rhoJ*(r12/5-r6/2)/r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=termVT*rhoI*(r12/5-r6/2)/r2;
			}
			else{
				r2=rc2;
				r6=rc6;
				r12=rc12;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=termVT*rhoI*(r12/5*(6/r2-5/(r*r))-r6/2*(3/r2-2/(r*r)));	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=termVT*rhoJ*(r12/5*(6/r2-5/(r*r))-r6/2*(3/r2-2/(r*r)));
			}
			uLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termU*rhoJ*(r12/5-r6/2)/r2;
			uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termU*rhoI*(r12/5-r6/2)/r2;
			vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVN*rhoJ*(r12-r6)/(r*r);
			vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVN*rhoI*(r12-r6)/(r*r);
			fLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += -termF*rhoJ*(r12-r6)/r;
			fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termF*rhoI*(r12-r6)/r;
		}
		// Calculation of the Periodic boundary 
		for (unsigned j=_slabs/2+i; j<_slabs; j++){
			r=sig/((_slabs-j+i)*delta);
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j <_slabs-cutoff_slabs+i){
				r2=r*r;
				r6=r2*r2*r2;
				r12=r6*r6;
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVT*rhoJ*(r12/5-r6/2)/r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVT*rhoI*(r12/5-r6/2)/r2;
			}
			else{
				r2=rc2;
				r6=rc6;
				r12=rc12;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVT*rhoI*(r12/5*(6/r2-5/(r*r))-r6/2*(3/r2-2/(r*r)));		
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVT*rhoJ*(r12/5*(6/r2-5/(r*r))-r6/2*(3/r2-2/(r*r)));
			}
			uLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termU*rhoJ*(r12/5-r6/2)/r2;
			uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termU*rhoI*(r12/5-r6/2)/r2;
			vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVN*rhoJ*(r12-r6)/(r*r);
			vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVN*rhoI*(r12-r6)/(r*r);
			fLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termF*rhoJ*(r12-r6)/r;
			fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += -termF*rhoI*(r12-r6)/r;
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()){
		rhoI=rho_l[i+si*_slabs+_slabs*numLJSum2[ci]];
		vTLJ[i+si*_slabs+_slabs*numLJSum2[cj]] += termVT*rhoI*(6*rc12/5-3*rc6/2)/rc2;
		uLJ[i+si*_slabs+_slabs*numLJSum2[cj]] += termU*rhoI*(rc12/5-rc6/2)/rc2;
		for (unsigned j=i+1; j<_slabs; j++){
			r=sig/((j-i)*delta);
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j> i+cutoff_slabs){
				r2=r*r;
				r6=r2*r2*r2;
				r12=r6*r6;
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVT*rhoJ*(r12/5-r6/2)/r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVT*rhoI*(r12/5-r6/2)/r2;
			}
			else{
				r2=rc2;
				r6=rc6;
				r12=rc12;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVT*rhoI*(r12/5*(6/r2-5/(r*r))-r6/2*(3/r2-2/(r*r)));
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVT*rhoJ*(r12/5*(6/r2-5/(r*r))-r6/2*(3/r2-2/(r*r)));
			}
			uLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termU*rhoJ*(r12/5-r6/2)/r2;
			uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termU*rhoI*(r12/5-r6/2)/r2;
			vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += termVN*rhoJ*(r12-r6)/(r*r);
			vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termVN*rhoI*(r12-r6)/(r*r);
			fLJ[i+si*_slabs+_slabs*numLJSum2[ci]] += -termF*rhoJ*(r12-r6)/r;
			fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]] += termF*rhoI*(r12-r6)/r;
		}
	}
  
}

void Planar::centerSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	double sig2=sig*sig;
	double sig3=sig2*sig;
	double t = eLong[numLJSum2[ci]+si] + eLong[numLJSum2[cj]+sj]; // one of them is equal to zero.
	double rcPt=sig/(cutoff+t);
	double rcPt3=rcPt*rcPt*rcPt;
	double rcPt4=rcPt3*rcPt;
	double rcPt9=rcPt3*rcPt3*rcPt3;
	double rcPt10=rcPt9*rcPt;
	double rcMt=sig/(cutoff-t);
	double rcMt3=rcMt*rcMt*rcMt;
	double rcMt4=rcMt3*rcMt;
	double rcMt9=rcMt3*rcMt3*rcMt3;
	double rcMt10=rcMt9*rcMt;
	double rc=cutoff;
	double rc2=rc*rc;
	double termURC=-2*3.1416*eps*delta*sig3/(3*t)*((rcPt9-rcMt9)/15-(rcPt3-rcMt3)/2);
	double termFRC=-2*3.1416*eps*delta*sig2/(t*rc)*((rcPt10-rcMt10)/5-(rcPt4-rcMt4)/2);
	double termVNRC=termFRC/2;
	double termVTRC1=-3.1416*eps*delta*sig2/(2*t*rc)*((rcPt10-rcMt10)/5-(rcPt4-rcMt4)/2);
	double termVTRC2=termURC/2;
	double r,r2;
	double rhoI,rhoJ;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()){
		rhoI=rho_l[i+si*_slabs+_slabs*numLJSum2[ci]];
		rhoJ=rho_l[i+sj*_slabs+_slabs*numLJSum2[cj]];
		vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j> i+cutoff_slabs){
				double rPt=sig/(r+t);
				double rPt3=rPt*rPt*rPt;
				double rPt4=rPt3*rPt;
				double rPt9=rPt3*rPt3*rPt3;
				double rPt10=rPt9*rPt;
				double rMt=sig/(r-t);
				double rMt3=rMt*rMt*rMt;
				double rMt4=rMt3*rMt;
				double rMt9=rMt3*rMt3*rMt3;
				double rMt10=rMt9*rMt;
				double termU=-2*3.1416*eps*delta*sig3/(3*t)*((rPt9-rMt9)/15-(rPt3-rMt3)/2);
				double termF=-2*3.1416*eps*delta*sig2/(t*r)*((rPt10-rMt10)/5-(rPt4-rMt4)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termU;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termU;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVN*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVN*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVT2;	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVT2;
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termF*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termF*r;
			}
			else{
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termURC;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVNRC*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVNRC*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termFRC*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termFRC*r;
			}
		}
		// Calculation of the Periodic boundary 
		for (unsigned j=_slabs/2+i; j<_slabs; j++){
			r=(_slabs-j+i)*delta;
			r2=r*r;
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j <_slabs-cutoff_slabs+i){
				double rPt=sig/(r+t);
				double rPt3=rPt*rPt*rPt;
				double rPt4=rPt3*rPt;
				double rPt9=rPt3*rPt3*rPt3;
				double rPt10=rPt9*rPt;
				double rMt=sig/(r-t);
				double rMt3=rMt*rMt*rMt;
				double rMt4=rMt3*rMt;
				double rMt9=rMt3*rMt3*rMt3;
				double rMt10=rMt9*rMt;
				double termU=-2*3.1416*eps*delta*sig3/(3*t)*((rPt9-rMt9)/15-(rPt3-rMt3)/2);
				double termF=-2*3.1416*eps*delta*sig2/(t*r)*((rPt10-rMt10)/5-(rPt4-rMt4)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termU;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termU;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVN*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVN*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVT2;	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVT2;
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termF*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=-rhoI*termF*r;
			}
			else{
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termURC;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVNRC*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVNRC*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termFRC*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=-rhoI*termFRC*r;
			}
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()){
		rhoI=rho_l[i+si*_slabs+_slabs*numLJSum2[ci]];
		rhoJ=rho_l[i+sj*_slabs+_slabs*numLJSum2[cj]];
		vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
		for (unsigned j=i+1; j<_slabs; j++){
			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j> i+cutoff_slabs){
				double rPt=sig/(r+t);
				double rPt3=rPt*rPt*rPt;
				double rPt4=rPt3*rPt;
				double rPt9=rPt3*rPt3*rPt3;
				double rPt10=rPt9*rPt;
				double rMt=sig/(r-t);
				double rMt3=rMt*rMt*rMt;
				double rMt4=rMt3*rMt;
				double rMt9=rMt3*rMt3*rMt3;
				double rMt10=rMt9*rMt;
				double termU=-2*3.1416*eps*delta*sig3/(3*t)*((rPt9-rMt9)/15-(rPt3-rMt3)/2);
				double termF=-2*3.1416*eps*delta*sig2/(t*r)*((rPt10-rMt10)/5-(rPt4-rMt4)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termU;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termU;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVN*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVN*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVT2;	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVT2;
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termF*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termF*r;
			}
			else{
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termURC;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVNRC*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVNRC*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termFRC*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termFRC*r;
			}
		}
	}
  
}

void Planar::siteSite(double sig, double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj){
	double sig2=sig*sig;
	double sig3=sig2*sig;
	double sig4=sig2*sig2;
	double t1 = eLong[numLJSum2[ci]+si];
	double t2 = eLong[numLJSum2[cj]+sj];
	double tP = t1 + t2; // tau+ 
	double tM = t1 - t2; // tau-
	double rcPtP=sig/(cutoff+tP);
	double rcPtP2=rcPtP*rcPtP;
	double rcPtP3=rcPtP2*rcPtP;
	double rcPtP8=rcPtP2*rcPtP3*rcPtP3;
	double rcPtP9=rcPtP8*rcPtP;
	double rcPtM=sig/(cutoff+tM);
	double rcPtM2=rcPtM*rcPtM;
	double rcPtM3=rcPtM2*rcPtM;
	double rcPtM8=rcPtM2*rcPtM3*rcPtM3;
	double rcPtM9=rcPtM8*rcPtM;
	double rcMtP=sig/(cutoff-tP);
	double rcMtP2=rcMtP*rcMtP;
	double rcMtP3=rcMtP2*rcMtP;
	double rcMtP8=rcMtP2*rcMtP3*rcMtP3;
	double rcMtP9=rcMtP8*rcMtP;
	double rcMtM=sig/(cutoff-tM);
	double rcMtM2=rcMtM*rcMtM;
	double rcMtM3=rcMtM2*rcMtM;
	double rcMtM8=rcMtM2*rcMtM3*rcMtM3;
	double rcMtM9=rcMtM8*rcMtM;
	double rc=cutoff;
	double rc2=rc*rc;
	double termURC=3.1416*eps*delta*sig4/(12*t1*t2)*((rcPtP8-rcPtM8-rcMtM8+rcMtP8)/30-(rcPtP2-rcPtM2-rcMtM2+rcMtP2));
	double termFRC=3.1416*eps*delta*sig3/(3*t1*t2*rc)*((rcPtP9-rcPtM9-rcMtM9+rcMtP9)/15-(rcPtP3-rcPtM3-rcMtM3+rcMtP3)/2);
	double termVNRC=termFRC/2;
	double termVTRC1=3.1416*eps*delta*sig3/(4*t1*t2*rc)*((rcPtP9-rcPtM9-rcMtM9+rcMtP9)/45-(rcPtP3-rcPtM3-rcMtM3+rcMtP3)/6);
	double termVTRC2=termURC/2;
	double r,r2;
	double rhoI,rhoJ;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()){
		rhoI=rho_l[i+si*_slabs+_slabs*numLJSum2[ci]];
		rhoJ=rho_l[i+sj*_slabs+_slabs*numLJSum2[cj]];
		vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j> i+cutoff_slabs){
				double rPtP=sig/(r+tP);
				double rPtP2=rPtP*rPtP;
				double rPtP3=rPtP2*rPtP;
				double rPtP8=rPtP2*rPtP3*rPtP3;
				double rPtP9=rPtP8*rPtP;
				double rPtM=sig/(r+tM);
				double rPtM2=rPtM*rPtM;
				double rPtM3=rPtM2*rPtM;
				double rPtM8=rPtM2*rPtM3*rPtM3;
				double rPtM9=rPtM8*rPtM;
				double rMtP=sig/(r-tP);
				double rMtP2=rMtP*rMtP;
				double rMtP3=rMtP2*rMtP;
				double rMtP8=rMtP2*rMtP3*rMtP3;
				double rMtP9=rMtP8*rMtP;
				double rMtM=sig/(r-tM);
				double rMtM2=rMtM*rMtM;
				double rMtM3=rMtM2*rMtM;
				double rMtM8=rMtM2*rMtM3*rMtM3;
				double rMtM9=rMtM8*rMtM;
				double termU=3.1416*eps*delta*sig4/(12*t1*t2)*((rPtP8-rPtM8-rMtM8+rMtP8)/30-(rPtP2-rPtM2-rMtM2+rMtP2));
				double termF=3.1416*eps*delta*sig3/(3*t1*t2*r)*((rPtP9-rPtM9-rMtM9+rMtP9)/15-(rPtP3-rPtM3-rMtM3+rMtP3)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termU;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termU;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVN*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVN*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVT2;	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVT2;
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termF*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termF*r;
			}
			else{
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termURC;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVNRC*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVNRC*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termFRC*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termFRC*r;
			}
		}
		// Calculation of the Periodic boundary 
		for (unsigned j=_slabs/2+i; j<_slabs; j++){
			r=(_slabs-j+i)*delta;
			r2=r*r;
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j <_slabs-cutoff_slabs+i){
				double rPtP=sig/(r+tP);
				double rPtP2=rPtP*rPtP;
				double rPtP3=rPtP2*rPtP;
				double rPtP8=rPtP2*rPtP3*rPtP3;
				double rPtP9=rPtP8*rPtP;
				double rPtM=sig/(r+tM);
				double rPtM2=rPtM*rPtM;
				double rPtM3=rPtM2*rPtM;
				double rPtM8=rPtM2*rPtM3*rPtM3;
				double rPtM9=rPtM8*rPtM;
				double rMtP=sig/(r-tP);
				double rMtP2=rMtP*rMtP;
				double rMtP3=rMtP2*rMtP;
				double rMtP8=rMtP2*rMtP3*rMtP3;
				double rMtP9=rMtP8*rMtP;
				double rMtM=sig/(r-tM);
				double rMtM2=rMtM*rMtM;
				double rMtM3=rMtM2*rMtM;
				double rMtM8=rMtM2*rMtM3*rMtM3;
				double rMtM9=rMtM8*rMtM;
				double termU=3.1416*eps*delta*sig4/(12*t1*t2)*((rPtP8-rPtM8-rMtM8+rMtP8)/30-(rPtP2-rPtM2-rMtM2+rMtP2));
				double termF=3.1416*eps*delta*sig3/(3*t1*t2*r)*((rPtP9-rPtM9-rMtM9+rMtP9)/15-(rPtP3-rPtM3-rMtM3+rMtP3)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termU;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termU;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVN*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVN*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVT2;	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVT2;
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termF*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=-rhoI*termF*r;
			}
			else{
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termURC;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVNRC*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVNRC*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termFRC*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=-rhoI*termFRC*r;
			}
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()){
		rhoI=rho_l[i+si*_slabs+_slabs*numLJSum2[ci]];
		rhoJ=rho_l[i+sj*_slabs+_slabs*numLJSum2[cj]];
		vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*rc2+termVTRC2);
		uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
		for (unsigned j=i+1; j<_slabs; j++){
			r=(j-i)*delta;
			r2=r*r;
			rhoJ=rho_l[j+sj*_slabs+_slabs*numLJSum2[cj]];
			if (j> i+cutoff_slabs){
				double rPtP=sig/(r+tP);
				double rPtP2=rPtP*rPtP;
				double rPtP3=rPtP2*rPtP;
				double rPtP8=rPtP2*rPtP3*rPtP3;
				double rPtP9=rPtP8*rPtP;
				double rPtM=sig/(r+tM);
				double rPtM2=rPtM*rPtM;
				double rPtM3=rPtM2*rPtM;
				double rPtM8=rPtM2*rPtM3*rPtM3;
				double rPtM9=rPtM8*rPtM;
				double rMtP=sig/(r-tP);
				double rMtP2=rMtP*rMtP;
				double rMtP3=rMtP2*rMtP;
				double rMtP8=rMtP2*rMtP3*rMtP3;
				double rMtP9=rMtP8*rMtP;
				double rMtM=sig/(r-tM);
				double rMtM2=rMtM*rMtM;
				double rMtM3=rMtM2*rMtM;
				double rMtM8=rMtM2*rMtM3*rMtM3;
				double rMtM9=rMtM8*rMtM;
				double termU=3.1416*eps*delta*sig4/(12*t1*t2)*((rPtP8-rPtM8-rMtM8+rMtP8)/30-(rPtP2-rPtM2-rMtM2+rMtP2));
				double termF=3.1416*eps*delta*sig3/(3*t1*t2*r)*((rPtP9-rPtM9-rMtM9+rMtP9)/15-(rPtP3-rPtM3-rMtM3+rMtP3)/2);
				double termVN=termF/2;
				double termVT2=termU/2;
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termU;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termU;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVN*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVN*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVT2;	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVT2;
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termF*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termF*r;
			}
			else{
				uLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termURC;
				uLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termURC;
				vNLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*termVNRC*r2;
				vNLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termVNRC*r2;
				vTLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*(termVTRC1*(rc2-r2)+termVTRC2);	
				vTLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=rhoJ*(termVTRC1*(rc2-r2)+termVTRC2);
				fLJ[i+si*_slabs+_slabs*numLJSum2[ci]]+=-rhoJ*termFRC*r;
				fLJ[j+sj*_slabs+_slabs*numLJSum2[cj]]+=rhoI*termFRC*r;
			}
		}
	}
  
}

void Planar::dipoleDipole(unsigned ci,unsigned cj,unsigned si,unsigned sj){
	double rc=cutoff;
	double rc2=rc*rc;
	double rc4=rc2*rc2;
	double rc6=rc4*rc2;
	double termU = 3.1416/4*muSquare[ci]*muSquare[cj]*delta / (3*temp);
	double termF = 3.1416 * muSquare[ci]*muSquare[cj]*delta / (3*temp);
	double termVN= 3.1416/2*muSquare[ci]*muSquare[cj]*delta / (3*temp);
	double termVT= termU;
	for (unsigned i=_domainDecomposition->getRank(); i<_slabs/2; i+=_domainDecomposition->getNumProcs()){
		double rhoI = rhoDipoleL[i+si*_slabs+_slabs*numDipoleSum2[ci]];
		for (unsigned j=i+1; j<i+_slabs/2; j++){
			double rhoJ = rhoDipoleL[j+sj*_slabs+_slabs*numDipoleSum2[cj]];
			double r=(j-i)*delta;
			double r2,r4,r6;
			if (j> i+cutoff_slabs){
			r2=r*r;
			r4=r2*r2;
			r6=r4*r2;
			}
			else{
			  r2=rc2;
			  r4=rc4;
			  r6=rc6;
			}
			fDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]] += termF*rhoJ/r6 * r;
			fDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]] -= termF*rhoI/r6 * r;
			uDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]] -= termU*rhoJ/r4;
			uDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]] -= termU*rhoI/r4;
			vNDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]]-= termVN*rhoJ/r6 *r*r;
			vNDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]]-= termVN*rhoI/r6 *r*r;
			vTDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]]-= termVT*rhoJ/r6 *(1.5*r2 - r*r);
			vTDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]]-= termVT*rhoI/r6 *(1.5*r2 - r*r);
			
		}
		// Calculation of the Periodic boundary 
		for (unsigned j=_slabs/2+i; j<_slabs; j++){
			double rhoJ = rhoDipoleL[j+sj*_slabs+_slabs*numDipoleSum2[cj]];
			double r=(_slabs-j+i)*delta;
			double r2,r4,r6;
			if (j <_slabs-cutoff_slabs+i){
			r2=r*r;
			r4=r2*r2;
			r6=r4*r2;
			}
			else{
			  r2=rc2;
			  r4=rc4;
			  r6=rc6;
			}
			fDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]] -= termF*rhoJ/r6 * r;
			fDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]] += termF*rhoI/r6 * r;
			uDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]] -= termU*rhoJ/r4;
			uDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]] -= termU*rhoI/r4;
			vNDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]]-= termVN*rhoJ/r6 *r*r;
			vNDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]]-= termVN*rhoI/r6 *r*r;
			vTDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]]-= termVT*rhoJ/r6 *(1.5*r2 - r*r);
			vTDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]]-= termVT*rhoI/r6 *(1.5*r2 - r*r);
		}
	}

	// Calculation of the Forces on the slabs of the right hand side
	for (unsigned i=_slabs/2+_domainDecomposition->getRank(); i<_slabs; i+=_domainDecomposition->getNumProcs()){
		double rhoI = rhoDipoleL[i+si*_slabs+_slabs*numDipoleSum2[ci]];
		for (unsigned j=i+1; j<_slabs; j++){
			double rhoJ = rhoDipoleL[j+sj*_slabs+_slabs*numDipoleSum2[cj]];
			double r=(j-i)*delta;
			double r2,r4,r6;
			if (j> i+cutoff_slabs){
			r2=r*r;
			r4=r2*r2;
			r6=r4*r2;
			}
			else{
			  r2=rc2;
			  r4=rc4;
			  r6=rc6;
			}
			fDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]] += termF*rhoJ/r6 * r;
			fDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]] -= termF*rhoI/r6 * r;
			uDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]] -= termU*rhoJ/r4;
			uDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]] -= termU*rhoI/r4;
			vNDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]]-= termVN*rhoJ/r6 *r*r;
			vNDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]]-= termVN*rhoI/r6 *r*r;
			vTDipole[i+si*_slabs+_slabs*numDipoleSum2[ci]]-= termVT*rhoJ/r6 *(1.5*r2 - r*r);
			vTDipole[j+sj*_slabs+_slabs*numDipoleSum2[cj]]-= termVT*rhoI/r6 *(1.5*r2 - r*r);
		}
	}  
}

double Planar::lrcLJ(Molecule* mol){
	double potentialEnergy = 0.;
	unsigned cid=mol->componentid();
	for (unsigned i=0; i<numLJ[cid]; i++){
		int loc=(mol->r(1)+mol->ljcenter_d(i)[1])/delta;
		if (loc < 0){
			loc=loc+_slabs;
		}
		else if (loc > sint-1){
			loc=loc-_slabs;
		}
		potentialEnergy += uLJ[loc+i*_slabs+_slabs*numLJSum2[cid]];
	}
	for (unsigned i=0;i<numDipole[cid]; i++){
		int loc=(mol->r(1)+mol->dipole_d(i)[1])/delta;
		if (loc < 0){
			loc=loc+_slabs;
		}
		else if (loc > sint-1){
			loc=loc-_slabs;
		}
		potentialEnergy += uDipole[loc+i*_slabs+_slabs*numDipoleSum2[cid]];
	}
	return potentialEnergy;
}

void Planar::directDensityProfile(){
	_smooth = false;
}
