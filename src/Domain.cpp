/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <iostream>
#include <string>
#include <cmath>

#include "Domain.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"
#include "ensemble/PressureGradient.h"
#include "CutoffCorrections.h"

#include "utils/Logger.h"
using Log::global_log;

using namespace std;


Domain::Domain(int rank, PressureGradient* pg){
	_localRank = rank;
	_localUpot = 0;
	_localVirial = 0;   
	_globalUpot = 0;
	_globalVirial = 0; 
	_globalRho = 0;

	this->_universalPG = pg;

	this->_componentToThermostatIdMap = map<int, int>();
	this->_localThermostatN = map<int, unsigned long>();
	this->_localThermostatN[-1] = 0;
	this->_localThermostatN[0] = 0;
	this->_universalThermostatN = map<int, unsigned long>();
	this->_universalThermostatN[-1] = 0;
	this->_universalThermostatN[0] = 0;
	this->_localRotationalDOF = map<int, unsigned long>();
	this->_localRotationalDOF[-1] = 0;
	this->_localRotationalDOF[0] = 0;
	this->_universalRotationalDOF = map<int, unsigned long>();
	this->_universalRotationalDOF[-1] = 0;
	this->_universalRotationalDOF[0] = 0;
	this->_globalLength[0] = 0;
	this->_globalLength[1] = 0;
	this->_globalLength[2] = 0;
	this->_universalBTrans = map<int, double>();
	this->_universalBTrans[0] = 1.0;
	this->_universalBRot = map<int, double>();
	this->_universalBRot[0] = 1.0;
	this->_universalTargetTemperature = map<int, double>();
	this->_universalTargetTemperature[0] = 1.0;
	this->_globalTemperatureMap = map<int, double>();
	this->_globalTemperatureMap[0] = 1.0;
	this->_local2KETrans[0] = 0.0;
	this->_local2KERot[0] = 0.0; 

	_currentTime = 0.0;

	this->_universalNVE = false;
	this->_globalUSteps = 0;
	this->_globalSigmaU = 0.0;
	this->_globalSigmaUU = 0.0;
#ifdef COMPLEX_POTENTIAL_SET
	this->_universalConstantAccelerationTimesteps = 30;
	if(!rank)
		for(int d = 0; d < 3; d++)
			this->_globalVelocitySum[d] = map<unsigned, long double>();
#endif
	this->_componentwiseThermostat = false;
#ifdef COMPLEX_POTENTIAL_SET
	this->_universalUndirectedThermostat = map<int, bool>();
	for(int d = 0; d < 3; d++)
	{
		this->_universalThermostatDirectedVelocity[d] = map<int, double>();
		this->_localThermostatDirectedVelocity[d] = map<int, double>();
	}
#endif
	this->_universalSelectiveThermostatCounter = 0;
	this->_universalSelectiveThermostatWarning = 0;
	this->_universalSelectiveThermostatError = 0;
	// default
	this->_universalCylindricalGeometry = false;
}

void Domain::setLocalUpot(double Upot) {_localUpot = Upot;}

double Domain::getLocalUpot() const {return _localUpot; }

void Domain::setLocalVirial(double Virial) {_localVirial = Virial;}

double Domain::getLocalVirial() const {return _localVirial; }

/* methods accessing thermostat info */
double Domain::getGlobalBetaTrans() { return _universalBTrans[0]; }
double Domain::getGlobalBetaTrans(int thermostat) { return _universalBTrans[thermostat]; }
double Domain::getGlobalBetaRot() { return _universalBRot[0]; }
double Domain::getGlobalBetaRot(int thermostat) { return _universalBRot[thermostat]; }

void Domain::setLocalSummv2(double summv2, int thermostat)
{
#ifndef NDEBUG
	global_log->debug() << "* local thermostat " << thermostat << ":  mvv = " << summv2 << endl;
#endif
	this->_local2KETrans[thermostat] = summv2;
}

void Domain::setLocalSumIw2(double sumIw2, int thermostat)
{
	_local2KERot[thermostat] = sumIw2;
} 

double Domain::getGlobalPressure()
{
	double globalTemperature = _globalTemperatureMap[0];
	return globalTemperature * _globalRho + _globalRho * getAverageGlobalVirial()/3.;
}

double Domain::getAverageGlobalVirial() const { return _globalVirial/_globalNumMolecules; }

double Domain::getAverageGlobalUpot() const { return _globalUpot/_globalNumMolecules; }

void Domain::setCurrentTime(double curtime){ _currentTime = curtime;}

void Domain::advanceTime(double timestep){ _currentTime += timestep;}

double Domain::getCurrentTime(){ return _currentTime;}

vector<Component>& Domain::getComponents(){
	return _components; 
}

void Domain::addComponent(Component component){
	_components.push_back(component);
}

Comp2Param& Domain::getComp2Params(){
	return _comp2params; 
}

void Domain::calculateGlobalValues(
		DomainDecompBase* domainDecomp,
		ParticleContainer* particleContainer,
		bool collectThermostatVelocities,
		double Tfactor
		) {
	double Upot = _localUpot;
	double Virial = _localVirial;

	// To calculate Upot, Ukin and Pressure, intermediate values from all      
	// processes are needed. Here the         
	// intermediate values of all processes are summed up so that the root    
	// process can calculate the final values. to be able to calculate all     
	// values at this point, the calculation of the intermediate value sum_v2  
	// had to be moved from Thermostat to upd_postF and the final calculations  
	// of m_Ukin, m_Upot and Pressure had to be moved from Thermostat / upd_F  
	// to this point           
	
	/* FIXME stuff for the ensemble class */
	domainDecomp->collCommInit(2);
	domainDecomp->collCommAppendDouble(Upot);
	domainDecomp->collCommAppendDouble(Virial);
	domainDecomp->collCommAllreduceSum();
	Upot = domainDecomp->collCommGetDouble();
	Virial = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();

	/* FIXME: why should process 0 do this alone? 
	 * we should keep symmetry of all proccesses! */
	// Process 0 has to add the dipole correction:
	// m_UpotCorr and m_VirialCorr already contain constant (internal) dipole correction
	_globalUpot = Upot + _UpotCorr;
	_globalVirial = Virial + _VirialCorr;

	/*
	 * thermostat ID 0 represents the entire system
	 */

	map<int, unsigned long>::iterator thermit;
	if( _componentwiseThermostat )
	{
#ifndef NDEBUG
		global_log->debug() << "* applying a componentwise thermostat" << endl;
#endif
		this->_localThermostatN[0] = 0;
		this->_localRotationalDOF[0] = 0;
		this->_local2KETrans[0] = 0;
		this->_local2KERot[0] = 0;
		for(thermit = _localThermostatN.begin(); thermit != _localThermostatN.end(); thermit++)
		{
			if(thermit->first == 0) continue;
			this->_localThermostatN[0] += thermit->second;
			this->_localRotationalDOF[0] += this->_localRotationalDOF[thermit->first];
			this->_local2KETrans[0] += this->_local2KETrans[thermit->first];
			this->_local2KERot[0] += this->_local2KERot[thermit->first];
		}
	}

	for(thermit = _universalThermostatN.begin(); thermit != _universalThermostatN.end(); thermit++)
	{
		// number of molecules on the local process. After the reduce operation
		// num_molecules will contain the global number of molecules
		unsigned long numMolecules = _localThermostatN[thermit->first];
		double summv2 = _local2KETrans[thermit->first];
		assert(summv2 >= 0.0);
		unsigned long rotDOF = _localRotationalDOF[thermit->first];
		double sumIw2 = (rotDOF > 0)? _local2KERot[thermit->first]: 0.0;

		domainDecomp->collCommInit(4);
		domainDecomp->collCommAppendDouble(summv2);
		domainDecomp->collCommAppendDouble(sumIw2);
		domainDecomp->collCommAppendUnsLong(numMolecules);
		domainDecomp->collCommAppendUnsLong(rotDOF);
		domainDecomp->collCommAllreduceSum();
		summv2 = domainDecomp->collCommGetDouble();
		sumIw2 = domainDecomp->collCommGetDouble();
		numMolecules = domainDecomp->collCommGetUnsLong();
		rotDOF = domainDecomp->collCommGetUnsLong();
		domainDecomp->collCommFinalize();
		global_log->debug() << "[ thermostat ID " << thermit->first << "]\tN = " << numMolecules << "\trotDOF = " << rotDOF 
			<< "\tmv2 = " <<  summv2 << "\tIw2 = " << sumIw2 << endl;

		this->_universalThermostatN[thermit->first] = numMolecules;
		this->_universalRotationalDOF[thermit->first] = rotDOF;

		/* calculate the temperature of the entire system */
		if(numMolecules > 0)
			_globalTemperatureMap[thermit->first] =
				(summv2 + sumIw2) / (double)(3*numMolecules + rotDOF);
		else
			_globalTemperatureMap[thermit->first] = _universalTargetTemperature[thermit->first];

		double Ti = Tfactor * _universalTargetTemperature[thermit->first];
		if((Ti > 0.0) && !_universalNVE)
		{
			_universalBTrans[thermit->first] = pow(3.0*numMolecules*Ti / summv2, 0.4);
			if( sumIw2 == 0.0 ) 
				_universalBRot[thermit->first] = 1.0;
			else 
				_universalBRot[thermit->first] = pow(rotDOF*Ti / sumIw2, 0.4);
		}
		else
		{
			this->_universalBTrans[thermit->first] = 1.0;
			this->_universalBRot[thermit->first] = 1.0;
		}

		// heuristice handling of the unfortunate special case of an explosion in the system
		if( ( (_universalBTrans[thermit->first] < MIN_BETA) || (_universalBRot[thermit->first] < MIN_BETA) )
				&& (0 >= _universalSelectiveThermostatError) )
		{
			global_log->warning() << "Explosion warning (time t=" << _currentTime << ")." << endl;
			global_log->debug() << "Selective thermostat will be applied to set " << thermit->first
				<< " (beta_trans = " << this->_universalBTrans[thermit->first]
				<< ", beta_rot = " << this->_universalBRot[thermit->first] << "!)" << endl;
			int rot_dof;
			double Utrans, Urot;
			double limit_energy =  KINLIMIT_PER_T * Ti;
			double limit_rot_energy;
			double vcorr, Dcorr;
			Molecule* tM;
			for( tM = particleContainer->begin();
					tM != particleContainer->end();
					tM = particleContainer->next() )
			{
				Utrans = tM->U_trans();
				if(Utrans > limit_energy)
				{
					vcorr = sqrt(limit_energy / Utrans);
					global_log->debug() << ": v(m" << tM->id() << ") *= " << vcorr << endl;
					tM->scale_v(vcorr);
					tM->scale_F(vcorr);
				}

				rot_dof = _components[tM->componentid()].getRotationalDegreesOfFreedom();
				if(rot_dof > 0)
				{
					limit_rot_energy = 3.0*rot_dof * Ti;
					Urot = tM->U_rot();
					if(Urot > limit_rot_energy)
					{
						Dcorr = sqrt(limit_rot_energy / Urot);
						global_log->debug() << "D(m" << tM->id() << ") *= " << Dcorr << endl;
						tM->scale_D(Dcorr);
						tM->scale_M(Dcorr);
					}
				}
			}

			/* FIXME: Unnamed constant 3960... */
			if(3960 >= _universalSelectiveThermostatCounter)
			{
				if( _universalSelectiveThermostatWarning > 0 )
					_universalSelectiveThermostatError = _universalSelectiveThermostatWarning;
				if( _universalSelectiveThermostatCounter > 0 )
					_universalSelectiveThermostatWarning = _universalSelectiveThermostatCounter;
				_universalSelectiveThermostatCounter = 4000;
			}
			_universalBTrans[thermit->first] = 1.0;
			_universalBRot[thermit->first] = pow(this->_universalBRot[thermit->first], 0.0091);
		}
#ifdef NDEBUG
		if( (_universalSelectiveThermostatCounter > 0) &&
				((_universalSelectiveThermostatCounter % 20) == 10) )
#endif
			/* FIXME: why difference counters? */
			global_log->debug() << "counter " << _universalSelectiveThermostatCounter
				<< ",\t warning " << _universalSelectiveThermostatWarning
				<< ",\t error " << _universalSelectiveThermostatError << endl;

		if(collectThermostatVelocities && _universalUndirectedThermostat[thermit->first])
		{
			double sigv[3];
			for(int d=0; d < 3; d++)
				sigv[d] = _localThermostatDirectedVelocity[d][thermit->first];

			domainDecomp->collCommInit(3);
			for(int d=0; d < 3; d++) domainDecomp->collCommAppendDouble(sigv[d]);
			domainDecomp->collCommAllreduceSum();
			for(int d=0; d < 3; d++) sigv[d] = domainDecomp->collCommGetDouble();
			domainDecomp->collCommFinalize();

			for(int d=0; d < 3; d++)
			{
				_localThermostatDirectedVelocity[d][thermit->first] = 0.0;
				if(numMolecules > 0)
					_universalThermostatDirectedVelocity[d][thermit->first] = sigv[d] / numMolecules;
				else 
					_universalThermostatDirectedVelocity[d][thermit->first] = 0.0;
			}

#ifndef NDEBUG
			global_log->debug() << "* thermostat " << thermit->first
				<< " directed velocity: ("
				<< _universalThermostatDirectedVelocity[0][thermit->first]
				<< " / " << _universalThermostatDirectedVelocity[1][thermit->first]
				<< " / " << _universalThermostatDirectedVelocity[2][thermit->first] 
				<< ")" << endl;
#endif
		}

#ifndef NDEBUG
		global_log->debug() << "* Th" << thermit->first << " N=" << numMolecules
			<< " DOF=" << rotDOF + 3.0*numMolecules
			<< " Tcur=" << _globalTemperatureMap[thermit->first]
			<< " Ttar=" << _universalTargetTemperature[thermit->first]
			<< " Tfactor=" << Tfactor
			<< " bt=" << _universalBTrans[thermit->first]
			<< " br=" << _universalBRot[thermit->first] << "\n";
#endif
	}

	if(this->_universalSelectiveThermostatCounter > 0)
		this->_universalSelectiveThermostatCounter--;
	if(this->_universalSelectiveThermostatWarning > 0)
		this->_universalSelectiveThermostatWarning--;
	if(this->_universalSelectiveThermostatError > 0)
		this->_universalSelectiveThermostatError--;
}

void Domain::calculateThermostatDirectedVelocity(ParticleContainer* partCont)
{
	Molecule* tM;
	if(this->_componentwiseThermostat)
	{
		for( map<int, bool>::iterator thit = _universalUndirectedThermostat.begin();
				thit != _universalUndirectedThermostat.end();
				thit ++ )
		{
			if(thit->second)
				for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][thit->first] = 0.0;
		}
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			int cid = tM->componentid();
			int thermostat = this->_componentToThermostatIdMap[cid];
			if(this->_universalUndirectedThermostat[thermostat])
			{
				for(int d=0; d < 3; d++)
					_localThermostatDirectedVelocity[d][thermostat] += tM->v(d);
			}
		}
	}
	else if(this->_universalUndirectedThermostat[0])
	{
		for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][0] = 0.0;
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			for(int d=0; d < 3; d++)
				_localThermostatDirectedVelocity[d][0] += tM->v(d);
		}
	}
}

void Domain::calculateVelocitySums(ParticleContainer* partCont)
{
	Molecule* tM;
	if(this->_componentwiseThermostat)
	{
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			int cid = tM->componentid();
			int thermostat = this->_componentToThermostatIdMap[cid];
			this->_localThermostatN[thermostat]++;
			this->_localRotationalDOF[thermostat] += _components[cid].getRotationalDegreesOfFreedom();
			if(this->_universalUndirectedThermostat[thermostat])
			{
				tM->calculate_mv2_Iw2( this->_local2KETrans[thermostat],
						this->_local2KERot[thermostat],
						this->_universalThermostatDirectedVelocity[0][thermostat],
						this->_universalThermostatDirectedVelocity[1][thermostat],
						this->_universalThermostatDirectedVelocity[2][thermostat]  );
			}
			else
			{
				tM->calculate_mv2_Iw2(_local2KETrans[thermostat], _local2KERot[thermostat]);
			}
		}
	}
	else
	{
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			this->_localThermostatN[0]++;
			this->_localRotationalDOF[0] += _components[ tM->componentid() ].getRotationalDegreesOfFreedom();
			if(this->_universalUndirectedThermostat[0])
			{
				tM->calculate_mv2_Iw2( this->_local2KETrans[0],
						this->_local2KERot[0],
						this->_universalThermostatDirectedVelocity[0][0],
						this->_universalThermostatDirectedVelocity[1][0],
						this->_universalThermostatDirectedVelocity[2][0]  );
			}
			else
			{
				tM->calculate_mv2_Iw2(_local2KETrans[0], _local2KERot[0]);
			}
		}
		global_log->debug() << "      * N = " << this->_localThermostatN[0]
			<< "rotDOF = " << this->_localRotationalDOF[0] << "   mv2 = "
			<< _local2KETrans[0] << " Iw2 = " << _local2KERot[0] << endl;
	}
}

void Domain::writeCheckpoint( string filename, 
		ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp )
{
	domainDecomp->assertDisjunctivity(particleContainer);
	if(!this->_localRank)
	{
		ofstream checkpointfilestream(filename.c_str());
		checkpointfilestream << "mardyn trunk " << VERSION;
		checkpointfilestream << "\n";
		checkpointfilestream << " currentTime\t"  << this->_currentTime << "\n";
		checkpointfilestream << " Length\t" << setprecision(9) << _globalLength[0] << " " << _globalLength[1] << " " << _globalLength[2] << "\n";
		if(this->_componentwiseThermostat)
		{
			for( map<int, int>::iterator thermit = this->_componentToThermostatIdMap.begin();
					thermit != this->_componentToThermostatIdMap.end();
					thermit++ )
			{
				if(0 >= thermit->second) continue;
				checkpointfilestream << " CT\t" << 1+thermit->first
					<< "\t" << thermit->second << "\n";
			}
			for( map<int, double>::iterator Tit = this->_universalTargetTemperature.begin();
					Tit != this->_universalTargetTemperature.end();
					Tit++ )
			{
				if((0 >= Tit->first) || (0 >= Tit->second)) continue;
				checkpointfilestream << " ThT " << Tit->first << "\t" << Tit->second << "\n";
			}
		}
		else
		{
			checkpointfilestream << " Temperature\t" << _universalTargetTemperature[0] << endl;
		}
#ifndef NDEBUG
		checkpointfilestream << "# rho\t" << this->_globalRho << "\n";
		checkpointfilestream << "# rc\t" << particleContainer->getCutoff() << "\n";
		checkpointfilestream << "# rcT\t" << particleContainer->getTersoffCutoff() << "\n";
#endif
		/* by Stefan Becker: the output line "I ..." causes an error: the restart run does not start!!!
		if(this->_globalUSteps > 1)
		{
			checkpointfilestream << setprecision(13);
			checkpointfilestream << " I\t" << this->_globalUSteps << " "
				<< this->_globalSigmaU << " " << this->_globalSigmaUU << "\n";
			checkpointfilestream << setprecision(8);
		}
		*/
		checkpointfilestream << " NumberOfComponents\t" << _components.size() << endl;
		for(vector<Component>::const_iterator pos=_components.begin();pos!=_components.end();++pos){
			pos->write(checkpointfilestream);
		}
		unsigned int numperline=_components.size();
		unsigned int iout=0;
		for(vector<double>::const_iterator pos=_mixcoeff.begin();pos!=_mixcoeff.end();++pos){
			checkpointfilestream << *pos;
			iout++;
			// 2 parameters (xi and eta)
			if(iout/2>=numperline) {
				checkpointfilestream << endl;
				iout=0;
				--numperline;
			}
			else if(!(iout%2)) {
				checkpointfilestream << "\t";
			}
			else {
				checkpointfilestream << " ";
			}
		}
		checkpointfilestream << _epsilonRF << endl;
		map<unsigned, unsigned> componentSets = this->_universalPG->getComponentSets();
		for( map<unsigned, unsigned>::const_iterator uCSIDit = componentSets.begin();
				uCSIDit != componentSets.end();
				uCSIDit++ )
		{ 
			if(uCSIDit->first > 100) continue;
			checkpointfilestream << " S\t" << 1+uCSIDit->first << "\t" << uCSIDit->second << "\n";
		}
		map<unsigned, double> tau = this->_universalPG->getTau();
		for( map<unsigned, double>::const_iterator gTit = tau.begin();
				gTit != tau.end();
				gTit++ )
		{
			unsigned cosetid = gTit->first;
			double* ttargetv = this->_universalPG->getTargetVelocity(cosetid);
			double* tacc = this->_universalPG->getAdditionalAcceleration(cosetid);
			checkpointfilestream << " A\t" << cosetid << "\t"
				<< ttargetv[0] << " " << ttargetv[1] << " " << ttargetv[2] << "\t"
				<< gTit->second << "\t"
				<< tacc[0] << " " << tacc[1] << " " << tacc[2] << "\n";
			delete ttargetv;
			delete tacc;
		}
		for( map<int, bool>::iterator uutit = this->_universalUndirectedThermostat.begin();
				uutit != this->_universalUndirectedThermostat.end();
				uutit++ )
		{
			if(0 > uutit->first) continue;
			if(uutit->second) checkpointfilestream << " U\t" << uutit->first << "\n";
		}
		checkpointfilestream << " NumberOfMolecules\t" << _globalNumMolecules << endl;

		checkpointfilestream << " MoleculeFormat\t" << "ICRVQD" << endl;
		checkpointfilestream.close();
	}

	domainDecomp->writeMoleculesToFile(filename, particleContainer); 
}

void Domain::initParameterStreams(double cutoffRadius, double cutoffRadiusLJ){
	_comp2params.initialize(_components, _mixcoeff, _epsilonRF, cutoffRadius, cutoffRadiusLJ); 
}

void Domain::initFarFieldCorr(double cutoffRadius, double cutoffRadiusLJ) {
	double UpotCorrLJ=0.;
	double VirialCorrLJ=0.;
	double MySelbstTerm=0.;
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
							UpotCorrLJ+=fac*(TICCu(-6,cutoffRadiusLJ,sig2)-TICCu(-3,cutoffRadiusLJ,sig2));
							VirialCorrLJ+=fac*(TICCv(-6,cutoffRadiusLJ,sig2)-TICCv(-3,cutoffRadiusLJ,sig2));
						}
						else if(tau1!=0. && tau2!=0.)
						{
							UpotCorrLJ += fac*( TISSu(-6,cutoffRadiusLJ,sig2,tau1,tau2)
									- TISSu(-3,cutoffRadiusLJ,sig2,tau1,tau2) );
							VirialCorrLJ += fac*( TISSv(-6,cutoffRadiusLJ,sig2,tau1,tau2)
									- TISSv(-3,cutoffRadiusLJ,sig2,tau1,tau2) );
						}
						else {
							if(tau2==0.) 
								tau2=tau1;
							UpotCorrLJ+=fac*(TICSu(-6,cutoffRadiusLJ,sig2,tau2)-TICSu(-3,cutoffRadiusLJ,sig2,tau2));
							VirialCorrLJ+=fac*(TICSv(-6,cutoffRadiusLJ,sig2,tau2)-TICSv(-3,cutoffRadiusLJ,sig2,tau2));
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

void Domain::setupProfile(unsigned xun, unsigned yun, unsigned zun)
{
	this->_universalNProfileUnits[0] = xun;
	this->_universalNProfileUnits[1] = yun;
	this->_universalNProfileUnits[2] = zun;
	// comment by Stefan Becker: the next operation is always carried out.
	// In case the profile is to be cylindircal, the method "void SesDrop()"is called afterwards and "_universalInvProfileUnit[d]" is changed accordingly.
	for(unsigned d = 0; d < 3; d++)
	{
		_universalInvProfileUnit[d] = _universalNProfileUnits[d] / _globalLength[d];
	}
	this->resetProfile();
}

void Domain::collectProfile(DomainDecompBase* dode)
{
	unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
		* this->_universalNProfileUnits[2];
	dode->collCommInit(6*unIDs);
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		dode->collCommAppendLongDouble(this->_localNProfile[unID]);
		for(int d=0; d<3; d++)
			dode->collCommAppendLongDouble(_localvProfile[d][unID]);
		dode->collCommAppendLongDouble(this->_localDOFProfile[unID]);
		dode->collCommAppendLongDouble(_localKineticProfile[unID]);
	}
	dode->collCommAllreduceSum();
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		_globalNProfile[unID] = (double)dode->collCommGetLongDouble();
		for(int d=0; d<3; d++)
			this->_globalvProfile[d][unID]
				= (double)dode->collCommGetLongDouble();
		this->_globalDOFProfile[unID]
			= (double)dode->collCommGetLongDouble();
		this->_globalKineticProfile[unID]
			= (double)dode->collCommGetLongDouble();
	}
	dode->collCommFinalize();
}

void Domain::outputProfile(const char* prefix)
{
	if(this->_localRank) return;

	string vzpryname(prefix);
	string Tpryname(prefix);
	string rhpryname(prefix);
	rhpryname += ".rhpry";
	vzpryname += ".vzpry";
	Tpryname += ".Tpry";
	ofstream rhpry(rhpryname.c_str());
	ofstream vzpry(vzpryname.c_str());
	ofstream Tpry(Tpryname.c_str());
	if (!(vzpry && Tpry && rhpry))
	{
		return;
	}
	rhpry.precision(4);
	rhpry << "# y\trho\ttotal DOF\n# \n";
	vzpry.precision(4);
	vzpry << "# y\tvz\tv\n# \n";
	Tpry.precision(5);
	Tpry << "# y\t2Ekin/#DOF\n# \n";

	double layerVolume = this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]
		/ this->_universalNProfileUnits[1];
	for(unsigned y = 0; y < this->_universalNProfileUnits[1]; y++)
	{
		double yval = (y + 0.5) / this->_universalInvProfileUnit[1];

		long double Ny = 0.0;
		long double DOFy = 0.0;
		long double twoEkiny = 0.0;
		long double velocitysumy[3];
		for(unsigned d = 0; d < 3; d++) velocitysumy[d] = 0.0;
		for(unsigned x = 0; x < this->_universalNProfileUnits[0]; x++)
		{
			for(unsigned z = 0; z < this->_universalNProfileUnits[2]; z++)
			{
				unsigned unID = x * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
					+ y * this->_universalNProfileUnits[2] + z;
				Ny += this->_globalNProfile[unID];
				DOFy += this->_globalDOFProfile[unID];
				twoEkiny += this->_globalKineticProfile[unID];
				for(unsigned d = 0; d < 3; d++) velocitysumy[d] += this->_globalvProfile[d][unID];
			}
		}

		if(Ny >= 64.0)
		{
			double vvdir = 0.0;
			for(unsigned d = 0; d < 3; d++)
			{
				double vd = velocitysumy[d] / Ny;
				vvdir += vd*vd;
			}
			rhpry << yval << "\t" << (Ny / (layerVolume * this->_globalAccumulatedDatasets))
				<< "\t" << DOFy << "\n";
			vzpry << yval << "\t" << (velocitysumy[2] / Ny) << "\t" << sqrt(vvdir) << "\n";
			Tpry << yval << "\t" << (twoEkiny / DOFy) << "\n";
		}
		else
		{
			rhpry << yval << "\t0.000\t" << DOFy << "\n";
		}
	}

	rhpry.close();
	vzpry.close();
	Tpry.close();
}

void Domain::resetProfile()
{
	unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
		* this->_universalNProfileUnits[2];
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		this->_localNProfile[unID] = 0.0;
		this->_globalNProfile[unID] = 0.0;
		for(int d=0; d<3; d++)
		{
			this->_localvProfile[d][unID] = 0.0;
			this->_globalvProfile[d][unID] = 0.0;
		}
		this->_localDOFProfile[unID] = 0.0;
		this->_globalDOFProfile[unID] = 0.0;
		this->_localKineticProfile[unID] = 0.0;
		this->_globalKineticProfile[unID] = 0.0;
	}
	this->_globalAccumulatedDatasets = 0;
}


void Domain::Nadd(unsigned cid, int N, int localN)
{
	this->_components[cid].incNumMolecules(N);
	this->_globalNumMolecules += N;
	this->_localRotationalDOF[0] += localN * _components[cid].getRotationalDegreesOfFreedom();
	this->_universalRotationalDOF[0] += N * _components[cid].getRotationalDegreesOfFreedom();
	if( (this->_componentwiseThermostat)
			&& (this->_componentToThermostatIdMap[cid] > 0) )
	{
		int thid = this->_componentToThermostatIdMap[cid];
		this->_localThermostatN[thid] += localN;
		this->_universalThermostatN[thid] += N;
		this->_localRotationalDOF[thid] += localN * _components[cid].getRotationalDegreesOfFreedom();
		this->_universalRotationalDOF[thid] += N * _components[cid].getRotationalDegreesOfFreedom();
	}
	this->_localThermostatN[0] += localN;
	this->_universalThermostatN[0] += N;
	this->_localRotationalDOF[0] += localN * _components[cid].getRotationalDegreesOfFreedom();
	this->_universalRotationalDOF[0] += N * _components[cid].getRotationalDegreesOfFreedom();
}

void Domain::evaluateRho(
		unsigned long localN, DomainDecompBase* domainDecomp
		) {
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(localN);
	domainDecomp->collCommAllreduceSum();
	this->_globalNumMolecules = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();

	this->_globalRho = this->_globalNumMolecules /
		(this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]);
}

void Domain::setTargetTemperature(int thermostat, double targetT)
{
	if(thermostat < 0)
	{
		global_log->warning() << "Warning: thermostat \'" << thermostat << "\' (T = "
			<< targetT << ") will be ignored." << endl;
		return;
	}

	this->_universalTargetTemperature[thermostat] = targetT;
	if(!(this->_universalUndirectedThermostat[thermostat] == true))
		this->_universalUndirectedThermostat[thermostat] = false;

	/* FIXME: Substantial change in program behaviour! */
	if(thermostat == 0) {
		global_log->warning() << "Disabling the component wise thermostat!" << endl;
		disableComponentwiseThermostat();
	}
	if(thermostat >= 1) {
		if( ! _componentwiseThermostat ) {
			/* FIXME: Substantial change in program behaviour! */
			global_log->warning() << "Enabling the component wise thermostat!" << endl;
			_componentwiseThermostat = true;
			_universalTargetTemperature.erase(0);
			_universalUndirectedThermostat.erase(0);
			for(int d=0; d < 3; d++) this->_universalThermostatDirectedVelocity[d].erase(0);
			for( vector<Component>::iterator tc = this->_components.begin();
					tc != this->_components.end();
					tc ++ )
				if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0))
					this->_componentToThermostatIdMap[ tc->ID() ] = -1;
		}
	}
}

void Domain::enableComponentwiseThermostat()
{
	if(this->_componentwiseThermostat) return;

	this->_componentwiseThermostat = true;
	this->_universalTargetTemperature.erase(0);
	for( vector<Component>::iterator tc = this->_components.begin();
			tc != this->_components.end();
			tc ++ )
		if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0))
			this->_componentToThermostatIdMap[ tc->ID() ] = -1;
}

void Domain::enableUndirectedThermostat(int tst)
{
	this->_universalUndirectedThermostat[tst] = true;
	for(int d=0; d < 3; d++)
	{
		this->_universalThermostatDirectedVelocity[d][tst] = 0.0;
		this->_localThermostatDirectedVelocity[d][tst] = 0.0;
	}
}

void Domain::setGlobalTemperature(double temp)
{
	this->disableComponentwiseThermostat();
	this->_universalTargetTemperature[0] = temp;
}

vector<double> & Domain::getmixcoeff() { return _mixcoeff; }

double Domain::getepsilonRF() const { return _epsilonRF; }

void Domain::setepsilonRF(double erf) { _epsilonRF = erf; }

unsigned long Domain::getglobalNumMolecules() const { return _globalNumMolecules; }

void Domain::setglobalNumMolecules(unsigned long glnummol) { _globalNumMolecules = glnummol; }

unsigned long Domain::getinpversion(){ return _inpversion;}

void Domain::setinpversion(unsigned long inpv){ _inpversion = inpv;}

double Domain::getglobalRho(){ return _globalRho;}

void Domain::setglobalRho(double grho){ _globalRho = grho;}

unsigned long Domain::getglobalRotDOF()
{
	return this->_universalRotationalDOF[0]; 
}

void Domain::setglobalRotDOF(unsigned long grotdof)
{
	this->_universalRotationalDOF[0] = grotdof;
}

void Domain::setGlobalLength(int index, double length) {
	_globalLength[index] = length;
}

void Domain::record_cv()
{
	if(_localRank != 0) return;

	this->_globalUSteps ++;
	this->_globalSigmaU += this->_globalUpot;
	this->_globalSigmaUU += this->_globalUpot*_globalUpot;
}

double Domain::cv()
{
	if((_localRank != 0) || (_globalUSteps == 0)) return 0.0;

	double id = 1.5 + 0.5*_universalRotationalDOF[0]/_globalNumMolecules;
	double conf = (_globalSigmaUU - _globalSigmaU*_globalSigmaU/_globalUSteps)
		/ (_globalUSteps * _globalNumMolecules * _globalTemperatureMap[0] * _globalTemperatureMap[0]);

	return id + conf;
}

//! methods implemented by Stefan Becker <stefan.becker@mv.uni-kl.de>
// the following two methods are used by the MmspdWriter (writing the output file in a format used by MegaMol)
double Domain::getSigma(unsigned cid, unsigned nthSigma){
  return _components[cid].getSigma(nthSigma);
}
unsigned Domain::getNumberOfComponents(){
  return _components.size();
}


// author: Stefan Becker, counterpart of method setupProfile(xun,yun,zun) => to be used for profile in cylindrical coordinates
// the radial component is located in the x,z-plane, the axial component is parallel to the y-direction.
// In the radial direction the spacing is linear in R^2 in order to obtain equal areas with increasing radius.
// @TODO: in the input file the token "profile" immedeatly causes a call of "setupProfile(xun,yun,zun)". This in turn sets up a cubic grid
// by default.
// Proposal: an enquiry that checks wheter a cubic, or cylindrical, sperical, etc. grid has to be set up. This first requires a check of
// the input file (what additional token is set that controls the kind of grid beeing used). => code more modular
void Domain::sesDrop(){
	this->_universalCylindricalGeometry = true;

	// R^2_max: (squared) maximum radius up to which the profile is recorded
	// (otherwise erroneous densities recorded at the boundary of the RECTANGULAR simulation box)
	double minXZ = this->_globalLength[0];
	if(this->_globalLength[2]<minXZ){
		minXZ = this->_globalLength[2];
	}
	this->_universalR2max = 0.25*minXZ*minXZ;

	// origin of the cylindrical coordinate system
	this->_universalCentre[0] = 0.5*this->_globalLength[0];
	this->_universalCentre[1] = 0;
	this->_universalCentre[2] = 0.5*this->_globalLength[2];

	_universalInvProfileUnit[0] = this->_universalNProfileUnits[0]/(2*M_PI);                   // delta_phi
	_universalInvProfileUnit[1] = this->_universalNProfileUnits[1]/(this->_universalR2max);  // delta_R^2
	_universalInvProfileUnit[2] = this->_universalNProfileUnits[2]/(this->_globalLength[1]); // delta_H
	global_log->info() << "\nInv Profile unit for sessile drop: (phi,R2,H) = (" << _universalInvProfileUnit[0] <<", " <<_universalInvProfileUnit[1]<<", "<<_universalInvProfileUnit[2]<<") \n";
}

// author: Stefan Becker. Method called by Simulation::output() in order to decide wheter or not a cylindrical profile is to be written out,
//i.e. wheter the method outputCylProfile() (isCylindrical==true) or the method outputProfile() (isCylindrical==false) is called.
bool Domain::isCylindrical(){
	return this->_universalCylindricalGeometry;
}

void Domain::considerComponentInProfile(unsigned cidmax)
{
	// default: false
	for(unsigned  i = 0; i<=_components.size(); i++){
	this->_universalProfiledComponents[i] = false;
	}
	
	for(unsigned i = 0; i<=cidmax; i++){
	this->_universalProfiledComponents[i] = true;
#ifndef NDEBUG
	global_log->info() << "profiled Compoenent: "<< i << "\n";
#endif
	}
}

void Domain::considerComponentForYShift(unsigned cidMin, unsigned cidMax){
    for (unsigned i = 0; i <= cidMax; i++) _componentForYShift[i] = false;
    for (unsigned i = 0; i <= cidMax-cidMin; i++) _componentForYShift[cidMin+i] = true;
#ifndef NDEBUG
    global_log->info() << "Components used for determination of shift in y-direction:\n";
    for(unsigned i = 0; i <= cidMax; i++)
      global_log->info() << "component id "<< i << "\t value: " << _componentForYShift[i] << "\n";
#endif
}


void Domain::sYOffset(double in_yOff){
  _yOff = in_yOff;
}

// author: Stefan Becker, method to determine the integer value of unID, extra method in order to not inflating the method "recordProfile()"
// method only matched for the case of a cylindrical profile (i.e. sessile drop)
long int Domain::unID(double qx, double qy, double qz){
		int xun,yun,zun;// (xun,yun,zun): bin number in a special direction, e.g. yun==5 corresponds to the 5th bin in the radial direction, 
		long int unID;	// as usual
		double xc,yc,zc; // distance of a particle with respect to the origin of the cylindrical coordinate system

		unID = -1; // initialization, causes an error message, if unID is not calculated in this method but used in record profile

		xc = qx - this->_universalCentre[0];
	    yc = qy - this->_universalCentre[1];
	    zc = qz - this->_universalCentre[2];

	    // transformation in polar coordinates
	    double R2 = xc*xc + zc*zc;
	    double phi = asin(zc/sqrt(R2)) + ((xc>=0.0) ? 0:M_PI);
	    if(phi<0.0) {phi = phi + 2.0*M_PI;}

	    xun = (int)floor(phi * this->_universalInvProfileUnit[0]);   // bin no. in phi-direction
	    yun = (int)floor(R2 *  this->_universalInvProfileUnit[1]);   // bin no. in R-direction
	    zun = (int)floor(yc *  this->_universalInvProfileUnit[2]);   // bin no. in H-direction

	    if((xun >= 0) && (yun >= 0) && (zun >= 0) &&
	          (xun < (int)_universalNProfileUnits[0]) && (yun < (int)_universalNProfileUnits[1]) && (zun < (int)_universalNProfileUnits[2]))
	       {
	          unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
	               + yun * this->_universalNProfileUnits[2] + zun;
	       }
	       else 
	       {
	          global_log->error() << "Severe error!! Invalid profile unit (" << xun << " / " << yun << " / " << zun << ").\n\n";
	          global_log->error() << "Coordinates (" << qx << " / " << qy << " / " << qz << ").\n";
		  global_log->error() << "unID = " << unID << "\n";
	          //exit(707);
	       }
	       return unID;
}

void Domain::recordProfile(ParticleContainer* molCont)
{
	int cid;
	unsigned xun, yun, zun;
	long int unID;
	double mv2, Iw2;
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
		cid = thismol->componentid();
		if(this->_universalProfiledComponents[cid])
		{
// by Stefan Becker: enquiry if(_universalCylindricalGeometry) ... implemented
// possible???: if no cylindrical profile is recorded => permanent calculation (before the "if..." ) of "distFor_unID" slows down the code
// if the profile is however recorded in cylindrical coordinates, the current implementation is fast (?)   => solution?
			double distFor_unID = pow((thismol->r(0)- this->_universalCentre[0]),2.0) + pow((thismol->r(2)- this->_universalCentre[2]),2.0);
			if(this->_universalCylindricalGeometry && distFor_unID <= this->_universalR2max){
				unID = this->unID(thismol->r(0), thismol->r(1), thismol->r(2));
				if(unID < 0){
				  global_log->error() << "ERROR: in Domain.cpp/recordProfile: unID < 0!!! Was not calculated in unID() but initialized with -1!\n";
				}
			}
			else if(!this->_universalCylindricalGeometry){
			xun = (unsigned)floor(thismol->r(0) * this->_universalInvProfileUnit[0]);
			yun = (unsigned)floor(thismol->r(1) * this->_universalInvProfileUnit[1]);
			zun = (unsigned)floor(thismol->r(2) * this->_universalInvProfileUnit[2]);
			unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
				+ yun * this->_universalNProfileUnits[2] + zun;
			}
			else{
				continue;
			}

// @TODO: (by Stefan Becker)  differentiation of _localNProfile by the component number cid => _localNProfile[cid][unID]!!!
			this->_localNProfile[unID] += 1.0;
			for(int d=0; d<3; d++) this->_localvProfile[d][unID] += thismol->v(d);
			this->_localDOFProfile[unID] += 3.0 + (long double)(_components[cid].getRotationalDegreesOfFreedom());

			// record _twice_ the total (ordered + unordered) kinetic energy
			// by Stefan Becker: ?! What happens in the following 4 lines?! calculate_mv2_Iw2() does not return anything, and _localKineticProfile[unID] remains zero at all times?!
			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2);
			this->_localKineticProfile[unID] += mv2+Iw2;
		}
	}
	this->_globalAccumulatedDatasets++;
}

/*
 *By Stefan Becker: 
 *realign tool, borrowed from Martin Horsch: it effects that the center of mass of all the particles (with respect to the (x,z)-plane)
 *corresponds with the box center (in the (x,z)plane). With respect to the y-direction, the shift is carried out so 
 *that the wall remains at the initial position.
 *Basically the current centre of mass of all the particles is determined and then all the particles are moved towards x==0.5*Lx and z==0*Lz
 *The method is frequently applied when the specified span of time (i.e. timesteps) is reached. 
 *Part of this tool in Domain.cpp: (i) determineShift() => determines the vector by which the particles are shifted
 *                                 (ii) realign() => carries out the actual realignment / shift
 * when this method is called, the halo should NOT be present
 */
void Domain::determineXZShift( DomainDecompBase* domainDecomp, ParticleContainer* molCont,
			     double fraction) 
{
   double localBalance[3]; // localBalance[1] not considered, employed only for not confusing: index 0 => x-direction, index 1 => y-direction, index 2 => z-direction
   for(unsigned d = 0; d < 3; d ++) localBalance[d] = 0.0; // initialising the array by zeros
   double localMass = 0.0;
   int cid;
      
   for(Molecule* tm = molCont->begin(); tm != molCont->end(); tm = molCont->next())
   {
      cid = tm->componentid();
      if(_universalProfiledComponents[cid])
      {
	double tmass = tm->gMass();
	localMass += tmass;
	for(unsigned d = 0; d < 3; d=d+2){ // counting d=d+2 causes the value d==1 (i.e. y-direction) to be skipped, this is performed by the method determineYShift()
	  localBalance[d] += tm->r(d) * tmass;
	}
      }
   }
   
   for(unsigned d = 0; d < 3; d++) _globalRealignmentBalance[d] = localBalance[d];
   _globalRealignmentMass[0] = localMass;
   // determining the global values by the use of collectiveCommunication
   domainDecomp->collCommInit(4);
   for(unsigned d = 0; d < 3; d++)  domainDecomp->collCommAppendDouble(_globalRealignmentBalance[d]);
   domainDecomp->collCommAppendDouble(_globalRealignmentMass[0]);
   domainDecomp->collCommAllreduceSum();
   for(unsigned d = 0; d < 3; d++)  _globalRealignmentBalance[d] = domainDecomp->collCommGetDouble();
   _globalRealignmentMass[0] = domainDecomp->collCommGetDouble();
   domainDecomp->collCommFinalize();
   
   for(unsigned short d = 0; d < 3; d=d+2){
      _universalRealignmentMotion[d]
         = -fraction*((_globalRealignmentBalance[d] / _globalRealignmentMass[0]) - 0.5*_globalLength[d]);
   }  
}

void Domain::determineYShift( DomainDecompBase* domainDecomp, ParticleContainer* molCont,
			     double fraction){
			        // keep in mind: variables declared as static are initialized only ONCE!
   static unsigned timesCalled = 0;
   static double initialCentreOfMassY;
   double localBalance = 0.0;
   double localMass = 0.0;
   int cid;
   for(Molecule* tm = molCont->begin(); tm != molCont->end(); tm = molCont->next())
   {
     cid = tm->componentid();
     if(_componentForYShift[cid])
     {  
	double tmass = tm->gMass();
	localMass += tmass;
	localBalance += (tm->r(1) - _globalLength[1]*floor(2.0* tm->r(1)/ _globalLength[1]))*tmass; // PBC
     }
   }
   _globalRealignmentBalance[1] = localBalance;
   _globalRealignmentMass[1] = localMass;
   // determining the global values by the use of collectiveCommunication
   domainDecomp->collCommInit(2);
   domainDecomp->collCommAppendDouble(_globalRealignmentBalance[1]);
   domainDecomp->collCommAppendDouble(_globalRealignmentMass[1]);
   domainDecomp->collCommAllreduceSum();
   _globalRealignmentBalance[1] = domainDecomp->collCommGetDouble();
   _globalRealignmentMass[1] = domainDecomp->collCommGetDouble();
   domainDecomp->collCommFinalize();
   
   /* 
   The centre of mass of the wall (solely the wall!) is determined once at the beginning of the simulation. 
   Then the realignment with respect to the y-direction is always carried out so that the wall remains at the initial posistion 
   (initial y-coordinate)
   */
   if (!timesCalled){
     initialCentreOfMassY = _globalRealignmentBalance[1] / _globalRealignmentMass[1];
   }
   timesCalled++;
   //global_log->info() <<"initialCentreOfMassY = " << initialCentreOfMassY << endl;
   // the realignment motion in y-direction so that the wall is always at the bottom of the simulation box
   _universalRealignmentMotion[1] = -fraction*((_globalRealignmentBalance[1] / _globalRealignmentMass[1]) - initialCentreOfMassY);
}

/*
 * By Stefa Becker, see above comment on determineShift().
 * this method REQUIRES the presence of the halo
 */
void Domain::realign(
   ParticleContainer* molCont
) {
   if(!this->_localRank)
   {
      global_log->info() << "Centre of mass: (" << _globalRealignmentBalance[0]/_globalRealignmentMass[0] << " / " << _globalRealignmentBalance[1]/_globalRealignmentMass[1] << " / " << _globalRealignmentBalance[2]/_globalRealignmentMass[0] << ") "
			 << "=> adjustment: (" << _universalRealignmentMotion[0] << ", " << _universalRealignmentMotion[1] << ", " << _universalRealignmentMotion[2] << ").\n";
      // cout << "If required, shifting by +- (" << _globalLength[0] << ", " << _globalLength[1] << ", " << _globalLength[2] << ").\n";
   }
   for(Molecule* tm = molCont->begin(); tm != molCont->end(); tm = molCont->next())
   {
     for(unsigned short d=0; d<3; d++){
       tm->move(d, _universalRealignmentMotion[d]);
      // tm->move(_universalRealignmentMotion, _globalLength);
     }
   }
}


/* by Stefan Becker, borrowed from Matrin Horsch
 * method cancelling the net moment
 */
void Domain::cancelMomentum(
     DomainDecompBase* domainDecomp,
     ParticleContainer* molCont
) {
   double localMomentum[3];
   for(unsigned d = 0; d < 3; d++) localMomentum[d] = 0.0;
   for(Molecule* tm = molCont->begin(); tm != molCont->end(); tm = molCont->next())
   {
      double tmass = tm->gMass();
      for(unsigned d = 0; d < 3; d++)
         localMomentum[d] += tm->v(d) * tmass;
   }
   double globalMomentum[3];
   for(unsigned d = 0; d < 3; d++) globalMomentum[d] = localMomentum[d];
   domainDecomp->collCommInit(3);
   for(unsigned short d = 0; d<3; d++)domainDecomp->collCommAppendDouble(globalMomentum[d]);
   domainDecomp->collCommAllreduceSum();
   for(unsigned short d = 0; d < 3; d++) globalMomentum[d]  = domainDecomp->collCommGetDouble();
   for(unsigned short d = 0; d < 3; d++) globalMomentum[d] /= _globalNumMolecules;
   if(!this->_localRank)
      global_log->info() << "Average momentum: (" << globalMomentum[0] << ", " << globalMomentum[1] << ", " << globalMomentum[2] << ") => removing.\n";
   for(Molecule* tm = molCont->begin(); tm != molCont->end(); tm = molCont->next())
   {
      double tmass = tm->gMass();
      tm->vsub(globalMomentum[0]/tmass, globalMomentum[1]/tmass, globalMomentum[2]/tmass);
   }
}


// author: Stefan Becker, method called in the case of a density profile established in cylindrical coordinates. Counterpart of outputProfile(...).
// reason for a sperate method (in addition to "outputProfile(...)"): method neatly matched to the particular needs of the (cylindrical density) profile, otherwise outpuProfile would be inflated, structure became too compilcated.
void Domain::outputCylProfile(const char* prefix){

	if(this->_localRank) return;

	   unsigned IDweight[3];
	   IDweight[0] = this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2];
	   IDweight[1] = this->_universalNProfileUnits[2];
	   IDweight[2] = 1;

	   for( map<unsigned, bool>::iterator pcit = _universalProfiledComponents.begin();    // Loop ueber alle Komponenten
	   	           pcit != _universalProfiledComponents.end();
	   	           pcit++ )
	   	      {
	   	         if(!pcit->second) continue;  // ->second weist auf den key-value von map, d.h. den bool-Wert => falls falsche id, d.h. hiervon ist kein Profil zu erstellen => naechste Schleife

			 // density profile
	   	         string rhoProfName(prefix);
	   	         rhoProfName += ".rhpr";
			 // temperature profile
			 string tmpProfName(prefix);
			 tmpProfName += ".Tpr";

	   	         if(!this->_universalCylindricalGeometry)
	   	         {
	   	           global_log->error() << "Incorrect call of method \"outputCylProfile()\" !";
	   	         }

	   	         ofstream* rhoProf = new ofstream(rhoProfName.c_str());
			 ofstream* tmpProf = new ofstream(tmpProfName.c_str());

	   	         if (!(*rhoProf ) || !(*tmpProf)) // geaendert durch M. Horsch, by Stefan Becker: wozu?
	   	         {
	   	            return;
	   	         }
	   	         rhoProf->precision(6);
			 tmpProf->precision(6);

	    //##########################################################################################################################################			 
			 // density profile: actual writing procedure 
	   	         double segmentVolume; // volume of a single bin, in a0^3 (atomic units)
	   	      	 segmentVolume = M_PI/this->_universalInvProfileUnit[1]/this->_universalInvProfileUnit[2]/this->_universalNProfileUnits[0];  // adapted to cylindrical coordinates
			  
			 *rhoProf << "//Local profile of the number density. Output file generated by the \"outputCylProfile\" method, located in Domain.cpp. \n";
	   	         *rhoProf << "//local density profile: Each matrix corresponds to a single value of \"phi_i\", measured in [rad]\n";
	   	         *rhoProf << "//one single matrix of the local number density rho'(phi_i;r_i',h_i') \n//      | r_i'\n//---------------------\n//  h_i'| rho'(r_i',h_i')\n//      | \n";
	   	         *rhoProf << "// T' \t sigma_ii' \t eps_ii' \t yOffset \t DELTA_phi \t DELTA_r2' \t DELTA_h' \t quantities in atomic units are denoted by an apostrophe '\n";
	   	         *rhoProf << this->_universalTargetTemperature[_components.size()] <<"\t"<<_components[0].getSigma(0)<<"\t"<<_components[0].getEps(0)<<"\t";
	   	         *rhoProf << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	   	         // info: getSigma() und getEps() implementiert in Component.h
	   	         // end of header, start of the data-part of the density file
	   	         for(unsigned n_phi = 0; n_phi < this->_universalNProfileUnits[0]; n_phi++)
	   	         {
	   	        	 *rhoProf <<"> "<< (n_phi+0.5)/this->_universalInvProfileUnit[0] <<"\n0 \t";
	   	        	 for(unsigned n_r2 = 0; n_r2 < this->_universalNProfileUnits[1]; n_r2++){
	   	        		 *rhoProf << sqrt(n_r2/this->_universalInvProfileUnit[1])+0.5*sqrt(this->_universalInvProfileUnit[1])<<"  \t"; // Eintragen der radialen Koordinaten r_i in Header
	   	        	 }
	   	        	 *rhoProf << "\n";
	   	        	for(unsigned n_h = 0; n_h < this->_universalNProfileUnits[2]; n_h++)
	   	        	{

	   	        		double hval = (n_h + 0.5) / this->_universalInvProfileUnit[2];
	   	        		*rhoProf << hval<< "  \t";
	   	        		for(unsigned n_r2 = 0; n_r2< this->_universalNProfileUnits[1]; n_r2++)
	   	        		{
	   	        			unsigned unID = n_phi * IDweight[0] + n_r2 * IDweight[1]
	   	        				                                    + n_h * IDweight[2];
	   	        		   	double rho_loc = this->_globalNProfile[unID] / (segmentVolume * this->_globalAccumulatedDatasets);
	   	        		   	*rhoProf << rho_loc << "\t";
	   	        		}
	   	        		*rhoProf << "\n";
	   	        	}
	   	         }
	   	         rhoProf->close();
	   	         delete rhoProf;
 
	    //##########################################################################################################################################
	   	         // temperature profile: actual writing procedure
			 *tmpProf << "//Local temperature profile generated by the \"outputCylProfile\" method.\n";
			 *tmpProf << "//Temperature expressed by 2Ekin/#DOF\n";
			 *tmpProf << "//one single matrix of the local temperature T(r_i,h_i) \n//      | r_i\n//---------------------\n//  h_i| T(r_i,h_i)\n//      | \n";
	   	         *tmpProf << "// T \t sigma_ii \t eps_ii \t yOffset \t DELTA_phi \t DELTA_r2 \t DELTA_h \n";
	   	         *tmpProf << this->_universalTargetTemperature[_components.size()] <<"\t"<<_components[0].getSigma(0)<<"\t"<<_components[0].getEps(0)<<"\t";
	   	         *tmpProf << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	   	         *tmpProf <<"> "<< 0.5/this->_universalInvProfileUnit[0] <<"\n0 \t";
			 for(unsigned n_r2 = 0; n_r2 < this->_universalNProfileUnits[1]; n_r2++){
	   	        		 *tmpProf << sqrt(n_r2/this->_universalInvProfileUnit[1])+0.5*sqrt(this->_universalInvProfileUnit[1])<<"  \t"; // Eintragen der radialen Koordinaten r_i
	   	          }
	   	        	 *tmpProf << "\n";
	   	        	for(unsigned n_h = 0; n_h < this->_universalNProfileUnits[2]; n_h++)
	   	        	{
	   	        		double hval = (n_h + 0.5) / this->_universalInvProfileUnit[2];
	   	        		*tmpProf << hval<< "  \t";
	   	        		for(unsigned n_r2 = 0; n_r2< this->_universalNProfileUnits[1]; n_r2++)
	   	        		{
					  long double DOFc = 0.0;
					  long double twoEkinc = 0.0;
					  for(unsigned n_phi = 0; n_phi < this->_universalNProfileUnits[0]; n_phi++){
					     unsigned unID = n_phi * IDweight[0] + n_r2 * IDweight[1]
	   	        				                                    + n_h * IDweight[2];
					      DOFc += this->_globalDOFProfile[unID];
					      twoEkinc += this->_globalKineticProfile[unID];
					   }
					   if(DOFc == 0.0){
					     *tmpProf << 0 << "\t";
					   }
					   else{
					  *tmpProf << (twoEkinc/DOFc ) << "\t";
					   }
	   	        		}
	   	        		*tmpProf << "\n";
	   	        	}
			  tmpProf->close();
			  delete tmpProf;
	   	      }
}
