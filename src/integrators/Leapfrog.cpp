// Modification-Log
//
// mheinen_2013-08-12 --> FOCUS_SYSTEM_CENTER_OF_MASS

#include <map>

#include "Leapfrog.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "ensemble/PressureGradient.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "ensemble/EnsembleBase.h"
#include "Simulation.h"


using namespace std;
using Log::global_log;

Leapfrog::Leapfrog(double timestepLength) :	Integrator(timestepLength) {
	init();
}

void Leapfrog::init() {
	// set starting state
	_state = STATE_POST_FORCE_CALCULATION;
}

Leapfrog::~Leapfrog() {}

void Leapfrog::readXML(XMLfileUnits& xmlconfig) {
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	global_log->info() << "Timestep: " << _timestepLength << endl;
	assert(_timestepLength > 0);
}

void Leapfrog::eventForcesCalculated(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_PRE_FORCE_CALCULATION) {
		transition2to3(molCont, domain);
	}
}

void Leapfrog::eventNewTimestep(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		transition3to1(molCont, domain);
		transition1to2(molCont, domain);
	}
}

void Leapfrog::transition1to2(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_NEW_TIMESTEP) {
		Molecule* tempMolecule;
		double vcorr = 2. - 1. / domain->getGlobalBetaTrans();
		double Dcorr = 2. - 1. / domain->getGlobalBetaRot();

/*
		// begin --> mheinen_2013-08-16 --> FOCUS_SYSTEM_CENTER_OF_MASS
		double dPos[3] = {0.0, 0.0, 0.0};
		int nNumMols = molCont->getNumberOfParticles();
		global_log->info() << "number of particles: " << nNumMols << endl;
		// end <-- FOCUS_SYSTEM_CENTER_OF_MASS

		// begin --> mheinen_2013-08-12 --> FOCUS_SYSTEM_CENTER_OF_MASS
		double rSum[3];
		rSum[0] = 0.0;
		rSum[1] = 0.0;
		rSum[2] = 0.0;
		// end <-- FOCUS_SYSTEM_CENTER_OF_MASS
*/
		for (tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next())
		{
			tempMolecule->upd_preF(_timestepLength, vcorr, Dcorr);
/*
			// begin --> mheinen_2013-08-12 --> FOCUS_SYSTEM_CENTER_OF_MASS
			rSum[0] += tempMolecule->r(0);
			rSum[1] += tempMolecule->r(1);
			rSum[2] += tempMolecule->r(2);

			global_log->info() << "rSum x: " << rSum[0] << " " << "rSum y: " << rSum[1] << " " << "rSum z: " << rSum[2] << " " << endl;
			// end <-- FOCUS_SYSTEM_CENTER_OF_MASS
*/
		}
/*
		// begin --> mheinen_2013-08-12 --> FOCUS_SYSTEM_CENTER_OF_MASS
		domain->SetSystemCenterOfMass(0, rSum[0] / (double)nNumMols );
		domain->SetSystemCenterOfMass(1, rSum[1] / (double)nNumMols );
		domain->SetSystemCenterOfMass(2, rSum[2] / (double)nNumMols );
		// end <-- FOCUS_SYSTEM_CENTER_OF_MASS

		// begin --> mheinen_2013-08-16 --> FOCUS_SYSTEM_CENTER_OF_MASS

		// Systemschwerpunkt in die Mitte des Systems verschieben
		global_log->info() << "length x: " << domain->getGlobalLength(0) << " " << "length y: " << domain->getGlobalLength(1) << " " << "length z: " << domain->getGlobalLength(2) << " " << endl;

		for (tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next())
		{
			dPos[0] += domain->getGlobalLength(0) / 2.0 - domain->GetSystemCenterOfMass(0);
			dPos[1] += domain->getGlobalLength(1) / 2.0 - domain->GetSystemCenterOfMass(1);
			dPos[2] += domain->getGlobalLength(2) / 2.0 - domain->GetSystemCenterOfMass(2);
			tempMolecule->setr(0, dPos[0] );
			tempMolecule->setr(1, dPos[1] );
			tempMolecule->setr(2, dPos[2] );
		}
		// end <-- FOCUS_SYSTEM_CENTER_OF_MASS
*/

		this->_state = STATE_PRE_FORCE_CALCULATION;
	}
	else {
		global_log->error() << "Leapfrog::transition1to2(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::transition2to3(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_PRE_FORCE_CALCULATION) {
		Molecule* tM;
		map<int, unsigned long> N;
		map<int, unsigned long> rotDOF;
		map<int, double> summv2;
		map<int, double> sumIw2;
		double dt_half = 0.5 * this->_timestepLength;
		if (domain->severalThermostats()) {
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				int cid = tM->componentid();
				int thermostat = domain->getThermostat(cid);
				tM->upd_postF(dt_half, summv2[thermostat], sumIw2[thermostat]);
				N[thermostat]++;
				rotDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
			}
		}
		else {
			unsigned long Ngt = 0;
			unsigned long rotDOFgt = 0;
			double summv2gt = 0.0;
			double sumIw2gt = 0.0;
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				tM->upd_postF(dt_half, summv2gt, sumIw2gt);
				assert(summv2gt >= 0.0);
				Ngt++;
				rotDOFgt += tM->component()->getRotationalDegreesOfFreedom();
			}
			N[0] = Ngt;
			rotDOF[0] = rotDOFgt;
			summv2[0] = summv2gt;
			sumIw2[0] = sumIw2gt;
		}
		for (map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
			assert(thermit->second > 0);
			domain->setLocalSummv2(thermit->second, thermit->first);
			domain->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
			domain->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
		}

		this->_state = STATE_POST_FORCE_CALCULATION;
	}
	else {
		global_log->error() << "Leapfrog::transition2to3(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::transition3to1(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		this->_state = STATE_NEW_TIMESTEP;
	}
	else {
		global_log->error() << "Leapfrog::transition3to1(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::accelerateUniformly(ParticleContainer* molCont, Domain* domain) {
	map<unsigned, double>* additionalAcceleration = domain->getPG()->getUAA();
	vector<Component> comp = *(_simulation.getEnsemble()->components());
	vector<Component>::iterator compit;
	map<unsigned, double> componentwiseVelocityDelta[3];
	for (compit = comp.begin(); compit != comp.end(); compit++) {
		unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
		if (cosetid != 0)
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = _timestepLength * additionalAcceleration[d][cosetid];
		else
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = 0;
	}

	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}

void Leapfrog::accelerateInstantaneously(ParticleContainer* molCont, Domain* domain) {
	vector<Component> comp = *(_simulation.getEnsemble()->components());
	vector<Component>::iterator compit;
	map<unsigned, double> componentwiseVelocityDelta[3];
	for (compit = comp.begin(); compit != comp.end(); compit++) {
		unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
		if (cosetid != 0)
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = domain->getPG()->getMissingVelocity(cosetid, d);
		else
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = 0;
	}

	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}

void Leapfrog::init1D(unsigned zoscillator, ParticleContainer* molCont) {
	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
		if (!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->setXY();
}

void Leapfrog::zOscillation(unsigned zoscillator, ParticleContainer* molCont) {
	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
		/* TODO: use cid instead of complicated id + tersoff */
		if (!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->resetXY();
}
