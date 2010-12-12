#include "integrators/Leapfrog.h"
#include "datastructures/ParticleContainer.h"
#include "Domain.h"
#include "molecules/Molecule.h"

utils::Log integrators::Leapfrog::_log("Leapfrog");

integrators::Leapfrog::Leapfrog(double timestepLength){
  // set starting state
  this->_state = 3; 
  
  this->_timestepLength = timestepLength;
}

integrators::Leapfrog::~Leapfrog(){
}


void integrators::Leapfrog::eventForcesCalculated(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state == 2){
    transition2to3(molCont, domain);
  }
}

void integrators::Leapfrog::eventNewTimestep(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state == 3){
    transition3to1(molCont, domain);
    transition1to2(molCont, domain);
  }
}


void integrators::Leapfrog::transition1to2(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state==1){
    Molecule* tempMolecule;
    double vcorr=2.-1./domain->getGlobalBetaTrans();
    double Dcorr=2.-1./domain->getGlobalBetaRot();
    for(tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()){
#ifndef NDEBUG
      // if(!(tempMolecule->id() % 43690))
      // {
         // cout << "molecule " << tempMolecule->id() << " in leapfrog state 1:\n"
         //      << "r\t" << tempMolecule->r(0) << "\t" << tempMolecule->r(1) << "\t" << tempMolecule->r(2)
         //      << "\nv\t" << tempMolecule->v(0) << "\t" << tempMolecule->v(1) << "\t" << tempMolecule->v(2)
         //      << "\nF\t" << tempMolecule->F(0) << "\t" << tempMolecule->F(1) << "\t" << tempMolecule->F(2) << "\n";
         // cout << "vcorr: " << vcorr << "\tDcorr: " << Dcorr << "\n";
      // }
#endif
      tempMolecule->upd_preF(_timestepLength, vcorr, Dcorr);
#ifndef NDEBUG
      // if(!(tempMolecule->id() % 43690))
      // {
         // cout << "molecule " << tempMolecule->id() << " in leapfrog state 2:\n"
         //      << "r\t" << tempMolecule->r(0) << "\t" << tempMolecule->r(1) << "\t" << tempMolecule->r(2)
         //      << "\nv\t" << tempMolecule->v(0) << "\t" << tempMolecule->v(1) << "\t" << tempMolecule->v(2)
         //      << "\nF\t" << tempMolecule->F(0) << "\t" << tempMolecule->F(1) << "\t" << tempMolecule->F(2) << "\n";
      // }
#endif
    }

    this->_state = 2;
  }
  else { _log.error("transition1to2(...)","Wrong state for state transition");}
}


void integrators::Leapfrog::transition2to3(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state==2){
    /*
    Molecule* tempMolecule;
    double summv2 = 0.0;
    double sumIw2 = 0.0;
    double dt_halve=.5*_timestepLength;
    
    for(tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()){
      tM->upd_postF(dt_halve,summv2,sumIw2); 
    }
    
    domain->setLocalSummv2(summv2);
    domain->setLocalSumIw2(sumIw2);
    */

    Molecule* tM;
#ifdef COMPLEX_POTENTIAL_SET
    map<int, unsigned long> N;
    map<int, unsigned long> rotDOF;
#endif
    map<int, double> summv2;
    map<int, double> sumIw2;
    double dt_half = 0.5 * this->_timestepLength;
    if(domain->severalThermostats())
    {
      for(tM = molCont->begin(); tM != molCont->end(); tM = molCont->next())
      {
        int cid = tM->componentid();
        int thermostat = domain->getThermostat(cid);
        tM->upd_postF(dt_half, summv2[thermostat], sumIw2[thermostat]); 
#ifdef COMPLEX_POTENTIAL_SET
        N[thermostat]++;
        rotDOF[thermostat] += domain->getComponentRotDOF(cid);
#endif
      }
    }
    else
    {
#ifdef COMPLEX_POTENTIAL_SET
      unsigned long Ngt = 0;
      unsigned long rotDOFgt = 0;
#endif
      double summv2gt = 0.0;
      double sumIw2gt = 0.0;
      for(tM = molCont->begin(); tM != molCont->end(); tM = molCont->next())
      {
#ifdef COMPLEX_POTENTIAL_SET
         int cid = tM->componentid();
#endif
         tM->upd_postF(dt_half, summv2gt, sumIw2gt); 
         assert(summv2gt >= 0.0);
#ifdef COMPLEX_POTENTIAL_SET
         Ngt++;
         rotDOFgt += domain->getComponentRotDOF(cid);
#endif
      }
#ifdef COMPLEX_POTENTIAL_SET
      N[0] = Ngt;
      rotDOF[0] = rotDOFgt;
#endif
      summv2[0] = summv2gt;
      sumIw2[0] = sumIw2gt;
    }
    for(map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++)
    {
       assert(thermit->second > 0);
       domain->setLocalSummv2(thermit->second, thermit->first);
       domain->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
#ifdef COMPLEX_POTENTIAL_SET
       domain->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
#endif
    }
    this->_state = 3;
  }
  else { _log.error("transition2to3(...)","Wrong state for state transition");}
}


void integrators::Leapfrog::transition3to1(datastructures::ParticleContainer<Molecule>* molCont, Domain* domain){
  if(this->_state==3){ this->_state = 1;}
  else{ _log.error("transition3to1(...)","Wrong state for state transition");}
}

#ifdef COMPLEX_POTENTIAL_SET
/*
 * diese Version beschleunigt nur in z-Richtung
 *
void integrators::Leapfrog::accelerateUniformly
(
   datastructures::ParticleContainer<Molecule>* molCont,
   Domain* domain   )
{
   map<unsigned, double>* additionalAcceleration = domain->getUAA();
   vector<Component> comp = domain->getComponents();
   vector<Component>::iterator compit;
   map<unsigned, double> componentwiseVelocityDelta;
   for(compit = comp.begin(); compit != comp.end(); compit++)
   {
      unsigned cosetid = domain->getComponentSet(compit->ID());
      if(cosetid != 0)
            componentwiseVelocityDelta[compit->ID()]
               = _timestepLength * additionalAcceleration[2][cosetid];
      else
            componentwiseVelocityDelta[compit->ID()] = 0;
   }

   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
   {
      unsigned cid = thismol->componentid();
#ifndef NDEBUG
      assert(componentwiseVelocityDelta.find(cid) != componentwiseVelocityDelta.end());
#endif
      thismol->vadd(0.0, 0.0, componentwiseVelocityDelta[cid]);
   }
}
 *
 */

/*
 * diese Version beschleunigt in alle Raumrichtungen
 */
void integrators::Leapfrog::accelerateUniformly
(
   datastructures::ParticleContainer<Molecule>* molCont,
   Domain* domain   )
{
   map<unsigned, double>* additionalAcceleration = domain->getUAA();
   vector<Component> comp = domain->getComponents();
   vector<Component>::iterator compit;
   map<unsigned, double> componentwiseVelocityDelta[3];
   for(compit = comp.begin(); compit != comp.end(); compit++)
   {
      unsigned cosetid = domain->getComponentSet(compit->ID());
      if(cosetid != 0)
         for(unsigned d = 0; d < 3; d++)
            componentwiseVelocityDelta[d][compit->ID()]
               = _timestepLength * additionalAcceleration[d][cosetid];
      else
         for(unsigned d = 0; d < 3; d++)
            componentwiseVelocityDelta[d][compit->ID()] = 0;
#ifndef NDEBUG
      // cout << "acc. for comp. [" << compit->ID() << "]\t";
      // for(unsigned d = 0; d < 3; d++) cout << componentwiseVelocityDelta[d][compit->ID()] << " ";
      // cout << "coset " << cosetid << "\n";
#endif
   }

   Molecule* thismol;
   if(domain->aladin()) for(
      thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()
   ) {
      unsigned cid = thismol->componentid();
      if(domain->doAccelerate(thismol->r(0), thismol->r(1), thismol->r(2)))
         thismol->vadd( componentwiseVelocityDelta[0][cid],
                        componentwiseVelocityDelta[1][cid], componentwiseVelocityDelta[2][cid] );
   }
   else for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
   {
      unsigned cid = thismol->componentid();
      // assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
      thismol->vadd( componentwiseVelocityDelta[0][cid],
                     componentwiseVelocityDelta[1][cid], componentwiseVelocityDelta[2][cid] );
   }
}

/*
 * diese Version beschleunigt nur in z-Richtung
 *
void integrators::Leapfrog::accelerateInstantaneously
(
   datastructures::ParticleContainer<Molecule>* molCont,
   Domain* domain   )
{
   vector<Component> comp = domain->getComponents();
   vector<Component>::iterator compit;
   map<unsigned, double> componentwiseVelocityDelta;
   for(compit = comp.begin(); compit != comp.end(); compit++)
   {
      unsigned cosetid = domain->getComponentSet(compit->ID());
      if(cosetid != 0)
            componentwiseVelocityDelta[compit->ID()]
               = domain->getMissingVelocity(cosetid, 2);
      else
            componentwiseVelocityDelta[compit->ID()] = 0;
   }

   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
   {
      unsigned cid = thismol->componentid();
      assert(componentwiseVelocityDelta.find(cid) != componentwiseVelocityDelta.end());
      thismol->vadd(0.0, 0.0, componentwiseVelocityDelta[cid]);
   }
}
 *
 */

/*
 * diese Version beschleunigt in alle Raumrichtungen
 */
void integrators::Leapfrog::accelerateInstantaneously
(
   datastructures::ParticleContainer<Molecule>* molCont,
   Domain* domain   )
{
   vector<Component> comp = domain->getComponents();
   vector<Component>::iterator compit;
   map<unsigned, double> componentwiseVelocityDelta[3];
   for(compit = comp.begin(); compit != comp.end(); compit++)
   {
      unsigned cosetid = domain->getComponentSet(compit->ID());
      if(cosetid != 0)
         for(unsigned d = 0; d < 3; d++)
            componentwiseVelocityDelta[d][compit->ID()]
               = domain->getMissingVelocity(cosetid, d);
      else
         for(unsigned d = 0; d < 3; d++)
            componentwiseVelocityDelta[d][compit->ID()] = 0;
   }

   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
   {
      unsigned cid = thismol->componentid();
      assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
      thismol->vadd( componentwiseVelocityDelta[0][cid],
                     componentwiseVelocityDelta[1][cid], componentwiseVelocityDelta[2][cid] );
   }
}

void integrators::Leapfrog::init1D(
   unsigned zoscillator,
   datastructures::ParticleContainer<Molecule>* molCont
)
{
   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
      if(!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->setXY();
}

void integrators::Leapfrog::init0D(
   unsigned oscillator,
   datastructures::ParticleContainer<Molecule>* molCont
)
{
   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
      if(!(thismol->id() % oscillator) && thismol->numTersoff()) thismol->setXYZ();
}

void integrators::Leapfrog::zOscillation(
   unsigned zoscillator,
   datastructures::ParticleContainer<Molecule>* molCont
      )
{
   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
      if(!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->resetXY();
}

void integrators::Leapfrog::oscillation(
   unsigned oscillator,
   datastructures::ParticleContainer<Molecule>* molCont
      )
{
   Molecule* thismol;
   for(thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
      if(!(thismol->id() % oscillator) && thismol->numTersoff()) thismol->resetXYZ();
}
#endif

