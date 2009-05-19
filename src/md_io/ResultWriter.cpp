#include "md_io/ResultWriter.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

#ifdef GRANDCANONICAL
#include "ensemble/GrandCanonical.h"
#endif

md_io::ResultWriter::ResultWriter(unsigned wF, string outputPrefix)
{
   this->_writeFrequency = wF;
   this->_outputPrefix = outputPrefix;
}

md_io::ResultWriter::~ResultWriter(){}

void md_io::ResultWriter::initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
   
  // initialize result file
  string resultfile(_outputPrefix+".res");
  time_t now;
  time(&now);
  if(domainDecomp->getRank()==0){
    _resultStream.open(resultfile.c_str());
    _resultStream << "# LS1 _M_olekul_ARDYN_amik: simulation starting at " << ctime(&now) << endl;
    _resultStream << "#\tt\tU_pot\tp\t\t";
#ifdef GRANDCANONICAL
    _resultStream << "N(uVT)[1..i]\trho(uVT)[1..i]\t\t";
#endif
    for(unsigned i=0; domain->maxThermostat() >= i; i++)
      _resultStream << "T(" << i << ")\tbetaTrans(" << i << ")\tbetaRot(" << i << ")\t";
#ifdef COMPLEX_POTENTIAL_SET
    for(unsigned i=1; domain->maxCoset() >= i; i++)
      _resultStream << "\tN(" << i << ")\tvx vy vz\tv(" << i << ")\tax ay az\ta(" << i << ")\t";
#endif
    _resultStream << "\n";
    _resultStream.precision(6);
  }
}

void md_io::ResultWriter::doOutput
(
   datastructures::ParticleContainer<Molecule>* particleContainer,
   parallel::DomainDecompBase* domainDecomp,
   Domain* domain, unsigned long simstep
#ifdef GRANDCANONICAL
   ,
   list<ensemble::ChemicalPotential>* lmu
#endif
)
{
  if((simstep%_writeFrequency) || (domain->thermostatWarning())) return;

  if(domainDecomp->getRank()==0)
  {
    _resultStream << simstep << "\t" << domain->getCurrentTime()
                  << "\t" << domain->getAverageGlobalUpot() << "\t"
                  << domain->getGlobalPressure() << "\t\t";
#ifdef GRANDCANONICAL
    std::list< ensemble::ChemicalPotential >::iterator cpit;
    for(cpit = lmu->begin(); cpit != lmu->end(); cpit++)
    {
       _resultStream << cpit->getGlobalN() << "\t";
    }
    for(cpit = lmu->begin(); cpit != lmu->end(); cpit++)
    {
       _resultStream << cpit->getGlobalRho() << "\t";
    }
    _resultStream << "\t";
#endif
    for(unsigned i=0; domain->maxThermostat() >= i; i++)
      _resultStream << domain->T(i) << "\t" << domain->getGlobalBetaTrans(i)
                    << "\t" << domain->getGlobalBetaRot(i) << "\t";
#ifdef COMPLEX_POTENTIAL_SET
    for(unsigned i=1; domain->maxCoset() >= i; i++)
    {
      double a = domain->getDirectedVelocity(i);
      if(a > 0)
      {
        _resultStream << "\t" << domain->getCosetN(i)
                      << "\t" << domain->getDirectedVelocity(i, 0) << " "
                      << domain->getDirectedVelocity(i, 1) << " "
                      << domain->getDirectedVelocity(i, 2)
                      << "\t" << a
                      << "\t" << domain->getUniformAcceleration(i, 0) << " "
                      << domain->getUniformAcceleration(i, 1) << " "
                      << domain->getUniformAcceleration(i, 2) << " "
                      << "\t" << domain->getUniformAcceleration(i) << "\t";
      }
    }
#endif
    _resultStream << "\n";
    _resultStream.flush();
  }
}

void md_io::ResultWriter::finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
  _resultStream.close();
}
