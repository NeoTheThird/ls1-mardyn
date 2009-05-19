#ifndef VISITTWRITER_H_
#define VISITTWRITER_H_

#include "md_io/OutputBase.h"
#include <string>
#include <fstream>

namespace datastructures{
  template<class ParticleType>
  class ParticleContainer;
}

namespace parallel{
  class DomainDecompBase; 
}

class Domain;
class Molecule;

namespace md_io{
  class VisittWriter; 
}
using namespace std;

class md_io::VisittWriter : public md_io::OutputBase{
 public:
  VisittWriter(unsigned long numberOfTimesteps, unsigned long writeFrequency, string outputPathAndPrefix);
  ~VisittWriter();
  void initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
  void doOutput
  (
     datastructures::ParticleContainer<Molecule>* particleContainer,
     parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep
#ifdef GRANDCANONICAL
     ,
     list<ensemble::ChemicalPotential>* lmu
#endif
  );
  void finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
 private:
  unsigned long _numberOfTimesteps; 
  unsigned _writeFrequency;
  string _outputPathAndPrefix;
  bool _wroteVIS;
};

#endif /*VISITTWRITER_H_*/
