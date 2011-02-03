#ifndef XYZWRITER_H_
#define XYZWRITER_H_

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
  class XyzWriter; 
}
using namespace std;

class md_io::XyzWriter : public md_io::OutputBase{
 public:
  XyzWriter(unsigned long numberOfTimesteps, unsigned long writeFrequency, string outputPathAndPrefix);
  ~XyzWriter();
  void initOutput(
     datastructures::ParticleContainer<Molecule>* particleContainer,
     parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned coord
  );
  void doOutput(
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
  unsigned long _writeFrequency;
  string _outputPathAndPrefix;
};

#endif /*XYZWRITER_H_*/
