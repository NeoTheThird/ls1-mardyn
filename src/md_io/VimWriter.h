#ifndef VIMWRITER_H_
#define VIMWRITER_H_

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
  class VimWriter; 
}
using namespace std;

class md_io::VimWriter : public md_io::OutputBase{
 public:
  VimWriter(unsigned long numberOfTimesteps, unsigned long writeFrequency, string outputPathAndPrefix);
  ~VimWriter();
  void initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
  void doOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
  void finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
 private:
  unsigned long _numberOfTimesteps; 
  unsigned _writeFrequency;
  string _outputPathAndPrefix;
  bool _wroteVIM;
};

#endif /*VISITTWRITER_H_*/
