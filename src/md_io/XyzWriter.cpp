
#include "md_io/XyzWriter.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "molecules/Molecule.h"

md_io::XyzWriter::XyzWriter(unsigned long numberOfTimesteps, unsigned long writeFrequency, string outputPathAndPrefix){
  _writeFrequency = writeFrequency;
  _outputPathAndPrefix = outputPathAndPrefix;
}

md_io::XyzWriter::~XyzWriter(){}

void md_io::XyzWriter::initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned coord){
}

void md_io::XyzWriter::doOutput(
   datastructures::ParticleContainer<Molecule>* particleContainer,
   parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep
#ifdef GRANDCANONICAL
     ,
     list<ensemble::ChemicalPotential>* lmu
#endif
) {
  if(simstep%_writeFrequency == 0) {
    stringstream filenamestream;
    filenamestream << _outputPathAndPrefix;
    // unsigned long temp = simstep/_writeFrequency;
    // while(temp < floor((long double)_numberOfTimesteps/_writeFrequency)){
    //   filenamestream << "0";
    //   temp = temp*10;
    // }
    filenamestream.width(6);
    filenamestream.fill('0');
    filenamestream << right << (unsigned)(simstep / _writeFrequency);
    filenamestream << 'r';
    filenamestream.width(3);
    filenamestream.fill('0');
    filenamestream << right << domain->ownrank();
    filenamestream << ".buxyz";

    ofstream xyzfilestream(filenamestream.str().c_str());
    if(!domain->ownrank())
    {
       xyzfilestream << domain->N() << endl;
       xyzfilestream << "comment line" << endl;
    }
    Molecule* tempMol;
    for(tempMol = particleContainer->begin(); tempMol != particleContainer->end(); tempMol = particleContainer->next())
    {
      xyzfilestream << tempMol->componentid() << " ";
      xyzfilestream << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << endl;
    }
    xyzfilestream.close();
  }
}

void md_io::XyzWriter::finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
}
