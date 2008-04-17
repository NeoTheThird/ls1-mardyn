
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
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
}

void md_io::XyzWriter::doOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep){
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
    filenamestream << right << (simstep / _writeFrequency);
    filenamestream << 'r';
    filenamestream.width(3);
    filenamestream.fill('0');
    filenamestream << right << domain->ownrank();
    filenamestream << simstep/_writeFrequency << ".buxyz";

    ofstream xyzfilestream(filenamestream.str().c_str());
    xyzfilestream << particleContainer->getNumberOfParticles() << endl;
    xyzfilestream << "comment line" << endl;
    Molecule* tempMol;
    for(tempMol = particleContainer->begin(); tempMol != particleContainer->end(); tempMol = particleContainer->next())
    {
      // if(tempMol->componentid() == 0) { xyzfilestream << "Ar ";}
      // else if(tempMol->componentid() == 1) {xyzfilestream << "Xe ";}
      // else if(tempMol->componentid() == 2) {xyzfilestream << "C ";}
      // else if(tempMol->componentid() == 3) {xyzfilestream << "O ";}
      // else {xyzfilestream << "H ";}
      xyzfilestream << tempMol->componentid() << " ";
      xyzfilestream << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << endl;
    }
    xyzfilestream.close();
  }
}

void md_io::XyzWriter::finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain){
}
