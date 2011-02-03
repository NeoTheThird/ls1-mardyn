#include "md_io/VimWriter.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "molecules/Molecule.h"

#include <iostream>
#include <iomanip>
#include <sstream>

md_io::VimWriter::VimWriter( unsigned long numberOfTimesteps,
                                   unsigned long writeFrequency, string outputPathAndPrefix )
{
   this->_writeFrequency = writeFrequency;
   this->_outputPathAndPrefix = outputPathAndPrefix;
   this->_wroteVIM = false;
}

md_io::VimWriter::~VimWriter(){}

void md_io::VimWriter::initOutput( datastructures::ParticleContainer<Molecule>* pC,
                                   parallel::DomainDecompBase* domainDecomp,
                                   Domain* domain, unsigned coord )
{
   this->_wroteVIM = false;
   this->_coord = coord;
}

void md_io::VimWriter::doOutput( datastructures::ParticleContainer<Molecule>* pC,
                                    parallel::DomainDecompBase* domainDecomp,
                                    Domain* domain, unsigned long simstep
#ifdef GRANDCANONICAL
                                    ,
                                    list<ensemble::ChemicalPotential>* lmu
#endif
)
{
  if(simstep%_writeFrequency == 0)
  {
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
    filenamestream << ".vim_";
    ofstream vimfstrm(filenamestream.str().c_str());

    double lmax = 0.0;
    for(int d=0; d < 3; d++) if(domain->getGlobalLength(d) > lmax) lmax = domain->getGlobalLength(d);

    unsigned Ncomp = domain->getComponents().size();
    if(!domain->ownrank() && !this->_wroteVIM)
    {
      vector<Component> tcomp = domain->getComponents();
      for(vector<Component>::iterator tcit = tcomp.begin(); tcit != tcomp.end(); tcit++)
      {
        tcit->writeVIM(vimfstrm, 0);
      }
      if(_coord > 0)
      {
         for(vector<Component>::iterator tcit = tcomp.begin(); tcit != tcomp.end(); tcit++)
         {
           tcit->writeVIM(vimfstrm, Ncomp);
         }
      }
      this->_wroteVIM = true;
    }
    if(domain->ownrank() == 0)
      vimfstrm << "# " << setprecision(3) << setw(10) << lmax << "\n";

    for(Molecule* pCit = pC->begin(); pCit != pC->end(); pCit = pC->next())
    {
      /*
      visittfstrm << tempMol->componentid() << " ";
      visittfstrm << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << endl;
      */
      bool halo = false;
      for(unsigned short d = 0; d < 3; d++)
      {
         if((pCit->r(d) < pC->getBoundingBoxMin(d)) || (pCit->r(d) > pC->getBoundingBoxMax(d)))
         {
           halo = true;
           break;
         }
      }
      if(!halo)
      {
        vimfstrm << "! " << setprecision(3)
                 << setw(2) << pCit->componentid()+1 + (((_coord > 0) && (pCit->getCurTWFN() >= _coord))? Ncomp: 0)
                 << ' ';
        for(unsigned short d = 0; d < 3; d++)
          vimfstrm << setw(3) << floor(999.9999 * pCit->r(d) / lmax) << ' ';
        vimfstrm << setw(4) << floor(999.9999 * pCit->q().qw()) << ' '
                 << setw(4) << floor(999.9999 * pCit->q().qx()) << ' '
                 << setw(4) << floor(999.9999 * pCit->q().qy()) << ' '
                 << setw(4) << floor(999.9999 * pCit->q().qz()) << '\n';
      }
    }
    vimfstrm.close();
  }
}

void md_io::VimWriter::finishOutput(datastructures::ParticleContainer<Molecule>* pC,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain)
{
}
