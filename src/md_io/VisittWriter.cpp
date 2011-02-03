#include "md_io/VisittWriter.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "molecules/Molecule.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef GRANDCANONICAL
#include "ensemble/GrandCanonical.h"
#endif

md_io::VisittWriter::VisittWriter( unsigned long numberOfTimesteps,
                                   unsigned long writeFrequency, string outputPathAndPrefix )
{
   this->_writeFrequency = writeFrequency;
   this->_outputPathAndPrefix = outputPathAndPrefix;
   this->_wroteVIS = false;
}

md_io::VisittWriter::~VisittWriter(){}

void md_io::VisittWriter::initOutput( datastructures::ParticleContainer<Molecule>* pC,
                                      parallel::DomainDecompBase* domainDecomp,
                                      Domain* domain, unsigned coord )
{
   this->_wroteVIS = false;
   this->_coord = coord;
}

void md_io::VisittWriter::doOutput( datastructures::ParticleContainer<Molecule>* pC,
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
    filenamestream << ".vis.itt_";
    ofstream visittfstrm(filenamestream.str().c_str());

    if(!domain->ownrank() && !this->_wroteVIS)
    {
      visittfstrm << "      id t          x          y          z     q0     q1     q2     q3        c\n";
      this->_wroteVIS = true;
    }
    else if(!domain->ownrank()) visittfstrm << "#\n";

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
        visittfstrm << setiosflags(ios::fixed) << setw(8) << pCit->id() << setw(2)
                    << pCit->componentid() << setprecision(3);
        for(unsigned short d = 0; d < 3; d++) visittfstrm << setw(11) << pCit->r(d);
        visittfstrm << setprecision(3) << setw(7) << pCit->q().qw() << setw(7) << pCit->q().qx()
                    << setw(7) << pCit->q().qy()<< setw(7) << pCit->q().qz()
                    << setw(9) << right << (((_coord > 0) && (pCit->getCurTWFN() >= _coord))? 1: 0) << "\n";
      }
    }
    visittfstrm.close();
  }
}

void md_io::VisittWriter::finishOutput(datastructures::ParticleContainer<Molecule>* pC,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain)
{
}

