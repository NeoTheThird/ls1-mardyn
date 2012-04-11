#include "io/MmspdWriter.h"

#include "Common.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"

#include <fstream>
#include <sstream>
#include <iomanip>


using namespace std;

MmspdWriter::MmspdWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
	_filename = filename;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	_numberOfTimesteps = numberOfTimesteps;

	if (filename == "default")
		_filenameisdate = true;
	else
		_filenameisdate = false;
}

MmspdWriter::~MmspdWriter(){}

void MmspdWriter::initOutput(ParticleContainer* particleContainer,
			   DomainDecompBase* domainDecomp, Domain* domain){
	
  
	if (_filenameisdate) {
			_filename = _filename + "mardyn";
			_filename = _filename + gettimestring();
		} 
	_filename = _filename +  ".mmspd";
		/*
	// generating the *.mmspd file 
	stringstream filenamestream;
		if (_filenameisdate) {
			filenamestream << "mardyn" << gettimestring();
		}
		else {
			filenamestream << _filename;
		}

		if (_incremental) {
			// align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. 
			int num_digits = (int) ceil(log(double(_numberOfTimesteps / _writeFrequency)) / log(10.));
			filenamestream << aligned_number(simstep / _writeFrequency, num_digits, '0');
		}
		filenamestream << ".mmspd";*/
		
	ofstream mmspdfstream(_filename.c_str(), ios::binary|ios::out);
  
  /* writing the header of the mmspd file, i.e. writing the format marker (UTF-8),  the header line and defining the particle types */
  // format marker
  mmspdfstream << "MMSPDu 1.0" << endl;
  // header line
  mmspdfstream << "1 " << particleContainer->getBoundingBoxMin(0) << " " << particleContainer->getBoundingBoxMin(1) << " " 
		       << particleContainer->getBoundingBoxMin(2) << " " << particleContainer->getBoundingBoxMax(0) << " " 
		       << particleContainer->getBoundingBoxMax(1) << " " << particleContainer->getBoundingBoxMax(2) << " "
		       << _numberOfTimesteps / _writeFrequency    << " " << domain-> getNumberOfComponents() << " " << "0" << endl;
		       
  // particle definitions every single line specifies a particular particle type
  for(unsigned i = 0; i < domain->getNumberOfComponents() ; i++){
      if (i == 0){
	mmspdfstream << "s 4 3 cr b 255 cg b 255 cb b 0 r f "; 
      }
      else if (i == 1){
	mmspdfstream << "s 4 3 cr b 50 cg b 50 cb b 50 r f "; 
      }
      else if (i == 2){
	mmspdfstream << "s 4 3 cr b 0 cg b 255 cb b 255 r f "; 
      }
      else if(i == 3){
	mmspdfstream << "s 4 3 cr b 150 cg b 0 cb b 150 r f "; 
      }
      else if (i == 4){
	mmspdfstream << "s 4 3 cr b 100 cg b 100 cb b 100 r f "; 
      }
      else {
	mmspdfstream << "**************** Error: Unspecified component!*************\n Possible reason: more than 5 components?\n"; 
      }
      mmspdfstream<< setprecision(4) << domain->getSigma(i,0) << " x f x f x f" << endl;
  } // end of particle definitions		
  
  mmspdfstream.close();
} // end initOutput()

void MmspdWriter::doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* lmu){
  
	
	ofstream mmspdfstream(_filename.c_str(), ios::out|ios::app);
	int rank = domainDecomp->getRank();
	if (rank == 0) mmspdfstream << "> " << domain->getglobalNumMolecules() << endl;
	
	for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				mmspdfstream << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(2)
				            << pos->componentid() << setprecision(3);
				for (unsigned short d = 0; d < 3; d++) mmspdfstream << setw(7) << pos->r(d);
			}
		}
	
	mmspdfstream.close();
} // end doOutput

void MmspdWriter::finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain){
  
}