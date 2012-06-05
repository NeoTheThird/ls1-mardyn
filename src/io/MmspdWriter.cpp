#include "io/MmspdWriter.h"

#include "Common.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#ifdef ENABLE_MPI
#include <mpi.h>
#include <mpi.h>
#endif

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
	ofstream mmspdfstream(_filename.c_str(), ios::binary|ios::out);
	cout << "Ofstream zum 1. mal erzeugt.\n";
  
  /* writing the header of the mmspd file, i.e. writing the BOM, the format marker (UTF-8),  the header line and defining the particle types */
  // BOM
  short int bom1,bom2,bom3;
  bom1 = 0xef;
  bom2 = 0xbb;
  bom3 = 0xbf;
  
  mmspdfstream.write(reinterpret_cast<const char*>(& bom1), 1);
  mmspdfstream.write(reinterpret_cast<const char*>(& bom2), 1);
  mmspdfstream.write(reinterpret_cast<const char*>(& bom3), 1);
  mmspdfstream.close();
  cout << "BOM geschrieben und stream geschlossen\n";
  // format marker
  mmspdfstream.open(_filename.c_str(), ios::out|ios::app);
  cout << " ofstream wieder geoeffnet zum neuen Schreiben d. Headers\n";
  mmspdfstream << "MMSPDu 1.0" << "\n";
  // header line
  mmspdfstream << "1 " << particleContainer->getBoundingBoxMin(0) << " " << particleContainer->getBoundingBoxMin(1) << " " 
		       << particleContainer->getBoundingBoxMin(2) << " " << particleContainer->getBoundingBoxMax(0) << " " 
		       << particleContainer->getBoundingBoxMax(1) << " " << particleContainer->getBoundingBoxMax(2) << " "
		       << _numberOfTimesteps / _writeFrequency+1    << " " << domain-> getNumberOfComponents() << " " << "0" << "\n";
		       
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
      mmspdfstream<< setprecision(4) << domain->getSigma(i,0) << " x f y f z f" << "\n";
  } // end of particle definitions		
  
  mmspdfstream.close();
} // end initOutput()

void MmspdWriter::doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* lmu){
	if (simstep % _writeFrequency == 0) {
#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	int tag = 4711;
	if (rank == 0){
#endif
		ofstream mmspdfstream(_filename.c_str(), ios::out|ios::app);
		mmspdfstream << "> " << domain->getglobalNumMolecules() << "\n";
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				mmspdfstream << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(3)
					<< pos->componentid() << setprecision(3) << " ";
				for (unsigned short d = 0; d < 3; d++) mmspdfstream << setw(7) << pos->r(d) << " " ;
				mmspdfstream << "\n";
			}
		}
#ifdef ENABLE_MPI
		for(int fromrank = 1; fromrank < domainDecomp->getNumProcs(); fromrank++) {
			MPI_Status status_probe;
			MPI_Status status_recv;
			MPI_Probe(fromrank, tag, MPI_COMM_WORLD, &status_probe);
			int numchars;
			MPI_Get_count(&status_probe, MPI_CHAR, &numchars);
			char *recvbuff = new char[numchars];
			MPI_Recv(recvbuff, numchars, MPI_CHAR, fromrank, tag, MPI_COMM_WORLD, &status_recv);
			mmspdfstream << string(recvbuff);
			delete recvbuff;
		}
#endif
		mmspdfstream.close();
#ifdef ENABLE_MPI
	}
	else {
		stringstream mmspdfstream;
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				mmspdfstream << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(3)
					<< pos->componentid() << setprecision(3) << " ";
				for (unsigned short d = 0; d < 3; d++) mmspdfstream << setw(7) << pos->r(d) << " " ;
				mmspdfstream << "\n";
			}
		}
		
		string sendbuff;
		sendbuff = mmspdfstream.str();
		MPI_Send((char*)sendbuff.c_str(), sendbuff.length() + 1, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
	}
#endif
  }
} // end doOutput

void MmspdWriter::finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain){
  
}
