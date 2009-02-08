#include "Simulation.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "datastructures/LinkedCells.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"
#ifdef PARALLEL
#include "parallel/DomainDecomposition.h"
#endif
#include "datastructures/adapter/ParticlePairs2PotForceAdapter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "md_io/OutputBase.h"
#include "md_io/ResultWriter.h"
#include "md_io/XyzWriter.h"
#include "md_io/VisittWriter.h"
#include "md_io/VimWriter.h"

#include <fstream>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

utils::Log Simulation::_log("Simulation");



Simulation::Simulation(int argc, char **argv){
#ifdef PARALLEL
    _domainDecomposition = (parallel::DomainDecompBase*) new parallel::DomainDecomposition(&argc, &argv);
#else
    _domainDecomposition = (parallel::DomainDecompBase*) new parallel::DomainDecompDummy();
#endif
  int ownrank = _domainDecomposition->getRank();

  // two parameters are necessary: config-file and number of timesteps
  // if any other number of parameters is provied, print usage message and exit
  if((argc != 4) && (argc != 5)) {
    if(ownrank == 0) {
      cout << "Usage: " << argv[0] << " <configfilename> <number of timesteps> <outputprefix> [-Y2]" << endl;
      cout << "Detected argc: " << argc << "\n";
    }
    delete _domainDecomposition;
    exit(1);
  }
  
  // open filestream to the input file
  string inputfilename(argv[1]);  
  ifstream inputfilestream(inputfilename.c_str());
  
  // store number of timesteps to be simulated
  _numberOfTimesteps = atol(argv[2]);

  // store prefix for output files
  _outputPrefix = string(argv[3]);

  // used to store one token of the inputfilestream
  string token;
  string phaseSpaceFileName;
  double timestepLength;
#ifdef COMPLEX_POTENTIAL_SET
  unsigned cosetid = 0;
#endif

  this->_resultOutputTimesteps = 25;
#ifdef COMPLEX_POTENTIAL_SET
  this->_doRecordProfile = false;
  this->_profileRecordingTimesteps = 7;
  this->_profileOutputTimesteps = 12500;
  this->_profileOutputPrefix = "out";
  this->_collectThermostatDirectedVelocity = 100;
  this->_zoscillation = false;
  this->_zoscillator = 512;
#endif 
  this->_doRecordRDF = false;
  this->_RDFOutputTimesteps = 25000;
  this->_RDFOutputPrefix = "out";
  
#ifndef NDEBUG
  if(!ownrank) cout << "constructing domain:";
#endif
  this->_domain = new Domain(ownrank);
#ifndef NDEBUG
  if(!ownrank) cout << " done.";
#endif
  _particlePairsHandler = new datastructures::ParticlePairs2PotForceAdapter(*_domain);
  
  // The first line of the config file has to contain the token "MDProjectConfig"
  inputfilestream >> token;
  if((token != "MDProjectConfig") && (token != "mardynconfig")) {
    if(ownrank == 0) cerr << "Not a MDProject config file! " << token << endl;
    exit(1);
  } 
  if((ownrank == 0) && (token != "mardynconfig"))
  {
    cout << "Warning: initial token of plain input file should be 'mardynconfig'"
         << " instead of '" << token << "'.\n";
  }
  while(inputfilestream){
    token.clear();
    inputfilestream >> token;
    if(!ownrank) cout << " " << token << " \t";

    if(token.substr(0,1)=="#"){
      inputfilestream.ignore(1024,'\n');
    }
    else if(token=="phaseSpaceFile"){      
      inputfilestream >> phaseSpaceFileName;
      _domain->setPhaseSpaceFile(phaseSpaceFileName);
      if(timestepLength == 0.0)
      {
         cout << "timestep missing.\n";
         exit(1);
      }
      _domain->readPhaseSpaceHeader(
         timestepLength
#ifdef TRUNCATED_SHIFTED
            ,
         this->_cutoffRadius
#endif
#ifdef COMPLEX_POTENTIAL_SET
            ,
         this->_tersoffCutoffRadius
#endif
      );
      this->_numberOfComponents = _domain->getComponents().size();

      _domain->initParameterStreams(_cutoffRadius);
    }
    else if(token=="timestepLength") {
      inputfilestream >> timestepLength;
    }
    else if(token=="cutoffRadius")   {
      inputfilestream >> _cutoffRadius;
    } 
#ifdef COMPLEX_POTENTIAL_SET
    else if(token=="tersoffCutoffRadius")   {
      inputfilestream >> _tersoffCutoffRadius;
    } 
#endif
    else if(token=="datastructure")  {
      inputfilestream >> token;
      if(token=="LinkedCells"){

        int cellsInCutoffRadius;
        inputfilestream >> cellsInCutoffRadius;
        double bBoxMin[3];
        double bBoxMax[3];
        for(int i=0;i<3;i++) {
          bBoxMin[i] = _domainDecomposition->getCoords(i)*_domain->getGlobalLength(i)/_domainDecomposition->getGridSize(i);
          bBoxMax[i] = (_domainDecomposition->getCoords(i)+1)*_domain->getGlobalLength(i)/_domainDecomposition->getGridSize(i);
        }
        
        _moleculeContainer = new datastructures::LinkedCells<Molecule>(bBoxMin, bBoxMax, _cutoffRadius, 
#ifdef COMPLEX_POTENTIAL_SET
           _tersoffCutoffRadius,
#endif
           cellsInCutoffRadius, *_particlePairsHandler);
        
      }
    }
    else if(token=="output")  {
      inputfilestream >> token;
      if(token=="ResultWriter"){
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new md_io::ResultWriter(writeFrequency, outputPathAndPrefix));
      }
      else if(token=="XyzWriter"){
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        this->_restartOutputInterval = writeFrequency;
        // _outputPlugins.push_back(new md_io::XyzWriter(_numberOfTimesteps, writeFrequency, outputPathAndPrefix));
      }
      else if(token == "VisittWriter")
      {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new md_io::VisittWriter(_numberOfTimesteps, writeFrequency, outputPathAndPrefix));
        if(!ownrank) cout << "VisItt " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
      }
      else if(token == "VimWriter")
      {
        unsigned long writeFrequency;
        string outputPathAndPrefix;
        inputfilestream >> writeFrequency >> outputPathAndPrefix;
        _outputPlugins.push_back(new md_io::VimWriter(_numberOfTimesteps, writeFrequency, outputPathAndPrefix));
        if(!ownrank) cout << "Vim " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
      }
    }
#ifdef COMPLEX_POTENTIAL_SET
    else if(token == "accelerate")
    {
      cosetid++;
      inputfilestream >> token;
#ifndef NDEBUG
      if(!ownrank)
      {
         cout << token << "\n";
         cout << flush;
      } 
#endif
      if(token != "comp")
      {
         if(ownrank == 0) cout << "Expected 'comp' instead of '" << token << "'.\n";
         exit(1);
      }
      int cid = 0;
      while(cid >= 0)
      {
         inputfilestream >> cid;
         if(!ownrank && (cid > 0)) cout << "acc. for component " << cid << "\n";
         cid--;
         _domain->assignCoset((unsigned)cid, cosetid);
      }
      double v;
      inputfilestream >> v;
      if(!ownrank) cout << "velocity " << v << "\n";
      inputfilestream >> token;
#ifndef NDEBUG
      if(!ownrank)
      {
         cout << token << "\n";
         cout << flush;
      } 
#endif
      if(token != "towards")
      {
         if(ownrank == 0) cout << "Expected 'towards' instead of '" << token << "'.\n";
         exit(1);
      }
      double dir[3];
      double dirnorm = 0;
      for(unsigned d = 0; d < 3; d++)
      {
         inputfilestream >> dir[d];
         dirnorm += dir[d]*dir[d];
      }
      dirnorm = 1.0 / sqrt(dirnorm);
      for(unsigned d = 0; d < 3; d++) dir[d] *= dirnorm;
      inputfilestream >> token;
#ifndef NDEBUG
      if(!ownrank)
      {
         cout << token << "\n";
         cout << flush;
      } 
#endif
      if(token != "within")
      {
         if(ownrank == 0) cout << "Expected 'within' instead of '" << token << "'.\n";
         exit(1);
      }
      double tau;
      inputfilestream >> tau;
      inputfilestream >> token;
#ifndef NDEBUG
      if(!ownrank)
      {
         cout << token << "\n";
         cout << flush;
      } 
#endif
      if(token != "from")
      {
         if(ownrank == 0) cout << "Expected 'from' instead of '" << token << "'.\n";
         exit(1);
      }
      double ainit[3];
      for(unsigned d = 0; d < 3; d++) inputfilestream >> ainit[3];
      if(timestepLength == 0.0)
      {
         cout << "timestep missing.\n";
         exit(1);
      }
      this->_domain->specifyComponentSet(cosetid, dir, tau, ainit, timestepLength);
    }
    else if(token == "constantAccelerationTimesteps")
    {
       unsigned uCAT;
       inputfilestream >> uCAT;
       this->_domain->setUCAT(uCAT);
    }
    else if(token == "profile")
    {
       unsigned xun, yun, zun;
       inputfilestream >> xun >> yun >> zun;
       this->_domain->setupProfile(xun, yun, zun);
       this->_doRecordProfile = true;
    }
    else if(token == "profileRecordingTimesteps")
    {
       inputfilestream >> this->_profileRecordingTimesteps;
    }
    else if(token == "profileOutputTimesteps")
    {
       inputfilestream >> this->_profileOutputTimesteps;
    }
#endif
    else if(token == "RDF")
    {
       double interval;
       unsigned bins;
       inputfilestream >> interval >> bins;
       this->_domain->setupRDF(interval, bins);
       this->_doRecordRDF = true;
    }
    else if(token == "RDFOutputTimesteps")
    {
       inputfilestream >> this->_RDFOutputTimesteps;
    }
    else if(token == "RDFOutputPrefix")
    {
       inputfilestream >> this->_RDFOutputPrefix;
    }
    else if(token == "resultOutputTimesteps")
    {
       inputfilestream >> this->_resultOutputTimesteps;
    }
#ifdef COMPLEX_POTENTIAL_SET
    else if(token == "profiledComponent")
    {
       unsigned cid;
       inputfilestream >> cid;
       cid--;
       this->_domain->considerComponentInProfile(cid);
    }
    else if(token == "profileOutputPrefix")
    {
       inputfilestream >> this->_profileOutputPrefix;
    }
    else if(token == "collectThermostatDirectedVelocity")
    {
       inputfilestream >> this->_collectThermostatDirectedVelocity;
    }
    else if(token == "zOscillator")
    {
       this->_zoscillation = true;
       inputfilestream >> this->_zoscillator;
    }
#endif
    else if(token == "NVE")
    {
       this->_domain->thermostatOff();
    }
  }
  
  
  _domain->readPhaseSpaceData(
    _moleculeContainer
#ifdef TRUNCATED_SHIFTED
    ,
    this->_cutoffRadius
#endif
  );
  _domain->initFarFieldCorr(_cutoffRadius);
  
  // @todo comment
  _integrator = new integrators::Leapfrog(timestepLength);
  

}


void Simulation::initialize(){

  cout << "Initialising rank " << this->_domain->ownrank() << ":\n\n";

  // clear halo
  if(!this->_domain->ownrank()) cout << "   * clearing the halo\n";
  _moleculeContainer->deleteOuterParticles();
  
  if(!this->_domain->ownrank()) cout << "   * updating domain decomposition\n";
  updateParticleContainerAndDecomposition();

  // Force calculation
  if(!this->_domain->ownrank()) cout << "   * force calculation\n";
  _moleculeContainer->traversePairs();
 
  // clear halo
  if(!this->_domain->ownrank()) cout << "   * clearing the halo\n";
  _moleculeContainer->deleteOuterParticles();

  if(this->_doRecordRDF) this->_domain->resetRDF();

#ifdef COMPLEX_POTENTIAL_SET
    cout << "   [:[:]:]   ";
    if(this->_domain->isAcceleratingUniformly())
    {
      if(!this->_domain->ownrank()) cout << "   * initialising uniform acceleration.\n";
      unsigned long uCAT = this->_domain->getUCAT();
      if(!this->_domain->ownrank()) cout << "        uCAT: " << uCAT << " steps.\n";
      /*
      this->_domain->determineAdditionalAcceleration( this->_domainDecomposition,
                                                      this->_moleculeContainer,
                                                      0.000001 * uCAT
                                                         * _integrator->getTimestepLength() );
      this->_integrator->accelerateInstantaneously(this->_moleculeContainer, this->_domain);
      */
      cout << "   [-[-]-]   ";
      this->_domain->determineAdditionalAcceleration( this->_domainDecomposition,
                                                      this->_moleculeContainer,
                                                      uCAT * _integrator->getTimestepLength() );
      if(!this->_domain->ownrank()) cout << "        uniform acceleration initialised.\n";
    }
#endif

  //! @todo calculation of macroscopic values, so that the output in
  //!       step 1 is correct. This doesn't work yet as for the method
  //!       _domain->calculateGlobalValues(...), the iterator has
  //!       to be executed before (sets summv2 and sumIw2)
  cout << "   [.[.].]   ";
#ifdef COMPLEX_POTENTIAL_SET
  if(!this->_domain->ownrank()) cout << "   * calculating directed velocity (if necessary)\n";
  _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
#endif
  if(!this->_domain->ownrank()) cout << "   * calculating velocity sums\n";
  _domain->calculateVelocitySums(_moleculeContainer);
  
  if(!this->_domain->ownrank()) cout << "   * calculating global values\n";
#ifdef COMPLEX_POTENTIAL_SET
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, true);
#else
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer);
#endif

#ifdef COMPLEX_POTENTIAL_SET
  if(this->_zoscillation)
  {
#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * initialize z-oscillators\n";
#endif
    this->_integrator->init1D(this->_zoscillator, this->_moleculeContainer);
  }
#endif

  // initialize output
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain); 
  }

  if(!this->_domain->ownrank()) cout << "system is initialised.\n\n";
}


void Simulation::simulate(){

  Molecule* tM;
  cout << "Now ";
  int ownrank = _domainDecomposition->getRank();
  cout << "simulating rank " << ownrank << ".\n";
  if(ownrank==0) _log.info("simulation(...)", "Starting Simulation: ");
  //double runclock = clock();

  initialize();
  // MAIN LOOP

#ifdef COMPLEX_POTENTIAL_SET
  unsigned uCAT = this->_domain->getUCAT();
#endif

  for( unsigned long simstep = (unsigned long)(this->_domain->getCurrentTime() / _integrator->getTimestepLength());
       simstep <= this->_numberOfTimesteps;
       simstep++ )
  {
    
#ifndef NDEBUG
    cout << "\ntimestep " << simstep << " [rank " << this->_domain->ownrank() << "]:\n";
#endif

    _integrator->eventNewTimestep(_moleculeContainer, _domain);

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * container and decomposition\n";
#endif
    
    updateParticleContainerAndDecomposition();

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * traverse pairs\n";
#endif

    // Force calculation
    _moleculeContainer->traversePairs();

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * delete outer particles\n";
#endif

    // clear halo
    _moleculeContainer->deleteOuterParticles();

   if(this->_doRecordRDF)
   {
      this->_particlePairsHandler->recordRDF();
      this->_moleculeContainer->countParticles(this->_domain);
      this->_domain->tickRDF();
   }

#ifdef COMPLEX_POTENTIAL_SET
    if(!(simstep % _collectThermostatDirectedVelocity))
       _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
#endif

#ifdef COMPLEX_POTENTIAL_SET
    if(this->_domain->isAcceleratingUniformly())
    {
      if(!(simstep % uCAT))
      {
#ifndef NDEBUG
        if(!this->_domain->ownrank()) cout << "   * determine the additional acceleration\n";
#endif
        this->_domain->determineAdditionalAcceleration( this->_domainDecomposition,
                                                        this->_moleculeContainer,
                                                        uCAT * _integrator->getTimestepLength() );
      }
#ifndef NDEBUG
      if(!this->_domain->ownrank()) cout << "   * process the uniform acceleration\n";
#endif
      this->_integrator->accelerateUniformly(this->_moleculeContainer, this->_domain);
    }

    if(this->_zoscillation)
    {
#ifndef NDEBUG
      if(!this->_domain->ownrank()) cout << "   * alert z-oscillators\n";
#endif
      this->_integrator->zOscillation(this->_zoscillator, this->_moleculeContainer);
    }
#endif

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * inform the integrator\n";
#endif

    // Inform the integrator about the calculated forces
    _integrator->eventForcesCalculated(_moleculeContainer, _domain);

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * calculate macroscopic values\n";
#endif

    // calculate the global macroscopic values from the local values
#ifdef COMPLEX_POTENTIAL_SET
    _domain->calculateGlobalValues( _domainDecomposition, _moleculeContainer,
                                    (!(simstep % _collectThermostatDirectedVelocity)) );
#else
    _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer);
#endif

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * velocity scaling\n";
#endif

    if(!this->_domain->NVE())
    {
       // scale velocity and angular momentum
       // @todo why here? what about preF
       if(this->_domain->severalThermostats())
       {
         for( tM = _moleculeContainer->begin();
              tM != _moleculeContainer->end();
              tM = _moleculeContainer->next() )
         {
           int thermostat = this->_domain->getThermostat(tM->componentid());
           if(0 >= thermostat) continue;
#ifdef COMPLEX_POTENTIAL_SET
           if(this->_domain->thermostatIsUndirected(thermostat))
           {
              tM->scale_v( _domain->getGlobalBetaTrans(thermostat),
                           _domain->thermostatv(thermostat, 0),
                           _domain->thermostatv(thermostat, 1),
                           _domain->thermostatv(thermostat, 2)  );
           }
           else
           {
#endif
             tM->scale_v(_domain->getGlobalBetaTrans(thermostat));
#ifdef COMPLEX_POTENTIAL_SET
           }
#endif
           tM->scale_D(_domain->getGlobalBetaRot(thermostat));  
         }
       }
       else
       {
         for( tM = _moleculeContainer->begin();
              tM != _moleculeContainer->end();
              tM = _moleculeContainer->next() )
         {
           tM->scale_v(_domain->getGlobalBetaTrans());
           tM->scale_D(_domain->getGlobalBetaRot());      
         }
       }
   }

    _domain->advanceTime(_integrator->getTimestepLength());
    
    output(simstep);
  }
  
  _moleculeContainer->deleteOuterParticles();

  ostringstream osstrm;
  osstrm.str("");
  osstrm << this->_outputPrefix;
  // osstrm.width(3);
  // osstrm.fill('0');
  // osstrm << this->_domain->ownrank();
  osstrm << ".restart.xdr";
  _domain->writeCheckpoint(osstrm.str(), _moleculeContainer, _domainDecomposition, _integrator->getTimestepLength());
  
  // finish output
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain); 
    delete (*outputIter);
  }
  
  delete _domainDecomposition;
  delete _domain;
  delete _particlePairsHandler;
  delete _moleculeContainer;
  delete _integrator;
  // wait for all processes to reach the end of the program
  //_DomainDecomposition->barrier();
}

void Simulation::output(int simstep){
  int ownrank = _domainDecomposition->getRank();
  
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++)
  {
    (*outputIter)->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep); 
  }

  if(this->_doRecordRDF && !(simstep % this->_RDFOutputTimesteps))
  {
    this->_domain->collectRDF(this->_domainDecomposition);
    if(!ownrank)
    {
      this->_domain->accumulateRDF();
      for(unsigned i=0; i < this->_numberOfComponents; i++)
      {
        for(unsigned j=i; j < this->_numberOfComponents; j++)
        {
          ostringstream osstrm;
          osstrm << this->_RDFOutputPrefix << "_" << i << "-" << j << ".";
          osstrm.fill('0');
          osstrm.width(9);
          osstrm << right << simstep;
          this->_domain->outputRDF(osstrm.str().c_str(), i, j);
          osstrm.str(""); osstrm.clear();
        }
      }
    }
    this->_domain->resetRDF();
  }
#ifdef COMPLEX_POTENTIAL_SET
  if(this->_doRecordProfile && !(simstep % this->_profileRecordingTimesteps))
  {
    this->_domain->recordProfile(this->_moleculeContainer);
  }
  if(this->_doRecordProfile && !(simstep % this->_profileOutputTimesteps))
  {
    this->_domain->collectProfile(this->_domainDecomposition);
    if(!ownrank)
    {
      ostringstream osstrm;
      osstrm << this->_profileOutputPrefix << ".";
      osstrm.fill('0');
      osstrm.width(9);
      osstrm << right << simstep;
      this->_domain->outputProfile(osstrm.str().c_str());
      osstrm.str(""); osstrm.clear();
    }
    this->_domain->resetProfile();
  }
#endif

  if(!(simstep % this->_restartOutputInterval))
  {
     _moleculeContainer->deleteOuterParticles();
     ostringstream osstrm;
     osstrm.str("");
     osstrm << this->_outputPrefix;
     // osstrm.width(3);
     // osstrm.fill('0');
     // osstrm << this->_domain->ownrank();
     osstrm << ".restart.xdr";
     _domain->writeCheckpoint(osstrm.str(), _moleculeContainer, _domainDecomposition, _integrator->getTimestepLength());
  }

  if(ownrank==0)
  {
    cout << simstep << "\t" << _domain->getAverageGlobalUpot() << "\t"
                    << _domain->getGlobalPressure();
#ifdef COMPLEX_POTENTIAL_SET
    double a = _domain->getDirectedVelocity(1);
    if(a > 0)
    {
       cout.precision(4);
       cout << "\t\t" << _domain->getDirectedVelocity(1, 2)
            << "\t" << _domain->getDirectedVelocity(1)
            << "\t\t" << _domain->getUniformAcceleration(1, 2)
            << "\t" << _domain->getUniformAcceleration(1)
            << "\t\t" << _domain->getCosetN(1);
    }
#endif
    cout << "\n";
  }
}

void Simulation::updateParticleContainerAndDecomposition(){
  
  _domainDecomposition->exchangeMolecules(
    _moleculeContainer, _domain->getComponents(), _domain, _cutoffRadius
  );

  // The cache of the molecules must be updated/build after the exchange process,
  // as the cache itself isn't transferred
  Molecule* tempMolecule;
  for(tempMolecule = _moleculeContainer->begin(); tempMolecule != _moleculeContainer->end(); tempMolecule = _moleculeContainer->next()){
    tempMolecule->upd_cache();
  }

  // The particles have moved, so the neighbourhood relations have
  // changed and have to be adjusted
  _moleculeContainer->update();
  
}
