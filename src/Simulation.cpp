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

  this->_cutoffRadius = 0;
  this->_LJCutoffRadius = 0;
  this->_resultOutputTimesteps = 25;
#ifdef COMPLEX_POTENTIAL_SET
  this->_doRecordProfile = false;
  this->_profileRecordingTimesteps = 7;
  this->_profileOutputTimesteps = 12500;
  this->_profileOutputPrefix = "out";
  this->_collectThermostatDirectedVelocity = 100;
  this->_zoscillation = false;
  this->_zoscillator = 1536;
  this->_oscillation = false;
  this->_oscillator = 1536;
  this->_wallLJ = true;
#else
  this->_wallLJ = false;
#endif 
  this->_doRecordRDF = false;
  this->_RDFOutputTimesteps = 25000;
  this->_RDFOutputPrefix = "out";

  this->_initCanonical = 5000;
#ifdef GRANDCANONICAL
  this->_initGrandCanonical = 10000;
  this->h = 0;
#endif
  this->_initStatistics = 20000;
  this->_doAlignCentre = false;
  this->_doCancelMomentum = false;
  this->_stillinger = 0.0;
  this->_TWFcoordination = 0;
  
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
    if(ownrank == 0) cerr << "Not a mardynconfig file! " << token << endl;
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
      if(this->_LJCutoffRadius == 0.0) _LJCutoffRadius = this->_cutoffRadius;
      _domain->readPhaseSpaceHeader(
         timestepLength
#ifdef TRUNCATED_SHIFTED
            ,
         this->_LJCutoffRadius
#endif
#ifdef COMPLEX_POTENTIAL_SET
            ,
         this->_tersoffCutoffRadius
#endif
      );
      this->_numberOfComponents = _domain->getComponents().size();

      _domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius, _wallLJ);
    }
    else if(token=="timestepLength") {
      inputfilestream >> timestepLength;
    }
    else if(token=="cutoffRadius")   {
      inputfilestream >> _cutoffRadius;
    } 
    else if(token=="LJCutoffRadius")   {
      inputfilestream >> _LJCutoffRadius;
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
        
        if(_LJCutoffRadius == 0) _LJCutoffRadius = _cutoffRadius;
        cout.precision(5);
        if(!ownrank) cout << "Fluid cutoff radii: " << _LJCutoffRadius << " (potential), " << _cutoffRadius << " (RDF / electrostatics)\n";
        _moleculeContainer = new datastructures::LinkedCells<Molecule>(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius,
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
    else if(token == "zetaFlow")
    {
       double zeta;
       inputfilestream >> zeta;
       this->_domain->setZetaFlow(zeta);
    }
    else if(token == "tauPrimeFlow")
    {
       double tauPrime;
       inputfilestream >> tauPrime;
       if(timestepLength == 0.0)
       {
          cout << "timestep missing.\n";
          exit(1);
       }
       this->_domain->specifyTauPrime(tauPrime, timestepLength);
    }
    else if(token == "flowControl")
    {
       double q0[3];
       double q1[3];
       inputfilestream >> q0[0] >> q0[1] >> q0[2] >> token;
       if(token != "to")
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n";
             cout << "Expected 'to' instead of '" << token << "'.\n\n";
             cout << "Syntax: flowControl <x0> <y0> <z0> to <x1> <y1> <z1>\n";
          }
          exit(1);
       }
       inputfilestream >> q1[0] >> q1[1] >> q1[2];
       this->_domain->specifyFlowControl(q0, q1);
    }
    else if(token == "wallLJ")
    {
       inputfilestream >> token;
       if(token == "on") this->_wallLJ = true;
       else if(token == "off") this->_wallLJ = false;
       else
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n"
                  << "Expected 'on' or 'off' instead of '" << token << "'.\n\n";
             cout << "Syntax: wallLJ [on|off]\n";
          }
          exit(1);
       }
       _domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius, _wallLJ);
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
    else if(token == "esfera")
    {
       this->_domain->esfera();
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
    else if(token == "TWF")
    {
       inputfilestream >> this->_stillinger >> this->_TWFcoordination;
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
    else if(token == "oscillator")
    {
       this->_oscillation = true;
       inputfilestream >> this->_oscillator;
    }
#endif
#ifdef GRANDCANONICAL
    // chemicalPotential <mu> component <cid> [control <x0> <y0> <z0>
    // to <x1> <y1> <z1>] conduct <ntest> tests every <nstep> steps
    else if(token == "chemicalPotential")
    {
       double imu;
       inputfilestream >> imu;
       inputfilestream >> token;
       if(token != "component")
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n";
             cout << "Expected 'component' instead of '" << token << "'.\n\n";
             cout << "Syntax: chemicalPotential <mu> component <cid> "
	          << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		  << "conduct <ntest> tests every <nstep> steps\n";
          }
          exit(1);
       }
       unsigned icid;
       inputfilestream >> icid;
       icid--;
       inputfilestream >> token;
       double x0, y0, z0, x1, y1, z1;
       bool controlVolume = false;
       if(token == "control")
       {
	  controlVolume = true;
          inputfilestream >> x0 >> y0 >> z0;
          inputfilestream >> token;
          if(token != "to")
          {
             if(ownrank == 0)
             {
                cout << "Input failure.\n";
                cout << "Expected 'to' instead of '" << token << "'.\n\n";
                cout << "Syntax: chemicalPotential <mu> component <cid> "
	             << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		     << "conduct <ntest> tests every <nstep> steps\n";
             }
             exit(1);
          }
          inputfilestream >> x1 >> y1 >> z1;
          inputfilestream >> token;
       }
       if(token != "conduct")
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n";
             cout << "Expected 'conduct' instead of '" << token << "'.\n\n";
             cout << "Syntax: chemicalPotential <mu> component <cid> "
	          << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		  << "conduct <ntest> tests every <nstep> steps\n";
          }
          exit(1);
       }
       unsigned intest;
       inputfilestream >> intest;
       inputfilestream >> token;
       if(token != "tests")
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n";
             cout << "Expected 'tests' instead of '" << token << "'.\n\n";
             cout << "Syntax: chemicalPotential <mu> component <cid> "
	          << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		  << "conduct <ntest> tests every <nstep> steps\n";
          }
          exit(1);
       }
       inputfilestream >> token;
       if(token != "every")
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n";
             cout << "Expected 'every' instead of '" << token << "'.\n\n";
             cout << "Syntax: chemicalPotential <mu> component <cid> "
	          << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		  << "conduct <ntest> tests every <nstep> steps\n";
          }
          exit(1);
       }
       unsigned instep;
       inputfilestream >> instep;
       inputfilestream >> token;
       if(token != "steps")
       {
          if(ownrank == 0)
          {
             cout << "Input failure.\n";
             cout << "Expected 'steps' instead of '" << token << "'.\n\n";
             cout << "Syntax: chemicalPotential <mu> component <cid> "
	          << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
		  << "conduct <ntest> tests every <nstep> steps\n";
          }
          exit(1);
       }
       ensemble::ChemicalPotential tmu = ensemble::ChemicalPotential();
       tmu.setMu(icid, imu);
       tmu.setInterval(instep); 
       tmu.setInstances(intest);
       if(controlVolume) tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
       cout.precision(6);
       if(!ownrank)
       {
	  cout << "chemical Potential " << imu << " component "
	       << icid+1 << " (internally " << icid << ") conduct "
	       << intest << " tests every " << instep << " steps: ";
       }
       this->_lmu.push_back(tmu);
       if(!ownrank)
       {
	  cout << " pushed back.\n";
       }
    }
    else if(token == "planckConstant")
    {
       inputfilestream >> this->h;
    }
#else
    else if(token == "NVE")
    {
       this->_domain->thermostatOff();
    }
#endif
    else if(token == "initCanonical")
    {
       inputfilestream >> this->_initCanonical;
    }
#ifdef GRANDCANONICAL
    else if(token == "initGrandCanonical")
    {
       inputfilestream >> this->_initGrandCanonical;
    }
#endif
    else if(token == "initStatistics")
    {
       inputfilestream >> this->_initStatistics;
    }
    else if(token == "nomomentum")
    {
       this->_doCancelMomentum = true;
       inputfilestream >> this->_momentumInterval;
    }
    else if(token == "AlignCentre")
    {
       this->_doAlignCentre = true;
       inputfilestream >> _alignmentInterval >> _alignmentCorrection;
    }
  }
  
  if(this->_LJCutoffRadius == 0.0) _LJCutoffRadius = this->_cutoffRadius;
  unsigned long maxid = _domain->readPhaseSpaceData(
    _moleculeContainer
#ifdef TRUNCATED_SHIFTED
    ,
    this->_LJCutoffRadius
#endif
#ifdef GRANDCANONICAL
    ,
    &(this->_lmu)
#endif
  );
  if(!ownrank) cout << "Maximal molecule ID: " << maxid << ".\n";
  _domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);
  
#ifdef GRANDCANONICAL
  unsigned idi = this->_lmu.size();
  unsigned j = 0;
  std::list< ensemble::ChemicalPotential >::iterator cpit;
  for(cpit = this->_lmu.begin(); cpit != this->_lmu.end(); cpit++)
  {
     cpit->setIncrement(idi);
     double tmp_molecularMass = _domain->getComponents()[
        cpit->getComponentID()
     ].m();
     cpit->setSystem(
        _domain->getGlobalLength(0), _domain->getGlobalLength(1),
        _domain->getGlobalLength(2), tmp_molecularMass
     );
     cpit->setGlobalN(
        this->_domain->N(cpit->getComponentID())
     );
     cpit->setNextID(j + (int)(1.001 * (256 + maxid)));

     cpit->setSubdomain(
        ownrank, _moleculeContainer->getBoundingBoxMin(0),
        _moleculeContainer->getBoundingBoxMax(0),
        _moleculeContainer->getBoundingBoxMin(1),
        _moleculeContainer->getBoundingBoxMax(1),
        _moleculeContainer->getBoundingBoxMin(2),
        _moleculeContainer->getBoundingBoxMax(2)
     );
     double Tcur = this->_domain->T(0);
     double Ttar = _domain->severalThermostats()? _domain->targetT(1)
                                                : _domain->targetT(0);
     if((Tcur < 0.85*Ttar) || (Tcur > 1.15*Ttar)) Tcur = Ttar;
     cpit->submitTemperature(Tcur);
     if(this->h != 0.0) cpit->setPlanckConstant(this->h);
     
     j++;
  }
#endif

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

  if(this->_stillinger > 0.0) this->_particlePairsHandler->setStillinger(_stillinger);

#ifdef COMPLEX_POTENTIAL_SET
  if(this->_wallLJ)
  {
     this->_particlePairsHandler->enableWallLJ();
     cout << "LJ interaction between wall atoms of different components is ENABLED now.\n";
  }
  else
  {
     this->_particlePairsHandler->disableWallLJ();
     cout << "LJ interaction between wall atoms is now generally TURNED OFF.\n";
  }
#endif

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
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, true, 1.0);
#else
  _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, 1.0);
#endif

#ifdef GRANDCANONICAL
  double Tcur = this->_domain->T(0);
  double Ttar = this->_domain->severalThermostats()? this->_domain->targetT(1)
                                                   : this->_domain->targetT(0);
  if((Tcur < 0.85*Ttar) || (Tcur > 1.15*Ttar)) Tcur = Ttar;
  
  std::list< ensemble::ChemicalPotential >::iterator cpit;
  if(this->h == 0.0) this->h = sqrt(6.2831853 * Ttar);
  for(cpit = this->_lmu.begin(); cpit != this->_lmu.end(); cpit++)
  {
     cpit->submitTemperature(Tcur);
     cpit->setPlanckConstant(this->h);
  }
#endif

#ifdef COMPLEX_POTENTIAL_SET
  if(this->_zoscillation)
  {
#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * initialize z-oscillators\n";
#endif
    this->_integrator->init1D(this->_zoscillator, this->_moleculeContainer);
  }
  if(this->_oscillation)
  {
#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * initialize oscillators\n";
#endif
    this->_integrator->init0D(this->_oscillator, this->_moleculeContainer);
  }
#endif

  // initialize output
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++){
    (*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain, _TWFcoordination); 
  }
  if((this->_initSimulation > this->_initStatistics) && this->_doRecordRDF) this->_domain->tickRDF();

  if(!this->_domain->ownrank()) cout << "system is initialised.\n\n";
}


void Simulation::simulate()
{
  Molecule* tM;
  cout.precision(3);
  cout.width(7);
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

  this->_initSimulation = (unsigned long)(this->_domain->getCurrentTime()
                             / _integrator->getTimestepLength());

  for( unsigned long simstep = this->_initSimulation;
       simstep <= this->_numberOfTimesteps;
       simstep++ )
  {
#ifndef NDEBUG
    if(this->_domain->ownrank() < 3)
       cout << "\ntimestep " << simstep << " [rank " << this->_domain->ownrank() << "]:\n";
    this->_domainDecomposition->barrier();
#endif
    if((simstep == this->_initStatistics) && this->_doRecordRDF) this->_domain->tickRDF();

#ifdef GRANDCANONICAL
    if(simstep >= this->_initGrandCanonical)
    {
       unsigned j = 0;
       std::list< ensemble::ChemicalPotential >::iterator cpit;
       for( cpit = this->_lmu.begin();
	    cpit != this->_lmu.end();
	    cpit++ )
       {
	  if(!((simstep + 2*j + 3) % cpit->getInterval()))
	  {
#ifndef NDEBUG
             if(!this->_domain->ownrank()) cout << " preparing chemical potential ID " << j
                  << "(" << cpit->getMu() << " for c" << cpit->getComponentID()
                  << " with interval " << cpit->getInterval() << ".";
#endif
	     cpit->prepareTimestep(
	        this->_moleculeContainer, this->_domainDecomposition
             );
	  }
	  j++;
       }
    }
#endif

    _integrator->eventNewTimestep(_moleculeContainer, _domain);

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * container and decomposition\n";
    this->_domainDecomposition->barrier();
#endif

   /*
    * halo must be absent to facilitate correct centre of mass calculation
    */
   if(this->_doAlignCentre && !(simstep % _alignmentInterval))
   {
      this->_domain->determineShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
   }

   /*
    * includes creation of the halo
    */
   updateParticleContainerAndDecomposition();

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * traverse pairs\n";
#endif

    // Force calculation
    _moleculeContainer->traversePairs();

#ifdef GRANDCANONICAL
    if(simstep >= _initGrandCanonical)
    {
       unsigned j = 0;
       std::list< ensemble::ChemicalPotential >::iterator cpit;
       for( cpit = this->_lmu.begin();
	    cpit != this->_lmu.end();
	    cpit++ )
       {
	  if(!((simstep + 2*j + 3) % cpit->getInterval()))
	  {     
             if(!this->_domain->ownrank()) 
             {
#ifndef NDEBUG
                cout << "   * grand canonical ensemble(" << j
		     << "): test deletions and insertions\n";
#endif
             }
             this->_moleculeContainer->grandcanonicalStep(
	        &(*cpit), this->_domain->T(0)
             );
	     cpit->assertSynchronization(this->_domainDecomposition);
	     
	     int localBalance = _moleculeContainer->localGrandcanonicalBalance();
	     int balance = _moleculeContainer->grandcanonicalBalance(
	        this->_domainDecomposition
	     );
#ifndef NDEBUG
             if(!this->_domain->ownrank()) 
             {
                cout << "   b[" << ((balance > 0)? "+": "") << balance
                     << "(" << ((localBalance > 0)? "+": "")
		     << localBalance << ")" << " / c = "
		     << cpit->getComponentID() << "]   ";
             }
#endif	     
             this->_domain->Nadd(
	        cpit->getComponentID(), balance, localBalance
	     );
	  }
	  
	  j++;
       }
    }
#endif

   /*
    * halo must be present to permit shifting the particles
    */
   if(this->_doAlignCentre && !(simstep % _alignmentInterval))
   {
      this->_domain->realign(_moleculeContainer);
   }
    
#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * delete outer particles\n";
#endif

    // clear halo
    _moleculeContainer->deleteOuterParticles();

#ifdef GRANDCANONICAL
    if(simstep >= _initGrandCanonical)
    {
       this->_domain->evaluateRho
       (
          this->_moleculeContainer->getNumberOfParticles(),
          this->_domainDecomposition
       );
    }
#endif

#ifdef COMPLEX_POTENTIAL_SET
    if(!(simstep % _collectThermostatDirectedVelocity))
       _domain->calculateThermostatDirectedVelocity(_moleculeContainer);
#endif

   if(this->_doCancelMomentum && !(simstep % _momentumInterval))
   {
      this->_domain->cancelMomentum(_domainDecomposition, _moleculeContainer);
   }

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
      this->_domain->adjustTau(this->_integrator->getTimestepLength());
    }

    if(simstep >= this->_initStatistics)
    {
#ifdef GRANDCANONICAL
       if(this->_lmu.size() == 0)
       {      
#endif
          this->_domain->record_cv_and_sigp();
#ifdef GRANDCANONICAL
       }
#endif
       if(this->_doRecordRDF)
       {
          this->_domain->tickRDF();
          this->_particlePairsHandler->recordRDF();
          this->_moleculeContainer->countParticles(this->_domain);
       }
    }

    if(this->_zoscillation)
    {
#ifndef NDEBUG
      if(!this->_domain->ownrank()) cout << "   * alert z-oscillators\n";
#endif
      this->_integrator->zOscillation(this->_zoscillator, this->_moleculeContainer);
    }
    if(this->_oscillation)
    {
#ifndef NDEBUG
      if(!this->_domain->ownrank()) cout << "   * alert oscillators\n";
#endif
      this->_integrator->oscillation(this->_oscillator, this->_moleculeContainer);
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
                                    (!(simstep % _collectThermostatDirectedVelocity)),
                                    Tfactor(simstep) );
#else
    _domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, Tfactor(simstep));
#endif

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "   * velocity scaling\n";
    this->_domainDecomposition->barrier();
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

#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << "            => advancing time: ";
#endif
    _domain->advanceTime(_integrator->getTimestepLength());
#ifndef NDEBUG
    if(!this->_domain->ownrank()) cout << " done.\n";
    this->_domainDecomposition->barrier();
#endif
    
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
  if(!this->_domain->ownrank()) cout << "Writing XDR checkpoint.\n";
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

void Simulation::output(unsigned long simstep){
  int ownrank = _domainDecomposition->getRank();
  
  std::list<md_io::OutputBase*>::iterator outputIter;
  for(outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++)
  {
#ifndef NDEBUG
     if(!ownrank) cout << "            => doOutput: ";
#endif
     (*outputIter)->doOutput(
        _moleculeContainer, _domainDecomposition, _domain, simstep
#ifdef GRANDCANONICAL
        ,
        &(this->_lmu)
#endif
     ); 
#ifndef NDEBUG
     if(!ownrank) cout << "done.\n";
#endif
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
#ifndef NDEBUG
     if(!ownrank) cout << "            => profile operations: ";
#endif
  if((simstep >= this->_initStatistics) && this->_doRecordProfile
                                        && !(simstep % this->_profileRecordingTimesteps))
  {
    this->_domain->recordProfile(this->_moleculeContainer);
  }
  if((simstep >= this->_initStatistics) && this->_doRecordProfile 
                                        && !(simstep % this->_profileOutputTimesteps))
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
#ifndef NDEBUG
     if(!ownrank) cout << "done.\n";
#endif
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
     if(ownrank==0) cout << "Writing regular XDR output.\n";
     _domain->writeCheckpoint(osstrm.str(), _moleculeContainer, _domainDecomposition, _integrator->getTimestepLength());
  }

  if(ownrank==0)
  {
    cout << (_domain->thermostatWarning()? "*": "")
         << simstep << "\t" << _domain->T(0)
                    << "\t" << _domain->getAverageGlobalUpot() << "\t"
                    << _domain->getGlobalPressure();
#ifdef COMPLEX_POTENTIAL_SET
    double a = _domain->getDirectedVelocity(1);
    if(a > 0)
    {
       cout.precision(4);
       cout << "\t\t" << _domain->getDirectedVelocity(1, 2) << " " << _domain->getDirectedVelocity(1)
            << "\t" << _domain->getUniformAcceleration(1, 2)
            << " " << _domain->getUniformAcceleration(1);
    }
#endif
#ifdef GRANDCANONICAL
    cout << "\t";
    std::list< ensemble::ChemicalPotential >::iterator cpit;
    for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++)
    {
       cout << (this->_domain->thermostatWarning()? "*": "")
            << cpit->getGlobalN() << "  ";
    }
#else
    cout << (this->_domain->thermostatWarning()? "  *": "");
#endif
    cout << "\n";
  }
}

void Simulation::updateParticleContainerAndDecomposition(){
  
  _domainDecomposition->exchangeMolecules(
    _moleculeContainer, _domain->getComponents(), _domain, _cutoffRadius
  );

  // The cache of the molecules must be updated/rebuilt after the exchange process,
  // as the cache itself isn't transferred
  Molecule* tempMolecule;
  for(tempMolecule = _moleculeContainer->begin(); tempMolecule != _moleculeContainer->end(); tempMolecule = _moleculeContainer->next()){
    tempMolecule->upd_cache();
  }

  // The particles have moved, so the neighbourhood relations have
  // changed and have to be adjusted
  _moleculeContainer->update();
  
}

double Simulation::Tfactor(unsigned long simstep)
{
   double xi = (double)(simstep - this->_initSimulation) / (double)(this->_initCanonical - this->_initSimulation);
   if((xi < 0.1) || (xi > 0.9)) return 1.0;
   else if(xi < 0.3) return 15.0*xi - 0.5;
   else if(xi < 0.4) return 10.0 - 20.0*xi;
   else if(xi < 0.6) return 2.0;
   else return 4 - 10.0*xi/3.0;
}

