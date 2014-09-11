/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// Simulation.cpp
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>

#define SIMULATION_SRC
#include "Simulation.h"

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/AdaptiveSubCells.h"
#include "parallel/DomainDecompBase.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#include "parallel/KDDecomposition.h"
#include "parallel/KDDecomposition2.h"
#else
#include "parallel/DomainDecompDummy.h"
#endif

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/LJFlopCounter.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"

#include "io/io.h"
#include "io/xmlreader.h"
#include "io/GeneratorFactory.h"

#include "ensemble/GrandCanonical.h"
#include "ensemble/CanonicalEnsemble.h"
#include "ensemble/PressureGradient.h"

#include "RDF.h"

#ifdef STEEREO
#include "utils/SteereoIntegration.h"
#include <steereo/steereoSimSteering.h>
#include <steereo/steereoCouplingSim.h>
#endif

#include "utils/OptionParser.h"
#include "utils/Timer.h"
#include "utils/Logger.h"
#include "utils/Random.h"

#include "io/TcTS.h"

using Log::global_log;
using optparse::OptionParser;
using optparse::OptionGroup;
using optparse::Values;
using namespace std;

Simulation* global_simulation;

Simulation::Simulation() :
	_rdf(NULL),
	_ljFlopCounter(NULL),
	_domainDecomposition(NULL) {

	initialize();
}

Simulation::~Simulation() {
	if (_rdf)
		delete _rdf;
	if (_domainDecomposition)
		delete _domainDecomposition;
	if (_pressureGradient)
		delete _pressureGradient;
	if (_domain)
		delete _domain;
	if (_particlePairsHandler)
		delete _particlePairsHandler;
	if (_moleculeContainer)
		delete _moleculeContainer;
	if (_integrator)
		delete _integrator;
	if (_inputReader)
		delete _inputReader;
	if (_cellProcessor)
		delete _cellProcessor;
	if (_ljFlopCounter)
		delete _ljFlopCounter;
}

void Simulation::exit(int exitcode) {
#ifdef ENABLE_MPI
	// terminate all mpi processes and return exitcode
	MPI_Abort(MPI_COMM_WORLD, exitcode);
#else
	// call global exit
	::exit(exitcode);
#endif
}

void Simulation::readConfigFile(string filename) {
	if (filename.rfind(".xml") == filename.size() - 4) {
		global_log->info() << "command line config file type is XML (*.xml)" << endl;
		initConfigXML(filename);
	} else if (filename.rfind(".cfg") == filename.size() - 4) {
		global_log->info() << "command line config file type is oldstyle (*.cfg)" << endl;
		initConfigOldstyle(filename);
	} else {
		global_log->info() << "command line config file type is unknown: trying oldstyle" << endl;
		initConfigOldstyle(filename);
	}
}


void Simulation::initConfigXML(const string& inputfilename) {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif
	global_log->info() << "init XML config file: " << inputfilename << endl;
	XMLfileUnits inp(inputfilename);
	XMLReader xmlreader(inputfilename);
	
	global_log->debug() << "Input XML:" << endl << string(inp) << endl;
	inp.changecurrentnode("/mardyn");
	
	string version = xmlreader.getVersion();
	global_log->info() << "MarDyn XML config file version: " << version << endl;
	global_log->info() << "Reading simulation parameters ..." << endl;

	global_log->info() << "Reading integrator information..." << endl;
	xmlreader.getIntegrator( _integrator );	
	double timestepLength = 0.0;
	timestepLength = _integrator->getTimestepLength();
	global_log->info() << " timestep:             " << timestepLength << endl;
	if ( timestepLength == 0 ) {
		/* TODO: perform checks before simulation start as there are command line parameters, too */
		global_log->error() << "Undefined timestep length." << endl;
	}
	
	global_log->info() << "Reading parallelization type..." << endl;
	_domainDecompositionType = xmlreader.getDomainDecompositionType();
	switch( _domainDecompositionType) {
		case DUMMY_DECOMPOSITION:
#ifdef ENABLE_MPI
			global_log->warning() << "Dummy domain decomposition does not work in a parallel environment." << endl;
			exit(1);
#else
			_domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
#endif
			break;
#ifdef ENABLE_MPI
		case DOMAIN_DECOMPOSITION:
			_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
			break;
		case KD_DECOMPOSITION:
			/* TODO: remove parameters from KDDecomposition constructor */
			_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, 1.0, 100);
			break;
#endif
		case UNKNOWN_DECOMPOSITION:
			global_log->error() << "Unknown domain decomposition type." << endl;
			exit(1);
			break;
	}
#ifndef ENABLE_MPI
	if( _domainDecompositionType != DUMMY_DECOMPOSITION ) {
		global_log->warning() << "Input demands parallelization, but the current compilation doesn't support parallel execution." << endl;
	}
#endif	
	
	double simBoxLength[3];
	if ( xmlreader.getSimBoxSize( simBoxLength ) ) {
		global_log->info() << " simulation box size:  (" 
			<< simBoxLength[0] << ", "
			<< simBoxLength[1] << ", "
			<< simBoxLength[2] << ")" << endl;
			
		for( int d = 0; d < 3; d++ ) {
			_domain->setGlobalLength( d, simBoxLength[d] );
		}
	}
	
	double temperature = xmlreader.getTemperature();
	_domain->setGlobalTemperature( temperature );
	global_log->info() << "Temperature of the enseble: " << temperature << endl;
		
	if (inp.changecurrentnode("simulation")) {
		string siminpfile, siminptype;
		if (inp.getNodeValue("input", siminpfile)) {
			global_log->info() << "reading input file:\t" << siminpfile << endl;
			// input type="oldstyle" to include oldstyle input files for backward compatibility - only temporary!!!
			if (inp.getNodeValue("input@type", siminptype) && siminptype == "oldstyle") {
				global_log->info() << "         file type:\t" << siminptype << endl;
				initConfigOldstyle(siminpfile);
			}
		}
		
		if (inp.getNodeValueReduced("cutoffs/radiusLJ", _LJCutoffRadius)) {
			_cutoffRadius = _LJCutoffRadius;
			global_log->info() << "dimensionless LJ cutoff radius:\t" << _LJCutoffRadius << endl;
			global_log->info() << "dimensionless cutoff radius:\t" << _cutoffRadius << endl;
		} else {
			global_log->error() << "Cutoff section missing." << endl;
			exit(1);
		}
		
		string pspfile;
		if (inp.getNodeValue("ensemble/phasespacepoint/file", pspfile)) {
			pspfile.insert(0, inp.getDir());
			global_log->info() << "phasespacepoint description file:\t" << pspfile << endl;
			
			string pspfiletype("OldStyle");
			inp.getNodeValue("ensemble/phasespacepoint/file@type", pspfiletype);
			global_log->info() << "       phasespacepoint file type:\t" << pspfiletype << endl;
			if (pspfiletype == "OldStyle") {
				_inputReader = (InputBase*) new InputOldstyle();
				_inputReader->setPhaseSpaceFile(pspfile);
				_inputReader->setPhaseSpaceHeaderFile(pspfile);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else if (pspfiletype == "PartGen") {
				string mode;
				inp.getNodeValue("ensemble/phasespacepoint/file@mode", mode);
				global_log->error() << "NOT IMPLEMENTED YET";
				/*
				 _inputReader = (InputBase*) new PartGen();
				 _inputReader->setPhaseSpaceHeaderFile(pspfile);
				 // PartGen has to modes, a "Homogeneous" mode, where particles
				 // with a homogeneous distribution are created, and a "Cluster" mode,
				 // where droplets are created. Currently, only the cluster mode is supported,
				 // which needs another config line starting with "clusterFile ..." directly
				 // after the config line starting with "phaseSpaceFile"
				 double gasDensity;
				 double fluidDensity;
				 double volPercOfFluid;
				 string clusterFileName;
				 inputfilestream >> token >> gasDensity >> fluidDensity >> volPercOfFluid >> clusterFileName;
				 ((PartGen*) _inputReader)->setClusterFile(gasDensity, fluidDensity, volPercOfFluid, clusterFileName);
				 ((PartGen*) _inputReader)->readPhaseSpaceHeader(_domain, timestepLength);
				 */
			} else if (pspfiletype == "1CLJGen") {
				string mode;
				inp.getNodeValue("ensemble/phasespacepoint/file@mode", mode);
				global_log->error() << "NOT IMPLEMENTED YET";
				/*
				 int N;
				 double T;
				 string line;
				 getline(inputfilestream, line);
				 stringstream lineStream(line);
				 lineStream >> mode >> N >> T;
				 global_log->debug() << "read: mode " << mode << " N " << N << " T " << T << endl;

				 OneCLJGenerator* generator = (OneCLJGenerator*) new OneCLJGenerator(mode, N, T);
				 if (mode == "Homogeneous") {
				 double rho;
				 lineStream >> rho;
				 generator->setHomogeneuosParameter(rho);
				 } else if (mode == "Cluster") {
				 double rho_gas, rho_fluid, fluidVolumePercent, maxSphereVolume, numSphereSizes;
				 lineStream >> rho_gas >> rho_fluid >> fluidVolumePercent >> maxSphereVolume >> numSphereSizes;
				 generator->setClusterParameters(rho_gas, rho_fluid, fluidVolumePercent, maxSphereVolume, numSphereSizes);
				 } else {
				 global_log->error() << "Error in inputfile: OneCLJGenerator option \""<< mode << "\" not supported!" << endl;
				 global_log->error() << " Has to be  \"Homogeneous\"  or \"Cluster\"  " << endl;
				 exit(1);
				 }
				 generator->readPhaseSpaceHeader(_domain, timestepLength);

				 _inputReader = (OneCLJGenerator*) generator;
				 */
			}
			/* TODO: Check this part of the code */
			if (_LJCutoffRadius == 0.0)
				_LJCutoffRadius = _cutoffRadius;
			_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
		}
		
		if (inp.changecurrentnode("algorithm")) {

			string datastructype;
			if (inp.getNodeValue("datastructure@type", datastructype)) {
				global_log->info() << "datastructure to use:\t" << datastructype << endl;
				if (datastructype == "LinkedCells") {
					int cellsInCutoffRadius = 1;
					inp.getNodeValue("datastructure/cellsInCutoffRadius", cellsInCutoffRadius);
					global_log->info() << "LinkedCells cells in cutoff radius:\t" << cellsInCutoffRadius << endl;
					double bBoxMin[3];
					double bBoxMax[3];
					global_log->debug() << "Box dimensions:" << endl;
					for (int d = 0; d < 3; d++) {
						bBoxMin[d] = _domainDecomposition->getBoundingBoxMin(d, _domain);
						bBoxMax[d] = _domainDecomposition->getBoundingBoxMax(d, _domain);
						global_log->debug() << "[" << d << "] " << bBoxMin[d] << " - " << bBoxMax[d] << endl;
					}
					if (_LJCutoffRadius == 0.0)
						_LJCutoffRadius = _cutoffRadius;
					global_log->debug() << "Cutoff: " << _cutoffRadius << "\nCellsInCutoff: " << cellsInCutoffRadius << endl;
					_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius,
						_tersoffCutoffRadius, cellsInCutoffRadius);
				} else if (datastructype == "AdaptiveSubCells") {
					double bBoxMin[3];
					double bBoxMax[3];
					for (int i = 0; i < 3; i++) {
						bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
						bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
					}
					//creates a new Adaptive SubCells datastructure
					if (_LJCutoffRadius == 0.0)
						_LJCutoffRadius = _cutoffRadius;
					_moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius,
						_tersoffCutoffRadius);
				} else {
					global_log->error() << "UNKOWN DATASTRUCTURE: " << datastructype << endl;
					exit(1);
				}
			}
			inp.changecurrentnode("..");
		} // algo-section
		else {
			global_log->info() << "XML config file " << inputfilename << ": no algo section" << endl;
		}
		if (inp.changecurrentnode("output")) {
			string outputPathAndPrefix;
			if (inp.getNodeValue("Resultwriter", outputPathAndPrefix)) {
				// TODO: add the output frequency to the xml!
				_outputPlugins.push_back(new ResultWriter(1, outputPathAndPrefix));
				global_log->debug() << "ResultWriter '" << outputPathAndPrefix << "'.\n";
			}
			inp.changecurrentnode("..");
		} // output-section
		else {
			global_log->info() << "XML config file " << inputfilename << ": no output section" << endl;
		}
		inp.changecurrentnode("..");
	} // simulation-section
	else {
		global_log->error() << "XML config file " << inputfilename << ": no simulation section" << endl;
	}
	

	// read particle data
	unsigned long maxid = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain, _domainDecomposition);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;
	
	global_log->info() << " reading in components from xml ..." << endl;
	/* TODO: save components into domain:
	 * take care of already existing components
	 * take care about parameter streams
	 */
	std::vector<Component>& dcomponents = _domain->getComponents();
	std::vector<Component> components; 
	xmlreader.getComponents( components );
// 	for( unsigned int j = 0; j < components.size(); j++ ) {
// 		global_log->info() << components[j] << endl;
// 	}
	long numComponents = dcomponents.size();
	global_log->info() << " number of components from input: " << numComponents << endl;
	global_log->info() << " number of components in xml (not implemented): " << components.size() << endl;
	
	// TODO: Check this...
// 	_domain->initParameterStreams( _cutoffRadius, _LJCutoffRadius );
	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);
	
	// test new Decomposition
	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		cpit->setIncrement(idi);
		double tmp_molecularMass = _domain->getComponents()[cpit->getComponentID()].m();
		cpit->setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2),
		        tmp_molecularMass);
		cpit->setGlobalN(_domain->N(cpit->getComponentID()));
		cpit->setNextID(j + (int) (1.001 * (256 + maxid)));

		cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0), _moleculeContainer->getBoundingBoxMax(0),
		        _moleculeContainer->getBoundingBoxMin(1), _moleculeContainer->getBoundingBoxMax(1),
		        _moleculeContainer->getBoundingBoxMin(2), _moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _domain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(
		        0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if (h != 0.0)
			cpit->setPlanckConstant(h);

		j++;
	}

}

void Simulation::initConfigOldstyle(const string& inputfilename) {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

	global_log->info() << "init oldstyle config file: " << inputfilename << endl;

	// open filestream to the input file
	ifstream inputfilestream(inputfilename.c_str());
	if (!inputfilestream.is_open()) {
		global_log->error() << "Could not open file " << inputfilename << endl;
		exit(1);
	}

	//  std::string inputPath;
	//  unsigned int lastIndex = inputfilename.find_last_of('/',inputfilename.size()-1);
	//  if (lastIndex == string::npos)
	//    inputPath="";
	//  else
	//    inputPath = inputfilename.substr(0, lastIndex+1);


	// used to store one token of the inputfilestream
	string token;

	double timestepLength;
	unsigned cosetid = 0;

	// The first line of the config file has to contain the token "MDProjectConfig"
	inputfilestream >> token;
	if ((token != "mardynconfig") && (token != "MDProjectConfig")) {
		global_log->error() << "Not a mardynconfig file! First token: " << token << endl;
		exit(1);
	}

	while (inputfilestream) {
		token.clear();
		inputfilestream >> token;
		global_log->debug() << " [[" << token << "]]" << endl;

		if (token.substr(0, 1) == "#") {
			inputfilestream.ignore(std::numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		if (token == "phaseSpaceFile") {
			string phaseSpaceFileFormat;
			inputfilestream >> phaseSpaceFileFormat;

			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			if (phaseSpaceFileFormat == "OldStyle") {
				string phaseSpaceFileName;
				inputfilestream >> phaseSpaceFileName;
				_inputReader = (InputBase*) new InputOldstyle();
				_inputReader->setPhaseSpaceFile(phaseSpaceFileName);
				_inputReader->setPhaseSpaceHeaderFile(phaseSpaceFileName);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else if (phaseSpaceFileFormat == "Generator") {
				global_log->info() << "phaseSpaceFileFormat is Generator!" << endl;
				string generatorName; // name of the library to load
				string inputFile; // name of the input file for the generator

				string line;
				getline(inputfilestream, line);
				stringstream lineStream(line);
				lineStream >> generatorName >> inputFile;
				_inputReader = GeneratorFactory::loadGenerator(generatorName, inputFile);
				_inputReader->readPhaseSpaceHeader(_domain, timestepLength);
			} else {
				global_log->error() << "Don't recognize phasespaceFile reader " << phaseSpaceFileFormat << endl;
				exit(1);
			}
			if (_LJCutoffRadius == 0.0)
				_LJCutoffRadius = _cutoffRadius;
			_domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
		} else if (token == "timestepLength") {
			inputfilestream >> timestepLength;
		} 
		else if (token == "cutoffRadius") {
			inputfilestream >> _cutoffRadius;
		}
		else if (token == "LJCutoffRadius") {
			inputfilestream >> _LJCutoffRadius;
		} 
		else if ((token == "parallelization") || (token == "parallelisation")) {
#ifndef ENABLE_MPI
			global_log->warning()
			        << "Input file demands parallelization, but the current compilation doesn't\n\tsupport parallel execution.\n"
			        << endl;
			inputfilestream >> token;
#else
			inputfilestream >> token;
			if (token=="DomainDecomposition") {
				_domainDecomposition = (DomainDecompBase*) new DomainDecomposition();
			}
			else if(token=="KDDecomposition") {
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition(_cutoffRadius, _domain, 1.0, 100);
			}
			else if(token=="KDDecomposition2") {
				int updateFrequency = 100;
				int fullSearchThreshold = 3;
				string line;
				getline(inputfilestream, line);
				stringstream lineStream(line);
				lineStream >> updateFrequency >> fullSearchThreshold;
				_domainDecomposition = (DomainDecompBase*) new KDDecomposition2(_cutoffRadius, _domain, updateFrequency, fullSearchThreshold);
			}
#endif
		}
		else if (token == "datastructure") {

			if (_domainDecomposition == NULL) {
				global_log->error()
				        << "_domainDecomposition is NULL! Probably you compiled for MPI, but didn't specify line \"parallelization\" before line \"datastructure\"!"
				        << endl;
				exit(1);
			}

			inputfilestream >> token;
			if (token == "LinkedCells") {
				_particleContainerType = LINKED_CELL;
				int cellsInCutoffRadius;
				inputfilestream >> cellsInCutoffRadius;
				double bBoxMin[3];
				double bBoxMax[3];
				for (int i = 0; i < 3; i++) {
					bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
					bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
				}
				if (this->_LJCutoffRadius == 0.0)
					_LJCutoffRadius = this->_cutoffRadius;
				_moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius,
				        _tersoffCutoffRadius, cellsInCutoffRadius);
			} else if (token == "AdaptiveSubCells") {
				_particleContainerType = ADAPTIVE_LINKED_CELL;
				double bBoxMin[3];
				double bBoxMax[3];
				for (int i = 0; i < 3; i++) {
					bBoxMin[i] = _domainDecomposition->getBoundingBoxMin(i, _domain);
					bBoxMax[i] = _domainDecomposition->getBoundingBoxMax(i, _domain);
				}
				//creates a new Adaptive SubCells datastructure
				if (_LJCutoffRadius == 0.0)
					_LJCutoffRadius = _cutoffRadius;
					_moleculeContainer = new AdaptiveSubCells(bBoxMin, bBoxMax, _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius);
			} else {
				global_log->error() << "UNKOWN DATASTRUCTURE: " << token << endl;
				exit(1);
			}
		} 
		else if (token == "output") {
			inputfilestream >> token;
			if (token == "ResultWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new ResultWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "ResultWriter '" << outputPathAndPrefix << "'.\n";
			} else if (token == "XyzWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new XyzWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "XyzWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "CheckpointWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new CheckpointWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "CheckpointWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "PovWriter") {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new PovWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "POVWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "DecompWriter") {
				unsigned long writeFrequency;
				string mode;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> mode >> outputPathAndPrefix;
				_outputPlugins.push_back(new DecompWriter(writeFrequency, mode, outputPathAndPrefix,
				        _numberOfTimesteps, true));
				global_log->debug() << "DecompWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if ((token == "VisittWriter") || (token == "VISWriter")) {
				unsigned long writeFrequency;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VISWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
				global_log->debug() << "VISWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			} else if (token == "VTKWriter") {
#ifdef VTK
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;
				_outputPlugins.push_back(new VTKMoleculeWriter(writeFrequency, outputPathAndPrefix));
				global_log->debug() << "VTKWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
#else
				Log::global_log->error() << std::endl << "VKT-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "VTKGridWriter") {
#ifdef VTK
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (_particleContainerType == LINKED_CELL) {
					_outputPlugins.push_back(new VTKGridWriter(writeFrequency, outputPathAndPrefix, static_cast<LinkedCells&> (*_moleculeContainer)));
					global_log->debug() << "VTKGridWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
				} else {
					global_log->warning() << "VTKGridWriter only supported with LinkedCells!" << std::endl;
					global_log->warning() << "Generating no VTK output for the grid!" << std::endl;
				}
#else
				Log::global_log->error() << std::endl << "VKT-Plotting demanded, but programme compiled without -DVTK!" << std::endl << std::endl;
#endif
			} else if (token == "StatisticsWriter") {
				unsigned long writeFrequency = 0;
				string outputPathAndPrefix;
				inputfilestream >> writeFrequency >> outputPathAndPrefix;

				if (_particleContainerType == LINKED_CELL) {
					_outputPlugins.push_back(new StatisticsWriter(writeFrequency, outputPathAndPrefix, static_cast<LinkedCells&> (*_moleculeContainer)));
					global_log->debug() << "StatisticsWriter " << writeFrequency << " '" << outputPathAndPrefix
					        << "'.\n";
				} else {
					global_log->warning() << "StatisticsWriter only supported with LinkedCells!" << std::endl;
					global_log->warning() << "Generating no statistics output for the grid!" << std::endl;
				}
			}
			// by Stefan Becker <stefan.becker@mv.uni-kl.de>
			// output for the MegaMol Simple Particle Data File Format (*.mmspd)
			else if (token == "MmspdWriter"){
			      unsigned long writeFrequency = 0;
			      string outputPathAndPrefix;
			      inputfilestream >> writeFrequency >> outputPathAndPrefix;
			      _outputPlugins.push_back(new MmspdWriter(writeFrequency, outputPathAndPrefix, _numberOfTimesteps, true));
			      global_log->debug() << "MmspdWriter " << writeFrequency << " '" << outputPathAndPrefix << "'.\n";
			}
		}
		else if (token == "accelerate") {
			cosetid++;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "comp") {
				global_log->error() << "Expected 'comp' instead of '" << token << "'.\n";
				exit(1);
			}
			int cid = 0;
			while (cid >= 0) {
				inputfilestream >> cid;
				if (cid > 0)
					global_log->info() << "acc. for component " << cid << endl;
				cid--;
				_pressureGradient->assignCoset((unsigned) cid, cosetid);
			}
			double v;
			inputfilestream >> v;
			global_log->debug() << "velocity " << v << endl;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "towards") {
				global_log->error() << "Expected 'towards' instead of '" << token << "'.\n";
				exit(1);
			}
			double dir[3];
			double dirnorm = 0;
			for (unsigned d = 0; d < 3; d++) {
				inputfilestream >> dir[d];
				dirnorm += dir[d] * dir[d];
			}
			dirnorm = 1.0 / sqrt(dirnorm);
			for (unsigned d = 0; d < 3; d++)
				dir[d] *= dirnorm;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "within") {
				global_log->error() << "Expected 'within' instead of '" << token << "'.\n";
				exit(1);
			}
			double tau;
			inputfilestream >> tau;
			inputfilestream >> token;
			global_log->debug() << "Found specifier '" << token << "'\n";

			if (token != "from") {
				global_log->error() << "Expected 'from' instead of '" << token << "'.\n";
				exit(1);
			}
			double ainit[3];
			for (unsigned d = 0; d < 3; d++)
				inputfilestream >> ainit[3];
			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			_pressureGradient->specifyComponentSet(cosetid, dir, tau, ainit, timestepLength);
		} 
		else if (token == "constantAccelerationTimesteps") {
			unsigned uCAT;
			inputfilestream >> uCAT;
			_pressureGradient->setUCAT(uCAT);
		}
		else if (token == "zetaFlow") {
			double zeta;
			inputfilestream >> zeta;
			_pressureGradient->setZetaFlow(zeta);
		}
		else if (token == "tauPrimeFlow") {
			double tauPrime;
			inputfilestream >> tauPrime;
			if (timestepLength == 0.0) {
				global_log->error() << "timestep missing." << endl;
				exit(1);
			}
			_pressureGradient->specifyTauPrime(tauPrime, timestepLength);
		}
		else if (token == "profile") {
			unsigned xun, yun, zun;
			inputfilestream >> xun >> yun >> zun;
			_domain->setupProfile(xun, yun, zun);
			_doRecordProfile = true;
		}
		//! by Stefan Becker, token determining a profile recorded in a slab in the x-y surface, in the middle of the box, i.e. 0.5*Lz
		else if (token == "slabProfile"){
			// slab width in which the profile is recorded
			double width;
			// number of increments in the x,y-direction
			unsigned xun, yun;
			inputfilestream >> xun >> yun >> width;
			this->_domain->setupSlabProfile(xun, yun, width);
			_doRecordSlabProfile = true;
		}
		else if (token == "yOffset") {
			double yOffset;
			inputfilestream >> yOffset;
			_domain->sYOffset(yOffset);
		}
		else if (token == "profileRecordingTimesteps") {
			inputfilestream >> _profileRecordingTimesteps;
		}
		else if (token == "profileOutputTimesteps") {
			inputfilestream >> _profileOutputTimesteps;
		}
		// by Stefan Becker <stefan.becker@mv.uni-kl.de>, token determining the corrdinate system of the density profile
		else if (token == "SessileDrop"){
			this->_domain->sesDrop();
		}
		else if (token == "RDF") {
			double interval;
			unsigned bins;
			inputfilestream >> interval >> bins;
			if (_domain->getComponents().size() <= 0) {
				global_log->error() << "PhaseSpaceFile-Specifiation has to occur befor RDF-Token!" << endl;
				exit(-1);
			}
			_rdf = new RDF(interval, bins, _domain->getComponents().size());
			//_domain->setupRDF(interval, bins);
		}
		else if (token == "RDFOutputTimesteps") {
			unsigned int RDFOutputTimesteps;
			inputfilestream >> RDFOutputTimesteps;
			_rdf->setOutputTimestep(RDFOutputTimesteps);
		}
		else if (token == "RDFOutputPrefix") {
			std::string RDFOutputPrefix;
			inputfilestream >> RDFOutputPrefix;
			_rdf->setOutputPrefix(RDFOutputPrefix);
		} 
		else if (token == "profiledComponent") {
			unsigned cid;
			inputfilestream >> cid;
			cid--;
			_domain->considerComponentInProfile(cid);
		} 
		else if (token == "profileOutputPrefix") {
			inputfilestream >> _profileOutputPrefix;
		}
		else if (token == "collectThermostatDirectedVelocity") {
			inputfilestream >> _collectThermostatDirectedVelocity;
		}
		else if (token == "thermostat"){
			inputfilestream >> _thermostatType;
			if(_thermostatType == ANDERSEN_THERMOSTAT){
			  inputfilestream >>_nuAndersen;
			  global_log->info() << "Using the Andersen Thermostat with nu = " << _nuAndersen << "\n";
			}
		}
		else if (token == "zOscillator") {
			_zoscillation = true;
			inputfilestream >> _zoscillator;
		}
		else if(token == "AlignCentre"){	
			_doAlignCentre = true;
			inputfilestream >> _alignmentInterval >> _alignmentCorrection;
		}
		else if(token == "ComponentForYShift"){
		    _componentSpecificAlignment = true;
		    unsigned cidMin, cidMax;
		    inputfilestream  >> cidMin >> cidMax;
		    cidMin--; // since internally the component number is reduced by one, i.e. cid == 1 in the input file corresponds to the internal cid == 0
		    cidMax--;
		    _domain->considerComponentForYShift(cidMin, cidMax);
		}
		else if(token == "nomomentum"){
		      this->_doCancelMomentum = true;
		      inputfilestream >> this->_momentumInterval;
		}
		// chemicalPotential <mu> component <cid> [control <x0> <y0> <z0>
		// to <x1> <y1> <z1>] conduct <ntest> tests every <nstep> steps
		else if (token == "chemicalPotential") {
			double imu;
			inputfilestream >> imu;
			inputfilestream >> token;
			if (token != "component") {
				global_log->error() << "Expected 'component' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
				        << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
				        << endl;
				exit(1);
			}
			unsigned icid;
			inputfilestream >> icid;
			icid--;
			inputfilestream >> token;
			double x0, y0, z0, x1, y1, z1;
			bool controlVolume = false;
			if (token == "control") {
				controlVolume = true;
				inputfilestream >> x0 >> y0 >> z0;
				inputfilestream >> token;
				if (token != "to") {
					global_log->error() << "Expected 'to' instead of '" << token << "'.\n";
					global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
					        << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
					        << "conduct <ntest> tests every <nstep> steps\n";
					exit(1);
				}
				inputfilestream >> x1 >> y1 >> z1;
				inputfilestream >> token;
			}
			if (token != "conduct") {
				global_log->error() << "Expected 'conduct' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
				        << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] "
				        << "conduct <ntest> tests every <nstep> steps\n";
				exit(1);
			}
			unsigned intest;
			inputfilestream >> intest;
			inputfilestream >> token;
			if (token != "tests") {
				global_log->error() << "Expected 'tests' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
				        << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
				        << endl;
				exit(1);
			}
			inputfilestream >> token;
			if (token != "every") {
				global_log->error() << "Expected 'every' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
				        << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
				        << endl;
				exit(1);
			}
			unsigned instep;
			inputfilestream >> instep;
			inputfilestream >> token;
			if (token != "steps") {
				global_log->error() << "Expected 'steps' instead of '" << token << "'.\n";
				global_log->debug() << "Syntax: chemicalPotential <mu> component <cid> "
				        << "[control <x0> <y0> <z0> to <x1> <y1> <z1>] " << "conduct <ntest> tests every <nstep> steps"
				        << endl;
				exit(1);
			}
			ChemicalPotential tmu = ChemicalPotential();
			tmu.setMu(icid, imu);
			tmu.setInterval(instep);
			tmu.setInstances(intest);
			if (controlVolume)
				tmu.setControlVolume(x0, y0, z0, x1, y1, z1);
			global_log->info() << setprecision(6) << "chemical Potential " << imu << " component " << icid + 1
			        << " (internally " << icid << ") conduct " << intest << " tests every " << instep << " steps: ";
			_lmu.push_back(tmu);
			global_log->info() << " pushed back." << endl;
		}
		else if (token == "planckConstant") {
			inputfilestream >> h;
		}
		else if (token == "NVE") {
			_domain->thermostatOff();
		}
		else if (token == "initCanonical") {
			inputfilestream >> _initCanonical;
		}
		else if (token == "initGrandCanonical") {
			inputfilestream >> _initGrandCanonical;
		}
		else if (token == "initStatistics") {
			inputfilestream >> _initStatistics;
		}
		else if (token == "cutoffRadius") {
                        double rc;
                        inputfilestream >> rc;
                        this->setcutoffRadius(rc);
                }
                else if (token == "LJCutoffRadius") {
                        double rc;
                        inputfilestream >> rc;
                        this->setLJCutoff(rc);
                }
                else if (token == "tersoffCutoffRadius") {
                        double rc;
                        inputfilestream >> rc;
                        this->setTersoffCutoff(rc);
                }
                else if (token == "WallFun_LJ_9_3"){
		  double rho_w, sig_w, eps_w, y_off, y_cut, yMirr;
		  unsigned numComponents;
		  _applyWallFun = true;
		  inputfilestream  >> numComponents >> rho_w >> sig_w >> eps_w >> y_off >> y_cut >> yMirr;
		  double *xi_sf = new double[numComponents];
		  double *eta_sf = new double[numComponents];
		  for(unsigned nc = 0; nc < numComponents; nc++ ){
		    inputfilestream >> xi_sf[nc] >> eta_sf[nc];
		  }
		  _wall.initialize(_domain->getComponents(),  rho_w, sig_w, eps_w, xi_sf, eta_sf, y_off, y_cut, yMirr);
		  delete[] xi_sf;
		  delete[] eta_sf;
		}
		else if (token == "NumberOfFluidComponents"){
		    double numFluidComp;
		    inputfilestream >> numFluidComp;
		    _domain->setNumFluidComponents(numFluidComp);
		}
                else {
			if(token != "") global_log->warning() << "Did not process unknown token " << token << endl;
		}
	}

	// read particle data
	maxid = _inputReader->readPhaseSpace(_moleculeContainer, &_lmu, _domain, _domainDecomposition);

	if (this->_LJCutoffRadius == 0.0)
		_LJCutoffRadius = this->_cutoffRadius;
	_domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

	bool lj_present = false;
	bool charge_present = false;
	bool dipole_present = false;
	bool quadrupole_present = false;
	bool tersoff_present = false;

	const vector<Component> components = _domain->getComponents();
	for (size_t i = 0; i < components.size(); i++) {
		lj_present |= components[i].numLJcenters() != 0;
		charge_present |= components[i].numCharges() != 0;
		dipole_present |= components[i].numDipoles() != 0;
		quadrupole_present |= components[i].numQuadrupoles() != 0;
		tersoff_present |= components[i].numTersoff() != 0;
	}
	global_log->info() << "xx lj present: " << lj_present << endl;
	global_log->info() << "xx charge present: " << charge_present << endl;
	global_log->info() << "xx dipole present: " << dipole_present << endl;
	global_log->info() << "xx quadrupole present: " << quadrupole_present << endl;
	global_log->info() << "xx tersoff present: " << tersoff_present << endl;

#if 1
	if (charge_present || dipole_present || quadrupole_present || tersoff_present) {
		_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius, _particlePairsHandler);
	} else {
		_cellProcessor = new VectorizedCellProcessor( *_domain,_LJCutoffRadius);
	}
#else
	_cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius, _particlePairsHandler);
#endif

	// @todo comment
	_integrator = new Leapfrog(timestepLength);

	// test new Decomposition
	_moleculeContainer->update();
	_moleculeContainer->deleteOuterParticles();

	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		cpit->setIncrement(idi);
		double tmp_molecularMass = _domain->getComponents()[cpit->getComponentID()].m();
		cpit->setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2),
		        tmp_molecularMass);
		cpit->setGlobalN(_domain->N(cpit->getComponentID()));
		cpit->setNextID(j + (int) (1.001 * (256 + maxid)));

		cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0), _moleculeContainer->getBoundingBoxMax(0),
		        _moleculeContainer->getBoundingBoxMin(1), _moleculeContainer->getBoundingBoxMax(1),
		        _moleculeContainer->getBoundingBoxMin(2), _moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _domain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if (h != 0.0)
			cpit->setPlanckConstant(h);

		j++;
	}
}

void Simulation::prepare_start() {
	global_log->info() << "Initializing simulation" << endl;
	
	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();
	global_log->info() << "Updating domain decomposition" << endl;
	updateParticleContainerAndDecomposition();
	global_log->info() << "Performing inital force calculation" << endl;

    if(!_cellProcessor) {
        global_log->warning() << "No cell processor initialised. Using Legacy Cell Processor." << endl;
        _cellProcessor = new LegacyCellProcessor( _cutoffRadius, _LJCutoffRadius, _tersoffCutoffRadius, _particlePairsHandler);
    }
	_moleculeContainer->traverseCells(*_cellProcessor);

	_ljFlopCounter = new LJFlopCounter(_LJCutoffRadius);
	_moleculeContainer->traverseCells(*_ljFlopCounter);

	// TODO:
	// here we have to call calcFM() manually, otherwise force and moment are not
	// updated inside the molecule (actually this is done in upd_postF)
	// or should we better call the integrator->eventForcesCalculated?
	for (Molecule* tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
		tM->calcFM();
	}
	// clear halo
	global_log->info() << "Clearing halos" << endl;
	_moleculeContainer->deleteOuterParticles();

	// initialize the radial distribution function
	// TODO this call should not be neccessary!
//	if (_rdf != NULL)
//		_rdf->reset();

	if (_pressureGradient->isAcceleratingUniformly()) {
		global_log->info() << "Initialising uniform acceleration." << endl;
		unsigned long uCAT = _pressureGradient->getUCAT();
		global_log->info() << "uCAT: " << uCAT << " steps." << endl;
		_pressureGradient->determineAdditionalAcceleration(_domainDecomposition, _moleculeContainer, uCAT
		        * _integrator->getTimestepLength());
		global_log->info() << "Uniform acceleration initialised." << endl;
	}

	global_log->info() << "Calculating global values" << endl;
	_domain->calculateThermostatDirectedVelocity(_moleculeContainer);

	_domain->calculateVelocitySums(_moleculeContainer);

	_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, true, 1.0);
	global_log->debug() << "Calculating global values finished." << endl;

	if (_lmu.size() > 0) {
		/* TODO: thermostat */
		double Tcur = _domain->getGlobalCurrentTemperature();
		/* FIXME: target temperature from thermostat ID 0 or 1? */
		double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(
		        0);
		if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;

		list<ChemicalPotential>::iterator cpit;
		if (h == 0.0)
			h = sqrt(6.2831853 * Ttar);
		for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
			cpit->submitTemperature(Tcur);
			cpit->setPlanckConstant(h);
		}
	}

	if (_zoscillation) {
		global_log->debug() << "Initializing z-oscillators" << endl;
		_integrator->init1D(_zoscillator, _moleculeContainer);
	}

	// initialize output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		(*outputIter)->initOutput(_moleculeContainer, _domainDecomposition, _domain);
	}

#ifdef STEEREO
	_steer = initSteereo (_domainDecomposition->getRank(), _domainDecomposition->getNumProcs());
#ifdef STEEREO_COUPLING
	_coupling = initCoupling(_steer, (long*) &_simstep);
#endif  /* STEEREO_COUPLING */
	registerSteereoCommands (_steer, this);
	startListeningSteereo (_steer);
#endif  /* STEEREO */

//	if ((_initSimulation > _initStatistics) && this->_rdf != NULL) {
//		this->_rdf->tickRDF();
//		this->_particlePairsHandler->setRDF(_rdf);
//	}

	global_log->info() << "System initialised\n" << endl;
	global_log->info() << "System contains " << _domain->getglobalNumMolecules() << " molecules." << endl;
}

void Simulation::simulate() {

	Molecule* tM;

	global_log->info() << "Started simulation" << endl;

	// (universal) constant acceleration (number of) timesteps
	unsigned uCAT = _pressureGradient->getUCAT();
	_initSimulation = (unsigned long) (_domain->getCurrentTime() / _integrator->getTimestepLength());

	/* demonstration for the usage of the new ensemble class */
	CanonicalEnsemble ensemble(_moleculeContainer, &(_domain->getComponents()));
	ensemble.updateGlobalVariable(NUM_PARTICLES);
	global_log->debug() << "Number of particles in the Ensemble: " << ensemble.N() << endl;
	ensemble.updateGlobalVariable(ENERGY);
	global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E() << endl;
	ensemble.updateGlobalVariable(TEMPERATURE);
	global_log->debug() << "Temperature of the Ensemble: " << ensemble.T() << endl;

	/***************************************************************************/
	/* BEGIN MAIN LOOP                                                         */
	/***************************************************************************/
	// all timers except the ioTimer messure inside the main loop
	Timer loopTimer;
	Timer decompositionTimer;
	Timer perStepIoTimer;
	Timer ioTimer;

#if defined(STEEREO) && defined(STEEREO_COUPLING)
	_simstep = _initSimulation;
	//SteereoLogger::setOutputLevel (4);
	_coupling->waitForConnection();
#endif

	loopTimer.start();
	for (_simstep = _initSimulation; _simstep <= _numberOfTimesteps; _simstep++) {
		if (_simstep >= _initGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
					cpit->prepareTimestep(_moleculeContainer, _domainDecomposition);
				}
				j++;
			}
		}
		global_log->debug() << "timestep " << _simstep << endl;

		_integrator->eventNewTimestep(_moleculeContainer, _domain);

#if defined(STEEREO) && defined(STEEREO_COUPLING)
		_steer -> processQueue (1);
		global_log->debug() << "molecules in simulation: " << _moleculeContainer->getNumberOfParticles() << std::endl;
#endif
		// activate RDF sampling
		if ((_simstep >= this->_initStatistics) && this->_rdf != NULL) {
			this->_rdf->tickRDF();
			this->_particlePairsHandler->setRDF(_rdf);
			this->_rdf->accumulateNumberOfMolecules(_domain->getComponents());
		}
		/*! by Stefan Becker <stefan.becker@mv.uni-kl.de> 
		 *realignment tools borrowed from Martin Horsch, for the determination of the centre of mass 
		 *the halo MUST NOT be present*/
#ifndef NDEBUG 
#ifndef ENABLE_MPI
			unsigned particleNoTest = 0;
                        for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next())
                        particleNoTest++;
                        global_log->info()<<"particles before determine shift-methods, halo not present:" << particleNoTest<< "\n";
#endif
#endif
		 if(_doAlignCentre && !(_simstep % _alignmentInterval))
		{
			if(_componentSpecificAlignment){
				//! !!! the sequence of calling the two methods MUST be: FIRST determineXZShift() THEN determineYShift() !!!
				_domain->determineXZShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
				_domain->determineYShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
			else{
				_domain->determineShift(_domainDecomposition, _moleculeContainer, _alignmentCorrection);
			}
#ifndef NDEBUG 
#ifndef ENABLE_MPI			
			particleNoTest = 0;
			for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) 
			particleNoTest++;
			global_log->info()<<"particles after determine shift-methods, halo not present:" << particleNoTest<< "\n";
#endif
#endif
		}

		// ensure that all Particles are in the right cells and exchange Particles
		global_log->debug() << "Updating container and decomposition" << endl;
		loopTimer.stop();
		decompositionTimer.start();
		updateParticleContainerAndDecomposition(); // includes the creation of the halo
		decompositionTimer.stop();
		loopTimer.start();
		
		// Force calculation
		global_log->debug() << "Traversing pairs" << endl;
		//_moleculeContainer->traversePairs(_particlePairsHandler);
		_moleculeContainer->traverseCells(*_cellProcessor);

		if(_applyWallFun){
		  _wall.calcTSLJ_9_3(_moleculeContainer, _domain);
		}

		// test deletions and insertions
		if (_simstep >= _initGrandCanonical) {
			unsigned j = 0;
			list<ChemicalPotential>::iterator cpit;
			for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
				if (!((_simstep + 2 * j + 3) % cpit->getInterval())) {
					global_log->debug() << "Grand canonical ensemble(" << j << "): test deletions and insertions"
					        << endl;
					/* TODO: thermostat */
					_moleculeContainer->grandcanonicalStep(&(*cpit), _domain->getGlobalCurrentTemperature());
#ifndef NDEBUG
					/* silly check if random numbers inserted are the same for all processes... */
					cpit->assertSynchronization(_domainDecomposition);
#endif

					int localBalance = _moleculeContainer->localGrandcanonicalBalance();
					int balance = _moleculeContainer->grandcanonicalBalance(_domainDecomposition);
					global_log->debug() << "   b[" << ((balance > 0) ? "+" : "") << balance << "(" << ((localBalance
					        > 0) ? "+" : "") << localBalance << ")" << " / c = " << cpit->getComponentID() << "]   "
					        << endl;
					_domain->Nadd(cpit->getComponentID(), balance, localBalance);
				}

				j++;
			}
		}
		/*! by Stefan Becker <stefan.becker@mv.uni-kl.de> 
		  * realignment tools borrowed from Martin Horsch
		  * For the actual shift the halo MUST be present!
		  */
		
		if(_doAlignCentre && !(_simstep % _alignmentInterval))
		{
			_domain->realign(_moleculeContainer);
#ifndef NDEBUG 
#ifndef ENABLE_MPI
			particleNoTest = 0;
			for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) 
			particleNoTest++;
			global_log->info()<<"particles after realign(), halo present:" << particleNoTest<< "\n";
#endif
#endif
		}
		
		// clear halo
		global_log->debug() << "Delete outer particles" << endl;
		_moleculeContainer->deleteOuterParticles();
#ifndef NDEBUG 
#ifndef ENABLE_MPI
			particleNoTest = 0;
			for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) 
			particleNoTest++;
			global_log->info()<<"particles after clearing the halo:" << particleNoTest<< "\n";
#endif				
#endif

		if (_simstep >= _initGrandCanonical) {
			_domain->evaluateRho(_moleculeContainer->getNumberOfParticles(), _domainDecomposition);
		}

		if(_doCancelMomentum && !(_simstep % _momentumInterval)) {
			_domain->cancelMomentum(_domainDecomposition, _moleculeContainer);
		}
		
		if (!(_simstep % _collectThermostatDirectedVelocity))
			_domain->calculateThermostatDirectedVelocity(_moleculeContainer);
		if (_pressureGradient->isAcceleratingUniformly()) {
			if (!(_simstep % uCAT)) {
				global_log->debug() << "Determine the additional acceleration" << endl;
				_pressureGradient->determineAdditionalAcceleration(_domainDecomposition, _moleculeContainer, uCAT
				        * _integrator->getTimestepLength());
			}
			global_log->debug() << "Process the uniform acceleration" << endl;
			_integrator->accelerateUniformly(_moleculeContainer, _domain);
			_pressureGradient->adjustTau(this->_integrator->getTimestepLength());
		}

		/*
		 * radial distribution function
		 */
		if (_simstep >= _initStatistics) {
			if (this->_lmu.size() == 0) {
				this->_domain->record_cv();
			}
//			if (this->_rdf != NULL) {
//				this->_rdf->tickRDF();
//
//			}
		}

		if (_zoscillation) {
			global_log->debug() << "alert z-oscillators" << endl;
			_integrator->zOscillation(_zoscillator, _moleculeContainer);
		}
		// Inform the integrator about the calculated forces
		global_log->debug() << "Inform the integrator" << endl;
		_integrator->eventForcesCalculated(_moleculeContainer, _domain);

		// calculate the global macroscopic values from the local values
		global_log->debug() << "Calculate macroscopic values" << endl;
		_domain->calculateGlobalValues(_domainDecomposition, _moleculeContainer, (!(_simstep
		        % _collectThermostatDirectedVelocity)), Tfactor(_simstep));

		// scale velocity and angular momentum
		if (!_domain->NVE()) {
		  if(_thermostatType == VELSCALE_THERMOSTAT){
			global_log->debug() << "Velocity scaling" << endl;
			if (_domain->severalThermostats()) {
				for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
					int thermostat = _domain->getThermostat(tM->componentid());
					if (0 >= thermostat)
						continue;
					if (_domain->thermostatIsUndirected(thermostat)) {
						/* TODO: thermostat */
						tM->scale_v(_domain->getGlobalBetaTrans(thermostat), _domain->getThermostatDirectedVelocity(
						        thermostat, 0), _domain->getThermostatDirectedVelocity(thermostat, 1),
						        _domain->getThermostatDirectedVelocity(thermostat, 2));
					} else {
						tM->scale_v(_domain->getGlobalBetaTrans(thermostat));
					}
					tM->scale_D(_domain->getGlobalBetaRot(thermostat));
				}
			} else {
				for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
					tM->scale_v(_domain->getGlobalBetaTrans());
					tM->scale_D(_domain->getGlobalBetaRot());
				}
			}
		  }
		  else if(_thermostatType == ANDERSEN_THERMOSTAT){ //! the Andersen Thermostat
//		    global_log->info() << "Andersen Thermostat" << endl;
		    double nuDt = _nuAndersen * _integrator->getTimestepLength();
//		    global_log->info() << "Timestep length = " << _integrator->getTimestepLength() << " nuDt = " << nuDt << "\n";
		    unsigned numPartThermo = 0; // for testing reasons
		    double tTarget;
		    double stdDevTrans, stdDevRot;
		    if(_domain->severalThermostats()) {
		      for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
			if (_rand.rnd() < nuDt){
			  numPartThermo++;
			  int thermostat = _domain->getThermostat(tM->componentid());
			  tTarget = _domain->getTargetTemperature(thermostat);
			  stdDevTrans = sqrt(tTarget/tM->gMass());
			  for(unsigned short d = 0; d < 3; d++){
			    stdDevRot = sqrt(tTarget*tM->getI(d));
			    tM->setv(d,_rand.gaussDeviate(stdDevTrans));
			    tM->setD(d,_rand.gaussDeviate(stdDevRot));
			  }
			}
		      }
		    }
		    else{
		      tTarget = _domain->getTargetTemperature(0);
		      for (tM = _moleculeContainer->begin(); tM != _moleculeContainer->end(); tM = _moleculeContainer->next()) {
			if (_rand.rnd() < nuDt){
			  numPartThermo++;
			  // action of the anderson thermostat: mimic a collision by assigning a maxwell distributed velocity
			  stdDevTrans = sqrt(tTarget/tM->gMass());
			  for(unsigned short d = 0; d < 3; d++){
			    stdDevRot = sqrt(tTarget*tM->getI(d));
			    tM->setv(d,_rand.gaussDeviate(stdDevTrans));
			    tM->setD(d,_rand.gaussDeviate(stdDevRot));
			  }
			}
		      }
		    }
//		    global_log->info() << "Andersen Thermostat: n = " << numPartThermo ++ << " particles thermostated\n";
		  }
		}

		_domain->advanceTime(_integrator->getTimestepLength());
#ifdef STEEREO
		_steer -> processQueue (0);
#endif

		/* BEGIN PHYSICAL SECTION:
		 * the system is in a consistent state so we can extract global variables 
		 */
		ensemble.updateGlobalVariable(NUM_PARTICLES);
		global_log->debug() << "Number of particles in the Ensemble: " << ensemble.N() << endl;
		ensemble.updateGlobalVariable(ENERGY);
		global_log->debug() << "Kinetic energy in the Ensemble: " << ensemble.E() << endl;
		ensemble.updateGlobalVariable(TEMPERATURE);
		global_log->debug() << "Temperature of the Ensemble: " << ensemble.T() << endl;
		/* END PHYSICAL SECTION */

		// measure per timestep IO
		loopTimer.stop();
		perStepIoTimer.start();
		output(_simstep);
		perStepIoTimer.stop();
		loopTimer.start();
	}
	loopTimer.stop();
	/***************************************************************************/
	/* END MAIN LOOP                                                           */
	/***************************************************************************/

	ioTimer.start();
	/* write final checkpoint */
	string cpfile(_outputPrefix + ".restart.xdr");
	_domain->writeCheckpoint(cpfile, _moleculeContainer, _domainDecomposition);
	// finish output
	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		(*outputIter)->finishOutput(_moleculeContainer, _domainDecomposition, _domain);
		delete (*outputIter);
	}
	ioTimer.stop();

	global_log->info() << "Computation in main loop took: " << loopTimer.get_etime() << " sec" << endl;
	global_log->info() << "Decomposition took. " << decompositionTimer.get_etime() << " sec" << endl;
	global_log->info() << "IO in main loop  took:         " << perStepIoTimer.get_etime() << " sec" << endl;
	global_log->info() << "Final IO took:                 " << ioTimer.get_etime() << " sec" << endl;

	unsigned long numTimeSteps = _numberOfTimesteps - _initSimulation + 1; // +1 because of <= in loop
	double elapsed_time = loopTimer.get_etime() + decompositionTimer.get_etime();
	double flop_rate = _ljFlopCounter->getTotalFlopCount() * numTimeSteps / elapsed_time / (1024*1024);
	global_log->info() << "LJ-FLOP-Count per Iteration: " << _ljFlopCounter->getTotalFlopCount() << " FLOPs" <<endl;
	global_log->info() << "FLOP-rate: " << flop_rate << " MFLOPS" << endl;
}

void Simulation::output(unsigned long simstep) {
	int ownrank = _domainDecomposition->getRank();

	std::list<OutputBase*>::iterator outputIter;
	for (outputIter = _outputPlugins.begin(); outputIter != _outputPlugins.end(); outputIter++) {
		(*outputIter)->doOutput(_moleculeContainer, _domainDecomposition, _domain, simstep, &(_lmu));
	}

	if (_rdf != NULL) {
		_rdf->doOutput(_domainDecomposition, _domain, simstep);
	}

	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileRecordingTimesteps)) {
		_domain->recordProfile(_moleculeContainer);
	}
	if ((simstep >= _initStatistics) && _doRecordProfile && !(simstep % _profileOutputTimesteps)) {
		_domain->collectProfile(_domainDecomposition);
		if (!ownrank) {
			ostringstream osstrm;
			osstrm << _profileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			if(this->_domain->isCylindrical()){
				this->_domain->outputCylProfile(osstrm.str().c_str());
			}
			else{
			_domain->outputProfile(osstrm.str().c_str());
			}
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetProfile();
	}
	//! the same profile procedure for the slab profile
	if ((simstep >= _initStatistics) && _doRecordSlabProfile && !(simstep % _profileRecordingTimesteps)) {
		_domain->recordSlabProfile(_moleculeContainer);
	}
	if ((simstep >= _initStatistics) && _doRecordSlabProfile && !(simstep % _profileOutputTimesteps)) {
		_domain->collectProfile(_domainDecomposition);
		if (!ownrank) {
			ostringstream osstrm;
			osstrm << _profileOutputPrefix << ".";
			osstrm.fill('0');
			osstrm.width(9);
			osstrm << right << simstep;
			_domain->outputSlabProfile(osstrm.str().c_str());
			osstrm.str("");
			osstrm.clear();
		}
		_domain->resetProfile();
	}

	if (_domain->thermostatWarning())
		global_log->warning() << "Thermostat!" << endl;
	/* TODO: thermostat */
	global_log->info() << "Simstep = " << simstep << "\tT = " << _domain->getGlobalCurrentTemperature() << "\tU_pot = "
	        << _domain->getAverageGlobalUpot() << "\tp = " << _domain->getGlobalPressure() << endl;
}

void Simulation::finalize() {
#ifdef ENABLE_MPI
	delete _domainDecomposition;
	_domainDecomposition = NULL;
#endif
}

void Simulation::updateParticleContainerAndDecomposition() {

	// The particles have moved, so the neighbourhood relations have
	// changed and have to be adjusted
	_moleculeContainer->update();
	//_domainDecomposition->exchangeMolecules(_moleculeContainer, _domain->getComponents(), _domain);
	_domainDecomposition->balanceAndExchange(true, _moleculeContainer, _domain->getComponents(), _domain);
	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();
}

/* FIXME: we shoud provide a more general way of doing this */
double Simulation::Tfactor(unsigned long simstep) {
	double xi = (double) (simstep - _initSimulation) / (double) (_initCanonical - _initSimulation);
	if ((xi < 0.1) || (xi > 0.9))
		return 1.0;
	else if (xi < 0.3)
		return 15.0 * xi - 0.5;
	else if (xi < 0.4)
		return 10.0 - 20.0 * xi;
	else if (xi < 0.6)
		return 2.0;
	else
		return 4 - 10.0 * xi / 3.0;
}

void Simulation::initialize() {
	int ownrank = 0;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

	global_simulation = this;

	_domainDecomposition = NULL;
	_domain = NULL;
	_particlePairsHandler = NULL;
	_moleculeContainer = NULL;
	_integrator = NULL;
	_inputReader = NULL;


#ifndef ENABLE_MPI
	global_log->info() << "Initializing the alibi domain decomposition ... " << endl;
	_domainDecomposition = (DomainDecompBase*) new DomainDecompDummy();
	global_log->info() << "Initialization done" << endl;
#endif

	/*
	 * default parameters
	 */
	_cutoffRadius = 0.0;
	_LJCutoffRadius = 0.0;
	_tersoffCutoffRadius = 3.0;
	_numberOfTimesteps = 1;
	_outputPrefix = string("mardyn");
	_outputPrefix.append(gettimestring());
	_doRecordProfile = false;
	_doRecordSlabProfile = false;
	_profileRecordingTimesteps = 7;
	_profileOutputTimesteps = 12500;
	_profileOutputPrefix = "out";
	_collectThermostatDirectedVelocity = 100;
	_zoscillation = false;
	_zoscillator = 512;
	_initCanonical = 5000;
	_initGrandCanonical = 10000000;
	_initStatistics = 20000;
	h = 0.0;
	
	_thermostatType = VELSCALE_THERMOSTAT;
	_nuAndersen = 0.0;
	_rand.init(8624);

	_doAlignCentre = false;
	_componentSpecificAlignment = false;
	_alignmentInterval = 25;
	_doCancelMomentum = false;
	_momentumInterval = 1000;
	_applyWallFun = false;

	_pressureGradient = new PressureGradient(ownrank);
	global_log->info() << "Constructing domain ..." << endl;
	_domain = new Domain(ownrank, this->_pressureGradient);
	global_log->info() << "Domain construction done." << endl;
	_particlePairsHandler = new ParticlePairs2PotForceAdapter(*_domain);
}

void Simulation::mkTcTS(int argc, char** argv)
{
   _particleContainerType = LINKED_CELL;

   TcTS(argc, argv, this->_domain, &(this->_domainDecomposition), &(this->_integrator), &(this->_moleculeContainer), &(this->_outputPlugins), this->_rdf, this);

   if (this->_LJCutoffRadius == 0.0)
      _LJCutoffRadius = this->_cutoffRadius;
   _moleculeContainer->update();
   _moleculeContainer->deleteOuterParticles();

   _domain->initParameterStreams(_cutoffRadius, _LJCutoffRadius);
   _domain->initFarFieldCorr(_cutoffRadius, _LJCutoffRadius);

   /*
   unsigned idi = _lmu.size();
   unsigned j = 0;
   std::list<ChemicalPotential>::iterator cpit;
   double Tcur = _domain->getCurrentTemperature(0);
   for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
      cpit->setIncrement(idi);
      double tmp_molecularMass = _domain->getComponents()[cpit->getComponentID()].m();
      cpit->setSystem(_domain->getGlobalLength(0), _domain->getGlobalLength(1), _domain->getGlobalLength(2),
         tmp_molecularMass);
      cpit->setGlobalN(_domain->N(cpit->getComponentID()));
      cpit->setNextID(j + (int) (1.001 * (256 + maxid)));
      cpit->setSubdomain(ownrank, _moleculeContainer->getBoundingBoxMin(0), _moleculeContainer->getBoundingBoxMax(0),
         _moleculeContainer->getBoundingBoxMin(1), _moleculeContainer->getBoundingBoxMax(1),
         _moleculeContainer->getBoundingBoxMin(2), _moleculeContainer->getBoundingBoxMax(2));
      double Ttar = _domain->severalThermostats() ? _domain->getTargetTemperature(1) : _domain->getTargetTemperature(0);
      if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar)) Tcur = Ttar;
      cpit->submitTemperature(Tcur);
      if (h != 0.0) cpit->setPlanckConstant(h);
      j++;
   }
   */
}

