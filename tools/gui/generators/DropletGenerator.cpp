/*
 * DropletGenerator.cpp
 *
 *  Created on: May 14, 2011
 *      Author: kovacevt
 */

#include "DropletGenerator.h"
#include "common/MardynConfigurationParameters.h"
#include "common/PrincipalAxisTransform.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Tokenize.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "common/DropletPlacement.h"
#include "utils/Logger.h"


#include <cmath>
#include <iostream>

extern "C" {

	Generator* create_generator() {
		return new DropletGenerator();
	}

	void destruct_generator(Generator* generator) {
		delete generator;
	}
}


DropletGenerator::DropletGenerator() :
	MDGenerator("DropletGenerator") {
	numOfMolecules = 50;
	_temperature = 300. / 315774.5;
	setClusterParameters(0.05, 0.8, 15, 5, 4);
	_components.resize(1);
	_components[0].addLJcenter(0, 0, 0, 1.0, 1.0, 1.0, 0.0, false);
}

DropletGenerator::~DropletGenerator() {
}

void DropletGenerator::setClusterParameters(double gas, double fluidDen,
		double fluidVol, double maxSphereVol, int numSphereSiz) {
	gasDensity = gas;
	fluidDensity = fluidDen;

	fluidVolume = fluidVol;
	maxSphereVolume = maxSphereVol;
	numSphereSizes = numSphereSiz;

	rho = fluidDensity * fluidVolume / 100.0 + gasDensity * (1 - fluidVolume
			/ 100.0);
	simBoxLength[0] = pow(numOfMolecules / rho, (1.0 / 3.0));
	simBoxLength[1] = pow(numOfMolecules / rho, (1.0 / 3.0));
	simBoxLength[2] = pow(numOfMolecules / rho, (1.0 / 3.0));

}

void DropletGenerator::readPhaseSpaceHeader(Domain* domain, double timestep) {
		domain->setCurrentTime(0);

		domain->disableComponentwiseThermostat();
		domain->setGlobalTemperature(_temperature);
		domain->setGlobalLength(0, simBoxLength[0]);
		domain->setGlobalLength(1, simBoxLength[1]);
		domain->setGlobalLength(2, simBoxLength[2]);

		_logger->debug() << "DropletGenerator: set global length=[" << simBoxLength[0]
		    << "," << simBoxLength[1] << "," <<  simBoxLength[2] << "]" << endl;

		for (unsigned int i = 0; i < _components.size(); i++) {
			Component component = _components[i];
			if (_configuration.performPrincipalAxisTransformation()) {
				principalAxisTransform(component);
			}
			domain->addComponent(component);
		}
		domain->setepsilonRF(1e+10);
}

unsigned long DropletGenerator::readPhaseSpace(
		ParticleContainer* particleContainer,
		std::list<ChemicalPotential>* lmu, Domain* domain,
		DomainDecompBase* domainDecomp) {

	vector<double> bBoxMin;
	vector<double> bBoxMax;

	bBoxMin.resize(3);
	bBoxMax.resize(3);
	for (int i = 0; i < 3; i++) {
		bBoxMin[i] = domainDecomp->getBoundingBoxMin(i, domain);
		bBoxMax[i] = domainDecomp->getBoundingBoxMax(i, domain);
	}

	_logger->info()
			<< "OneCLJGenerator  generating cluster distribution. " << " T "
			<< _temperature << " #molecules " << numOfMolecules << " rho_gas "
			<< gasDensity << " rho_fluid " << fluidDensity << endl;
	generateMoleculesCluster(particleContainer, bBoxMin, bBoxMax, domain,
			domainDecomp);

	vector<Component>& dcomponents = domain->getComponents();
	vector<unsigned long> partsPerComp;
	partsPerComp.resize(1);
	particleContainer->update();
	particleContainer->deleteOuterParticles();
	domain->setglobalNumMolecules(
			domainDecomp->countMolecules(particleContainer, partsPerComp));

	for (unsigned int i = 0; i < partsPerComp.size(); i++) {
		dcomponents[i].setNumMolecules(partsPerComp[i]);
		domain->setglobalRotDOF(
				partsPerComp[i]
						* dcomponents[i].getRotationalDegreesOfFreedom());
	}
	domain->setglobalRho(
			domain->getglobalNumMolecules() / (simBoxLength[0]
					* simBoxLength[1] * simBoxLength[2]));
	return domain->getglobalNumMolecules();
}

void DropletGenerator::readLocalClusters(Domain* domain,
		DomainDecompBase* domainDecomp) {

	localClusters.clear();
	DropletPlacement dropletPlacement(fluidVolume, maxSphereVolume,
			numSphereSizes, _logger);
	vector<DropletPlacement::Droplet> droplets =
			dropletPlacement.generateDroplets();

	vector<double> shiftedSphere;
	vector<double> sphere;
	double distanceToDomain;
	sphere.resize(4); // x y z r
	shiftedSphere.resize(4); // x y z r
	for (unsigned int count = 0; count < droplets.size(); count++) {
		for (int i = 0; i < 3; i++) {
			sphere[i] = droplets[count]._center[i] * simBoxLength[i];
		}
		sphere[3] = droplets[count]._radius * simBoxLength[0];
		shiftedSphere[3] = sphere[3];

		for (int ix = -1; ix <= 1; ix++) {
			for (int iy = -1; iy <= 1; iy++) {
				for (int iz = -1; iz <= 1; iz++) {
					shiftedSphere[0] = sphere[0] + ix * simBoxLength[0];
					shiftedSphere[1] = sphere[1] + iy * simBoxLength[1];
					shiftedSphere[2] = sphere[2] + iz * simBoxLength[2];
					distanceToDomain = domainDecomp->guaranteedDistance(
							shiftedSphere[0], shiftedSphere[1],
							shiftedSphere[2], domain);
					if (distanceToDomain <= shiftedSphere[3]) { // sphere[3] = radius
						// reduce number of spheres for domains with periodic boundary
						bool tooFar = false;
						for (int i = 0; i < 3; i++) {
							if (shiftedSphere[i] < -sphere[3]
									|| shiftedSphere[i] > simBoxLength[i]
											+ sphere[3])
								tooFar = true;
						}
						if (not (tooFar)) {
							localClusters.push_back(shiftedSphere);
						}
					}
				}
			}
		}
	}
	//clusterStream.close();
}

bool DropletGenerator::belongsToPreviousCluster(double x, double y, double z,
		int clusterid) {
	bool belongsToOther = false;
	for (int i = 0; i < clusterid; i++) {
		if (sqrt(
				pow(x - localClusters[i][0], 2.0) + pow(
						y - localClusters[i][1], 2.0) + pow(
						z - localClusters[i][2], 2.0)) <= localClusters[i][3]) {
			belongsToOther = true;
			break;
		}
	}
	return belongsToOther;
}

bool DropletGenerator::closeToAnyCluster(double x, double y, double z,
		double offset) {
	bool belongsToOther = false;
	for (unsigned int i = 0; i < localClusters.size(); i++) {
		if (sqrt(
				pow(x - localClusters[i][0], 2.0) + pow(
						y - localClusters[i][1], 2.0) + pow(
						z - localClusters[i][2], 2.0)) <= localClusters[i][3]
				+ offset) {
			belongsToOther = true;
			break;
		}
	}
	return belongsToOther;
}

vector<ParameterCollection*> DropletGenerator::getParameters() {
	vector<ParameterCollection*> parameters;
	parameters.push_back(new MardynConfigurationParameters(_configuration));

	ParameterCollection* tab = new ParameterCollection("DropletGenerator", "DropletGenerator",
			"Parameters of DropletGenerator", Parameter::BUTTON);
	parameters.push_back(tab);

	tab->addParameter(
			new ParameterWithDoubleValue("fluidDensity", "Fluid Density: ",
					"Density in the fluid phase ", Parameter::LINE_EDIT, false, fluidDensity));
	tab->addParameter(
			new ParameterWithDoubleValue("gasDensity", "Gas Density: ",
					"Density in the gas phase", Parameter::LINE_EDIT, false, gasDensity));
	tab->addParameter(
			new ParameterWithIntValue("numOfMolecules", "numOfMolecules",
					"Number of Molecules", Parameter::LINE_EDIT, false, numOfMolecules));
	tab->addParameter(
			new ParameterWithIntValue("fluidVolume", "fluidVolume",
					"Volume of the fluid, i.e. percentage of the volume covered by fluid (X percent)",Parameter::LINE_EDIT,  false, fluidVolume));
	tab->addParameter(
			new ParameterWithDoubleValue("maxSphereVolume", "maxSphereVolume",
					"The Volume covered by the largest drop (X Percent of the total volume", Parameter::LINE_EDIT, false, maxSphereVolume));
	tab->addParameter(
			new ParameterWithIntValue("numSphereSizes", "numSphereSizes",
					"determines how many different sizes of spheres exist, with \n 1) each class covering the same volume in total and \n 2) the size of a sphere is determined by pow(0.9, i) * maxSphereRadius; i in [1,maxSphereSize]",Parameter::SPINBOX,  false, numSphereSizes));
	tab->addParameter(
			new ParameterWithDoubleValue("temperature", "Temperature [K]",
					"Temperature in the domain in Kelvin", Parameter::LINE_EDIT,
					false, _temperature / MDGenerator::kelvin_2_mardyn ));
	tab->addParameter(
			new ComponentParameters("component1", "component1", "Set up the parameters of component 1",
					_components[0]));
	return parameters;
}



void DropletGenerator::generateMoleculesCluster(
		ParticleContainer* particleContainer, vector<double> &bBoxMin,
		vector<double> &bBoxMax, Domain* domain, DomainDecompBase* domainDecomp) {

	readLocalClusters(domain, domainDecomp);

	_components[0].updateMassInertia();

	vector<int> globalFccCells;
	vector<int> clusterFccCellsMin;
	vector<int> clusterFccCellsMax;
	vector<double> fccCellLength;

	globalFccCells.resize(3);
	clusterFccCellsMin.resize(3);
	clusterFccCellsMax.resize(3);
	fccCellLength.resize(3);

	// fluid properties
	for (int dim = 0; dim < 3; dim++) {
		globalFccCells[dim] = (int) ceil(
				pow(fluidDensity / 4., 1. / 3.) * simBoxLength[dim]);
		fccCellLength[dim] = simBoxLength[dim] / globalFccCells[dim];
	}
	vector<vector<double> > fccOffsets;
	fccOffsets.resize(4);
	for (int i = 0; i < 4; i++) {
		fccOffsets[i].resize(3);
		for (int j = 0; j < 3; j++) {
			fccOffsets[i][j] = 0.0;
		}
	}
	fccOffsets[1][0] = 0.5 * fccCellLength[0];
	fccOffsets[1][1] = 0.5 * fccCellLength[1];
	fccOffsets[2][0] = 0.5 * fccCellLength[0];
	fccOffsets[2][2] = 0.5 * fccCellLength[2];
	fccOffsets[3][1] = 0.5 * fccCellLength[1];
	fccOffsets[3][2] = 0.5 * fccCellLength[2];

	double securityOffset = 0.5 * fccCellLength[0];

	vector<double> clusterPos;
	clusterPos.resize(3);
	double radius;

	double r_[3];

	int molCount = 0;
	for (unsigned int cluster = 0; cluster < localClusters.size(); cluster++) {
		radius = localClusters[cluster][3];
		for (int dim = 0; dim < 3; dim++) {
			clusterPos[dim] = localClusters[cluster][dim];
			clusterFccCellsMin[dim] = (int) floor(
					(clusterPos[dim] - radius) / fccCellLength[dim]);
			clusterFccCellsMax[dim] = (int) floor(
					(clusterPos[dim] + radius) / fccCellLength[dim]);
		}

		for (int fcc = 0; fcc < 4; fcc++) {
			for (int iz = clusterFccCellsMin[2]; iz <= clusterFccCellsMax[2]; iz++) {
				for (int iy = clusterFccCellsMin[1]; iy
						<= clusterFccCellsMax[1]; iy++) {
					for (int ix = clusterFccCellsMin[0]; ix
							<= clusterFccCellsMax[0]; ix++) {
						// Position
						r_[0] = ix * (fccCellLength[0]) + fccOffsets[fcc][0]
								+ 0.000000001;
						r_[1] = iy * (fccCellLength[1]) + fccOffsets[fcc][1]
								+ 0.000000001;
						r_[2] = iz * (fccCellLength[2]) + fccOffsets[fcc][2]
								+ 0.000000001;

						if (sqrt(
								pow(r_[0] - clusterPos[0], 2.0) + pow(
										r_[1] - clusterPos[1], 2.0) + pow(
										r_[2] - clusterPos[2], 2.0)) > radius) {
							// molecule is too far away from the center of the cluster
							continue;
						}
						// TODO put them back there
						if (not domainDecomp->procOwnsPos(r_[0], r_[1], r_[2],
								domain)) {
							// if the position is not in the domain of this proc,
							// the molecule must not be created
							continue;
						}
						if (belongsToPreviousCluster(r_[0], r_[1], r_[2],
								cluster)) {
							// some other cluster already created this molecule
							continue;
						}

						if (isInsideDomain(domain, r_)
								 && domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)) {
							double I[3] = {0.,0.,0.};
							I[0] = _components[0].I11();
							I[1] = _components[0].I22();
							I[2] = _components[0].I33();
							/*****  Loop Copied from animake - initialize anular velocity *****/
							double w[3];
							for(int d=0; d < 3; d++) {
								w[d] = (I[d] == 0)? 0.0: ((randdouble(0,1) > 0.5)? 1: -1) *
										sqrt(2.0* randdouble(0,1)* _temperature / I[d]);
								w[d] = w[d] * MDGenerator::fs_2_mardyn;
							}

							vector<double> v = getRandomVelocity(_temperature);
							Molecule m(molCount, 0, r_[0], r_[1], r_[2], v[0], v[1], v[2],
									1, 0, 0, 0, w[0], w[1], w[2], &_components);
							particleContainer->addParticle(m);
							//std::cout << "XXXXXXX DropletGenerator: added particle with ID=" << molCount << std::endl;
							molCount++;
						}

						//molCount++;
					}
				}
			}
		}
	}

	// gas properties
	vector<int> fccCellsMin;
	vector<int> fccCellsMax;
	fccCellsMin.resize(3);
	fccCellsMax.resize(3);

	for (int dim = 0; dim < 3; dim++) {
		globalFccCells[dim] = (int) ceil(
				pow(gasDensity / 4., 1. / 3.) * simBoxLength[dim]);
		fccCellLength[dim] = simBoxLength[dim] / globalFccCells[dim];
		fccCellsMin[dim] = (int) floor(bBoxMin[dim] / fccCellLength[dim]);
		fccCellsMax[dim] = (int) floor(bBoxMax[dim] / fccCellLength[dim]);
	}

	fccOffsets[1][0] = 0.5 * fccCellLength[0];
	fccOffsets[1][1] = 0.5 * fccCellLength[1];
	fccOffsets[2][0] = 0.5 * fccCellLength[0];
	fccOffsets[2][2] = 0.5 * fccCellLength[2];
	fccOffsets[3][1] = 0.5 * fccCellLength[1];
	fccOffsets[3][2] = 0.5 * fccCellLength[2];

	// GAS
	for (int fcc = 0; fcc < 4; fcc++) {
		for (int iz = fccCellsMin[2]; iz <= fccCellsMax[2]; iz++) {
			for (int iy = fccCellsMin[1]; iy <= fccCellsMax[1]; iy++) {
				for (int ix = fccCellsMin[0]; ix <= fccCellsMax[0]; ix++) {
					// Position
					r_[0] = ix * (fccCellLength[0]) + fccOffsets[fcc][0]
							+ 0.000000001;
					r_[1] = iy * (fccCellLength[1]) + fccOffsets[fcc][1]
							+ 0.000000001;
					r_[2] = iz * (fccCellLength[2]) + fccOffsets[fcc][2]
							+ 0.000000001;

					if (not domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)) {
					 // if the position is not in the domain of this proc,
					 // the molecule must not be created
					 continue;
					 }
					 if (closeToAnyCluster(r_[0], r_[1], r_[2], securityOffset)) {
					 // some other cluster already created this molecule
					 continue;
					 }

					 if (isInsideDomain(domain, r_)
							 && domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)) {
						 double I[3] = {0.,0.,0.};
						 I[0] = _components[0].I11();
						 I[1] = _components[0].I22();
						 I[2] = _components[0].I33();
						 /*****  Loop Copied from animake - initialize anular velocity *****/
						 double w[3];
						 for(int d=0; d < 3; d++) {
							 w[d] = (I[d] == 0)? 0.0: ((randdouble(0,1) > 0.5)? 1: -1) *
									 sqrt(2.0* randdouble(0,1)* _temperature / I[d]);
							 w[d] = w[d] * MDGenerator::fs_2_mardyn;
						 }

						 vector<double> v = getRandomVelocity(_temperature);
						 Molecule m(molCount, 0, r_[0], r_[1], r_[2], v[0], v[1], v[2],
								 1, 0, 0, 0, w[0], w[1], w[2], &_components);
						 particleContainer->addParticle(m);
						 //std::cout << "XXXXXXX DropletGenerator: added particle with ID=" << molCount << std::endl;
						 molCount++;
					 }
				}
			}
		}
	}
}

void DropletGenerator::setParameter(Parameter* p) {

	string id = p->getNameId();
	if (id == "fluidDensity") {
		fluidDensity = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: fluidDensity: " << fluidDensity << endl;
		setClusterParameters(gasDensity, fluidDensity, fluidVolume,
				maxSphereVolume, numSphereSizes);
	} else if (id == "gasDensity") {
		gasDensity = static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: gasDensity: " << gasDensity << endl;
		setClusterParameters(gasDensity, fluidDensity, fluidVolume,
				maxSphereVolume, numSphereSizes);
	} else if (id == "numOfMolecules") {
		numOfMolecules = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: numOfMolecules: " << numOfMolecules
				<< endl;
		setClusterParameters(gasDensity, fluidDensity, fluidVolume,
				maxSphereVolume, numSphereSizes);
	} else if (id == "fluidVolume") {
		fluidVolume = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: fluidVolume: " << fluidVolume << endl;
		setClusterParameters(gasDensity, fluidDensity, fluidVolume,
				maxSphereVolume, numSphereSizes);
	} else if (id == "maxSphereVolume") {
		maxSphereVolume
				= static_cast<ParameterWithDoubleValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: maxSphereVolume: " << maxSphereVolume
				<< endl;
		setClusterParameters(gasDensity, fluidDensity, fluidVolume,
				maxSphereVolume, numSphereSizes);
	} else if (id == "numSphereSizes") {
		numSphereSizes = static_cast<ParameterWithIntValue*> (p)->getValue();
		cout << "OneCenterLJDroplet: numSphereSizes: " << numSphereSizes
				<< endl;
		setClusterParameters(gasDensity, fluidDensity, fluidVolume,
				maxSphereVolume, numSphereSizes);
	} else if (id == "temperature") {
		_temperature = static_cast<ParameterWithDoubleValue*> (p)->getValue() * MDGenerator::kelvin_2_mardyn;
	}  else if (id.find("component1") != std::string::npos) {
		std::string part = id.substr(11);
		ComponentParameters::setParameterValue(_components[0], p, part);
	} else if (firstSubString(".", id) == "ConfigurationParameters") {
		std::string part = remainingSubString(".", id);
		MardynConfigurationParameters::setParameterValue(_configuration, p, part);
	} else {
		std::cout << "UNKOWN Parameter: id = " << p->getNameId() << " value= " << p->getStringValue() << endl;
		exit(-1);
	}
}

bool DropletGenerator::validateParameters() {
	bool valid = true;

	if (_configuration.getScenarioName() == "") {
		valid = false;
		_logger->error() << "ScenarioName not set!" << endl;
	}

	if (_configuration.getOutputFormat() == MardynConfiguration::XML) {
		valid = false;
		_logger->error() << "OutputFormat XML not yet supported!" << endl;
	}

	for (int i = 0; i < 3; i++) {
		if (simBoxLength[i] < 2.0 * _configuration.getCutoffRadius()) {
			valid = false;
			_logger->error() << "Cutoff radius is too big (there would be only 1 cell in the domain!)" << endl;
			_logger->error() << "Cutoff radius=" << _configuration.getCutoffRadius()
					<< " domain size=" << simBoxLength[i] << endl;
		}
	}

	return valid;
}

