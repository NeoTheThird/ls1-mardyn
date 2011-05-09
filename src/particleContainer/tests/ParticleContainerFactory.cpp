/*
 * ParticleContainerFactory.cpp
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#include "particleContainer/tests/ParticleContainerFactory.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/ReorderedLinkedCells.h"
#include "particleContainer/BlockedReorderedLinkedCells.h"
#include "particleContainer/AdaptiveSubCells.h"
#include "particleContainer/ReorderedLinkedCells.h"

#include "ensemble/GrandCanonical.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

#include "io/InputOldstyle.h"
#include "utils/Logger.h"

#include <list>

using namespace Log;

ParticleContainer* ParticleContainerFactory::createEmptyParticleContainer(type type) {
	double bBoxMin[] = {0.0, 0.0, 0.0, 0.0};
	double bBoxMax[] = {2.0, 2.0, 2.0, 2.0};
	double cutoffRadius = 1.0;
	double LJCutoffRadius = 1.0;
	double tersoffCutoffRadius = 1.0;
	double cellsInCutoffRadius = 1.0;

	if (type == LinkedCell) {
		LinkedCells* container = new LinkedCells(bBoxMin, bBoxMax, cutoffRadius, LJCutoffRadius,
		                                        tersoffCutoffRadius, cellsInCutoffRadius, NULL);
		return container;
	} else if (type == ReorderedLinkedCell) {
		ReorderedLinkedCells* container = new ReorderedLinkedCells(bBoxMin, bBoxMax, cutoffRadius, LJCutoffRadius,
				                                        tersoffCutoffRadius, cellsInCutoffRadius, NULL);
				return container;
	} else if (type == BlockedReorderedLinkedCell) {
		BlockedReorderedLinkedCells* container = new BlockedReorderedLinkedCells(bBoxMin, bBoxMax, cutoffRadius, LJCutoffRadius,
				                                        tersoffCutoffRadius, cellsInCutoffRadius, NULL);
				return container;
	} else {
		global_log->error() << "ParticleContainerFactory: Unsupported type requested! " << std::endl;
		return NULL;
	}
}



ParticleContainer* ParticleContainerFactory::createInitializedParticleContainer(
		type type, Domain* domain, DomainDecompBase* domainDecomposition, double cutoff, const std::string& fileName) {

	InputOldstyle inputReader;
	inputReader.setPhaseSpaceHeaderFile(fileName.c_str());
	inputReader.setPhaseSpaceFile(fileName.c_str());
	inputReader.readPhaseSpaceHeader(domain, 1.0);

	double bBoxMin[3];
	double bBoxMax[3];
	for (int i = 0; i < 3; i++) {
		bBoxMin[i] = domainDecomposition->getBoundingBoxMin(i, domain);
		bBoxMax[i] = domainDecomposition->getBoundingBoxMax(i, domain);
	}

	ParticleContainer* moleculeContainer;
	if (type == LinkedCell) {
		moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, cutoff, cutoff, cutoff, 1.0, NULL);
	} else if (type == ReorderedLinkedCell) {
		moleculeContainer = new ReorderedLinkedCells(bBoxMin, bBoxMax, cutoff, cutoff, cutoff, 1.0, NULL);
	} else if (type == BlockedReorderedLinkedCell) {
		moleculeContainer = new BlockedReorderedLinkedCells(bBoxMin, bBoxMax, cutoff, cutoff, cutoff, 1.0, NULL);
	} else if (type == AdaptiveSubCell) {
		global_log->error() << "ParticleContainerFactory: Implement ParticleContainerFactory for AdaptiveSubCells! " << std::endl;
		//moleculeContainer = new ReorderedLinkedCells(bBoxMin, bBoxMax, cutoff, cutoff, cutoff, 1.0, NULL);
		return NULL;
	} else {
		global_log->error() << "ParticleContainerFactory: Unsupported type requested! " << std::endl;
		return NULL;
	}

	std::list<ChemicalPotential> chemPot;
	inputReader.readPhaseSpace(moleculeContainer, &chemPot, domain, domainDecomposition);
	moleculeContainer->update();
	moleculeContainer->deleteOuterParticles();

	domain->initParameterStreams(cutoff, cutoff);
	return moleculeContainer;
}
