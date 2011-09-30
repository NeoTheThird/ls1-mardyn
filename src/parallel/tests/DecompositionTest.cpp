/*
 * DecompositionTest.cpp
 *
 * @Date: 30.09.2011
 * @Author: eckhardw
 */

#include "DecompositionTest.h"
#include "parallel/DomainDecomposition.h"
#include "Domain.h"
#include "molecules/MoleculeTypes.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/BlockedReorderedLinkedCells.h"

using namespace std;

TEST_SUITE_REGISTRATION(DecompositionTest);

DecompositionTest::DecompositionTest() {
}

DecompositionTest::~DecompositionTest() {
}


void DecompositionTest::testExchangeBasicMolecules() {

#ifndef ENABLE_MPI
	Log::global_log->warning() << "SKIPPED DecompositionTest::testExchangeBasicMolecules()!" << endl;
	Log::global_log->warning() << "For running this test, compile with \"-DENABLE_MPI\"!" << endl;
	return;
#endif
	if (IsSame<Molecule, BasicMolecule>::Result::value == false) {
		Log::global_log->warning() << "SKIPPED DecompositionTest::testExchangeBasicMolecules()!" << endl;
		Log::global_log->warning() << "For running this test, \"MoleculeType\" has to be defined as "
				"\"BasicMolecule\" in molecules/MoleculeTypes.h!" << endl;
		return;
	}

	// initialize the two linked cells

	LinkedCells* linkedCells = static_cast<LinkedCells*> (
			initializeFromFile(ParticleContainerFactory::LinkedCell, "DecompositionTest.inp", 13));
	BlockedReorderedLinkedCells* blockedCells = static_cast<BlockedReorderedLinkedCells*>(
			initializeFromFile(ParticleContainerFactory::BlockedReorderedLinkedCell, "DecompositionTest.inp", 13));

	AssertEqualParticleContainers(linkedCells, blockedCells);

	DomainDecomposition* domainDecomposition = dynamic_cast<DomainDecomposition*>(_domainDecomposition);
	if (domainDecomposition == NULL) {
		ASSERT_FAIL("Error: _domainDecomposition has to be of type DomainDecomposition in order to test the exchange of basic molecules!");
	}

	// perform the exchange
	std::cout << "LinkedCells.size() = " << linkedCells->getNumberOfParticles() << endl;
	std::cout << "BlockedCells.size() = " << blockedCells->getNumberOfParticles() << endl;
	domainDecomposition->exchangeParticleData(linkedCells, _domain->getComponents(), _domain);
	domainDecomposition->exchangeBasicMolecules(blockedCells, _domain->getComponents(), _domain);
	std::cout << "BlockedCells.size() = " << blockedCells->getNumberOfParticles() << endl;
	std::cout << "LinkedCells.size() = " << linkedCells->getNumberOfParticles() << endl;

	// make sure their content is the same
	AssertEqualParticleContainers(linkedCells, blockedCells);
}


void DecompositionTest::AssertEqualParticleContainers(ParticleContainer* first, ParticleContainer* second) {
	ASSERT_EQUAL(first->getNumberOfParticles(), second->getNumberOfParticles());
	Molecule* lcMolecule = first->begin();
	while(lcMolecule != first->end()) {
		Molecule* blockedMolecule = second->begin();
		bool found = false;

		while (blockedMolecule != second->end()) {
			if (blockedMolecule->id() == lcMolecule->id()) {
				found = true;
				// we can't compare the positions as they may be changed du to the periodic boundaries
				ASSERT_DOUBLES_EQUAL(lcMolecule->v(0), blockedMolecule->v(0), 0.0000001);
				ASSERT_DOUBLES_EQUAL(lcMolecule->v(1), blockedMolecule->v(1), 0.0000001);
				ASSERT_DOUBLES_EQUAL(lcMolecule->v(2), blockedMolecule->v(2), 0.0000001);
				ASSERT_DOUBLES_EQUAL(lcMolecule->D(0), blockedMolecule->D(0), 0.0000001);
				ASSERT_DOUBLES_EQUAL(lcMolecule->D(1), blockedMolecule->D(1), 0.0000001);
				ASSERT_DOUBLES_EQUAL(lcMolecule->D(2), blockedMolecule->D(2), 0.0000001);
				break;
			}
			blockedMolecule = second->next();
		}
		if (!found) {
			ASSERT_FAIL("Molecule not found in Container Second! Something is wrong in initialisation of test.");
		}

		lcMolecule = first->next();
	}
}
