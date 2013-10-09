/*
 * LinkedCellsTest.h
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#ifndef LINKEDCELLSTEST_H_
#define LINKEDCELLSTEST_H_

#include "particleContainer/tests/ParticleContainerTest.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "molecules/Molecule.h"

#include <vector>

using namespace std;

class LinkedCellsTest: public ParticleContainerTest {

	TEST_SUITE(LinkedCellsTest);
	TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_METHOD(testSortParticles);
	TEST_SUITE_END();

public:

	LinkedCellsTest();

	virtual ~LinkedCellsTest();

	void testInsertion() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		LinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 1);
		this->ParticleContainerTest::testInsertion(&container);
	}

	void testMoleculeIteration() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		LinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 1);
		this->ParticleContainerTest::testMoleculeIteration(&container);
	}

	void testUpdateAndDeleteOuterParticles() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		LinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 1);
		this->ParticleContainerTest::testUpdateAndDeleteOuterParticles(&container);
	}

	void testSortParticles();


private:

	class TestCellHandler : public CellProcessor {

		virtual void initTraversal(const size_t numCells) { }

		virtual void preprocessCell(ParticleCell& cell) { }

		virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2) { }

		virtual void processCell(ParticleCell& cell) {
			vector<Molecule*>& molecules = cell.getParticlePointers();

			if (molecules.size() > 0) {
				vector<Molecule*>::iterator mol = molecules.begin();

				int nBuckets = (NBUCKET * NBUCKET);
				int numMoleculesPerBucket = 64 / nBuckets;
				for (int i = 0; i < nBuckets; i++) {
					double current_z = (*mol)->r(2);
					for (int j = 0; j < numMoleculesPerBucket-1; j++) {
						mol++;
						ASSERT_TRUE(current_z <= (*mol)->r(2));
						current_z = (*mol)->r(2);
					}
					mol++;
				}
				ASSERT_TRUE_MSG("Iterator not at end!", (mol == molecules.end()));
			}
		}

		virtual void postprocessCell(ParticleCell& cell) { }

		virtual void endTraversal() { }
	};

};

#endif /* LINKEDCELLSTEST_H_ */
