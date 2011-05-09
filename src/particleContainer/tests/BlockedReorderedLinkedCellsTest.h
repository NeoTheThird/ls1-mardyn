/*
 * BlockedReorderedLinkedCellsTest.h
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#ifndef BLOCKEDREORDEREDLINKEDCELLSTEST_H_
#define BLOCKEDREORDEREDLINKEDCELLSTEST_H_

#include "particleContainer/tests/ParticleContainerTest.h"
#include "particleContainer/BlockedReorderedLinkedCells.h"

class BlockedReorderedLinkedCellsTest: public ParticleContainerTest {

	TEST_SUITE(BlockedReorderedLinkedCellsTest);
	TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_SUITE_END();

public:

	BlockedReorderedLinkedCellsTest();

	virtual ~BlockedReorderedLinkedCellsTest();

	void testInsertion() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		BlockedReorderedLinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 2.5, 1, NULL);
		this->ParticleContainerTest::testInsertion(&container);
	}

	void testMoleculeIteration() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		BlockedReorderedLinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 2.5, 1, NULL);
		this->ParticleContainerTest::testMoleculeIteration(&container);
	}

	void testUpdateAndDeleteOuterParticles() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		BlockedReorderedLinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 2.5, 1, NULL);
		this->ParticleContainerTest::testUpdateAndDeleteOuterParticles(&container);
	}
};

#endif /* BLOCKEDREORDEREDLINKEDCELLSTEST_H_ */
