/*
 * ReorderedLinkedCellsTest.h
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#ifndef REORDEREDLINKEDCELLSTEST_H_
#define REORDEREDLINKEDCELLSTEST_H_

#include "particleContainer/tests/ParticleContainerTest.h"
#include "particleContainer/ReorderedLinkedCells.h"

class ReorderedLinkedCellsTest: public ParticleContainerTest {

	TEST_SUITE(ReorderedLinkedCellsTest);
	TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_SUITE_END();

public:

	ReorderedLinkedCellsTest();

	virtual ~ReorderedLinkedCellsTest();

	void testInsertion() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		ReorderedLinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 2.5, 1, NULL);
		this->ParticleContainerTest::testInsertion(&container);
	}

	void testMoleculeIteration() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		ReorderedLinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 2.5, 1, NULL);
		this->ParticleContainerTest::testMoleculeIteration(&container);
	}

	void testUpdateAndDeleteOuterParticles() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		ReorderedLinkedCells container(boundings_min, boundings_max, 2.5, 2.5, 2.5, 1, NULL);
		this->ParticleContainerTest::testUpdateAndDeleteOuterParticles(&container);
	}
};

#endif /* REORDEREDLINKEDCELLSTEST_H_ */
