/*
 * BlockedReorderedLinkedCellsTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#include "particleContainer/tests/BlockedReorderedLinkedCellsTest.h"

TEST_SUITE_REGISTRATION(BlockedReorderedLinkedCellsTest);

BlockedReorderedLinkedCellsTest::BlockedReorderedLinkedCellsTest() {

}

BlockedReorderedLinkedCellsTest::~BlockedReorderedLinkedCellsTest() {
}


void BlockedReorderedLinkedCellsTest::testLinearize() {
	BlockedReorderedLinkedCells* container = static_cast<BlockedReorderedLinkedCells*>(initializeFromFile(ParticleContainerFactory::BlockedReorderedLinkedCell, "1clj-regular-12x12x12.inp", 1.8));
	utils::DynamicArray<Molecule, true, false> molecules;
	std::vector<int> cellStartIndices;
	container->linearize(molecules, cellStartIndices);

	// 1728 molecules
	ASSERT_EQUAL((size_t) 1728, molecules.size());
	// (1 + 6 + 1)^3 = 512 cells
	ASSERT_EQUAL((size_t) 512, cellStartIndices.size());

	for (int i = 0; i < 74; i++) {
		ASSERT_EQUAL(0, cellStartIndices[0]);
	}

	// the first cell (index=73) contains 8 molecules, thus the cellStartIndex
	// of the second cell has to be 8
	ASSERT_EQUAL(8, cellStartIndices[74]);

	for (int i = 0; i < 1728; i++) {
		if (molecules[i].id() == 170) {
			double moment[] = {5.5, 6.6, 7.7};
			molecules[i].Madd(moment);
		}
	}

	container->delinearize(molecules);

	for (Molecule* tmp = container->begin(); tmp != container->end(); tmp = container->next()) {
		if (tmp->id() == 170) {
			ASSERT_DOUBLES_EQUAL(5.5, tmp->M(0), 10e-6);
			ASSERT_DOUBLES_EQUAL(6.6, tmp->M(1), 10e-6);
			ASSERT_DOUBLES_EQUAL(7.7, tmp->M(2), 10e-6);
		}
	}
}
