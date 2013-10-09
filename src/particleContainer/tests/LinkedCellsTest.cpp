/*
 * LinkedCellsTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#include "LinkedCellsTest.h"

TEST_SUITE_REGISTRATION(LinkedCellsTest);

LinkedCellsTest::LinkedCellsTest() {

}

LinkedCellsTest::~LinkedCellsTest() {
}

void LinkedCellsTest::testSortParticles() {
	// loading the inp like this results in two cells per dimension with 8 molecules per dimension
	ParticleContainer* moleculeContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, "Madelung512.inp", 22.65);
	LinkedCells* lc = dynamic_cast<LinkedCells*> (moleculeContainer);
	TestCellHandler cellHandler;
	lc->traverseCells(cellHandler);
	delete lc;
}

