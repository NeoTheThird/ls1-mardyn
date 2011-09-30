/*
 * MoleculeTest.h
 *
 * @Date: 04.02.2011
 * @Author: eckhardw
 */

#ifndef MOLECULETEST_H_
#define MOLECULETEST_H_

#include "utils/Testing.h"

class MoleculeTest : public utils::Test {

	TEST_SUITE(MoleculeTest);
	TEST_METHOD(testIsLessThan);
	TEST_METHOD(testRotationZero);
	TEST_METHOD(testRotationPI);
	TEST_METHOD(testForceMoment);
	TEST_SUITE_END();

public:

	MoleculeTest();

	virtual ~MoleculeTest();

	void testIsLessThan();

	void testForceMoment();

	// test the calculation of the site coordinates, if the orientation of the
	// molecule is the same like the global coordinate system (i.e. rotation is angle zero)
	void testRotationZero();

	// test a rotation of the positions around the z-axis
	// through an angle of pi.
	void testRotationPI();
};

#endif /* MOLECULETEST_H_ */
