/*
 * DecompositionTest.h
 *
 * @Date: 30.09.2011
 * @Author: eckhardw
 */

#ifndef DECOMPOSITIONTEST_H_
#define DECOMPOSITIONTEST_H_

#include "utils/TestWithSimulationSetup.h"

class DecompositionTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(DecompositionTest);
	TEST_METHOD(testExchangeBasicMolecules);
	TEST_SUITE_END();

public:

	DecompositionTest();

	virtual ~DecompositionTest();

	/**
	 * test the particle exchange of BasicMolecules vs. the reference implementation
	 * "exchangeParticleData".
	 */
	void testExchangeBasicMolecules();

	/**
	 * utility method to make sure that the molecules contained are equal
	 */
	void AssertEqualParticleContainers(ParticleContainer* first, ParticleContainer* second);
};

#endif /* DECOMPOSITIONTEST_H_ */
