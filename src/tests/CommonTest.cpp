/*
 * CommonTest.cpp
 *
 * @Date: 21.05.2010
 * @Author: Wolfgang Eckhardt
 */

#include "CommonTest.h"
#include "Common.h"
#include "utils/Logger.h"

#include <xmmintrin.h>
#include <malloc.h>

#include <iostream>

TEST_SUITE_REGISTRATION(CommonTest);

CommonTest::CommonTest() {
}

CommonTest::~CommonTest() {
}

void CommonTest::testGetTimeString() {
	std::string time = gettimestring();
	ASSERT_EQUAL((size_t)13, time.size());
	ASSERT_EQUAL(time[6], 'T');
}


void CommonTest::testAlignedNumber() {
//	ASSERT_EQUAL_MSG("One should be zero!", 1, 0);

	std::string result = aligned_number(123, 7, '.');
	std::string expected("....123");
	ASSERT_EQUAL_MSG("Align number 123 to 7 digits, filling char is .", expected, result);

	result = aligned_number(-2, 6, ' ');
	expected = "    -2";
	ASSERT_EQUAL_MSG("Align number -2 to 6 digits, filling char is ' '(space)", expected, result);
}

void CommonTest::testCalculateDistances() {
	typedef float fp_type;

	fp_type* valuesA[3];
	fp_type* valuesB[3];
	for (int i = 0; i < 3; i++) {
		valuesA[i] = new fp_type[3];
		valuesB[i] = new fp_type[1];
		valuesB[i][0] = 4.;
		for (int j = 0; j < 3; j++) {
			valuesA[i][j] = j+1;
		}
	}

	/*for (int i = 0; i < 3; i++) {
		std::cout << "A Points: " << valuesA[0][i] << "," << valuesA[1][i] << "," << valuesA[2][i] << std::endl;
	}
	for (int i = 0; i < 1; i++) {
			std::cout << "B Points: " << valuesB[0][i] << "," << valuesB[1][i] << "," << valuesB[2][i] << std::endl;
	}*/


	fp_type** distances;
	distances = new fp_type*[3];
	for (int i = 0; i < 3; i++) {
		distances[i] = new fp_type[1];
	}
	fp_type** distanceVectors[3];
	for (int i = 0; i < 3; i++) {
		distanceVectors[i] = new fp_type*[3];
		for (int j = 0; j < 3; j++) {
			distanceVectors[i][j] = new fp_type[1];
		}
	}

	calculateDistances(valuesA, valuesB, 3, 1, distances, distanceVectors);

	ASSERT_DOUBLES_EQUAL(distanceVectors[0][0][0], -3, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[1][0][0], -3, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[2][0][0], -3, 0.000001);
	ASSERT_DOUBLES_EQUAL(distances[0][0], 27, 0.000001);

	ASSERT_DOUBLES_EQUAL(distanceVectors[0][1][0], -2, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[1][1][0], -2, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[2][1][0], -2, 0.000001);
	ASSERT_DOUBLES_EQUAL(distances[1][0], 12, 0.000001);

	ASSERT_DOUBLES_EQUAL(distanceVectors[0][2][0], -1, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[1][2][0], -1, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[2][2][0], -1, 0.000001);
	ASSERT_DOUBLES_EQUAL(distances[2][0], 3, 0.000001);
}

void CommonTest::testCalculateDistancesFloat() {
#ifdef __INTEL_COMPILER
	__declspec( align(16) ) float x[4]; // = {1, 2, 2, 2};
#else
	float x[4] __attribute__ ((aligned (16))) = {1, 2, 3, 4};
	float y[4] __attribute__ ((aligned (16))) = {1, 2, 3, 4};
	float z[4] __attribute__ ((aligned (16))) = {1, 2, 3, 4};

	float x2[5] __attribute__ ((aligned (16))) = {1, 2, 3, 4, 5};
	float y2[5] __attribute__ ((aligned (16))) = {1, 2, 3, 4, 5};
	float z2[5] __attribute__ ((aligned (16))) = {1, 2, 3, 4, 5};

	float** distances = new float*[4];
	for (int i = 0; i < 4; i++) {
		distances[i] = static_cast<float*>(_mm_malloc(8 * sizeof(float), 16));
	}

	float* coordsA[3];
	coordsA[0] = x;
	coordsA[1] = y;
	coordsA[2] = z;

	float* coordsB[3];
	coordsB[0] = x2;
	coordsB[1] = y2;
	coordsB[2] = z2;

	calculateDistances(4, coordsA, 5, coordsB, distances);

	double expectedDistances[4][5] = { {0, 3, 12, 27, 48},
	                                   {3, 0, 3, 12, 27},
	                                   {12, 3, 0, 3, 12},
	                                   {27, 12, 3, 0, 3}};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++) {
			std::stringstream msg;
			msg << "expectedDistance[" << i << "]["<<j<<"]=" << expectedDistances[i][j] << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(msg.str().c_str(), expectedDistances[i][j], distances[i][j], 0.000001);
		}
	}

	for (int i = 0; i < 4; i++) {
		_mm_free(distances[i]);
	}

#endif
}
