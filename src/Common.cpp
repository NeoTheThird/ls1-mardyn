// Common.cpp

#include "Common.h"

#include <sstream>
#include <ctime>
#include <cmath>
#include <iostream>

#include <xmmintrin.h>

using namespace std;

string gettimestring(const char* fmt) {
	time_t rawtime;
	struct tm* timeinfo;
	char buffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, fmt, timeinfo);
	string date(buffer);

	return date;
}

/**
 * Align a number to the right by padding with a filling character.
 *
 * @param number the number value to be aligned
 * @param num_digits the number of digits
 * @param c character used for filling up
 */
string aligned_number(int number, int num_digits, char c) {
	stringstream numstream;
	numstream.fill(c);
	numstream.width(num_digits);
	numstream << number;
	string numstr(numstream.str());
	return numstr;
}


void calculateDistances( float * __restrict__ valuesA[3], float* __restrict__ const valuesB[3], int numValuesA, int numValuesB,
		float ** __restrict__ distances, float *** __restrict__ distanceVectors ) {

/*	for (int i = 0; i < numValuesA; i++) {
		std::cout << "A Points: " << valuesA[0][i] << "," << valuesA[1][i] << "," << valuesA[2][i] << endl;
	}
	for (int i = 0; i < numValuesB; i++) {
			std::cout << "B Points: " << valuesB[0][i] << "," << valuesB[1][i] << "," << valuesB[2][i] << endl;
	}
*/

	for (int i = 0; i < numValuesA; i++) {
		float xA = valuesA[0][i];
		float yA = valuesA[1][i];
		float zA = valuesA[2][i];
#ifdef __INTEL_COMPILER
		#pragma ivdep
#endif
		for (int j = 0; j < numValuesB; j++) {
			distanceVectors[0][i][j] = xA - valuesB[0][j];
			distanceVectors[1][i][j] = yA - valuesB[1][j];
			distanceVectors[2][i][j] = zA - valuesB[2][j];
		}
	}

	for (int i = 0; i < numValuesA; i++) {
#ifdef __INTEL_COMPILER
		#pragma ivdep
#endif
		for (int j = 0; j < numValuesB; j++) {
			distances[i][j] = distanceVectors[0][i][j] * distanceVectors[0][i][j] +
			                  distanceVectors[1][i][j] * distanceVectors[1][i][j] +
			                  distanceVectors[2][i][j] * distanceVectors[2][i][j];
		}
	}
}


/**
 * @param valuesA positions of molecules in set A: [x1, y1, z1, x2, y2, z2, ... , xN, yN, zN]; Length = 3 * numValuesA
 * @param valuesB positions of molecules in set B
 * @param distanceVectors: linearized array, innermost x,y,z, then molecules b; length = numValuesA * numValuesB * 3
 *        access distance [a][b][0] as a * numValuesB * 3 + b + 0
 * @param distances length = numValuesA * numValuesB; access distance [a][b] as a * numValuesB + b
 */
/*void calculateDistances( double valuesA[], double const valuesB[], int numValuesA, int numValuesB,
		double * __restrict__ distances, double * __restrict__ distanceVectors ) {

	for (int i = 0; i < numValuesA; i++) {
		double xA = valuesA[i];
		double yA = valuesA[i+1];
		double zA = valuesA[i+2];
		for (int j = 0; j < numValuesB; j+=3) {
			distanceVectors[i * numValuesB + j ] = xA - valuesB[j];
			distanceVectors[i * numValuesB + j + 1] = yA - valuesB[j + 1];
			distanceVectors[i * numValuesB + j + 2] = zA - valuesB[j + 2];
		}
	}

	for (int i = 0; i < numValuesA * numValuesB; i+=3) {
		distances[i] = distanceVectors[i] + distanceVectors[i];
		distances[i+1] = distanceVectors[i+1] + distanceVectors[i+1];
		distances[i+2] = distanceVectors[i+2] + distanceVectors[i+2];
	}
}*/


void calculateDistances( int numA, float* A[3], int numB, float* B[3], float** distances) {

	//int numIterationsA = numA / 4;
	int numIterationsB = numB / 4 + (numB % 4 == 0 ? 0 : 1) ;

	//std::cout << "SSE Distances: numA=" << numA << " numIterationsB=" << numIterationsB << std::endl;
	for (int i = 0; i < numA; i++) {
		__m128 xCoordsA = _mm_load1_ps(A[0] + i);
		__m128 yCoordsA = _mm_load1_ps(A[1] + i);
		__m128 zCoordsA = _mm_load1_ps(A[2] + i);

		for (int j = 0; j < numIterationsB; j++) {


			__m128 xCoordsB = _mm_load_ps(B[0] + j * 4);
			xCoordsB = _mm_sub_ps(xCoordsA, xCoordsB);
			__m128 yCoordsB = _mm_load_ps(B[1] + j * 4);
			yCoordsB = _mm_sub_ps(yCoordsA, yCoordsB);
			__m128 zCoordsB = _mm_load_ps(B[2] + j * 4);
			zCoordsB = _mm_sub_ps(zCoordsA, zCoordsB);

			xCoordsB = _mm_mul_ps(xCoordsB, xCoordsB);
			yCoordsB = _mm_mul_ps(yCoordsB, yCoordsB);
			zCoordsB = _mm_mul_ps(zCoordsB, zCoordsB);

			xCoordsB = _mm_add_ps(xCoordsB, yCoordsB);
			xCoordsB = _mm_add_ps(xCoordsB, zCoordsB);

			_mm_store_ps(distances[i] + j*4, xCoordsB);
		}
	}

}


void calculateInteractionTable( int numA, float* A[3], int numB, float* B[3], float** distances, float cutoff) {

	__m128 cutoffRegister = _mm_load1_ps(&cutoff);

	//int numIterationsA = numA / 4;
	int numIterationsB = numB / 4 + (numB % 4 == 0 ? 0 : 1) ;

	//std::cout << "SSE Distances: numA=" << numA << " numIterationsB=" << numIterationsB << std::endl;
	for (int i = 0; i < numA; i++) {
		__m128 xCoordsA = _mm_load1_ps(A[0] + i);
		__m128 yCoordsA = _mm_load1_ps(A[1] + i);
		__m128 zCoordsA = _mm_load1_ps(A[2] + i);

		for (int j = 0; j < numIterationsB; j++) {


			__m128 xCoordsB = _mm_load_ps(B[0] + j * 4);
			xCoordsB = _mm_sub_ps(xCoordsA, xCoordsB);
			__m128 yCoordsB = _mm_load_ps(B[1] + j * 4);
			yCoordsB = _mm_sub_ps(yCoordsA, yCoordsB);
			__m128 zCoordsB = _mm_load_ps(B[2] + j * 4);
			zCoordsB = _mm_sub_ps(zCoordsA, zCoordsB);

			xCoordsB = _mm_mul_ps(xCoordsB, xCoordsB);
			yCoordsB = _mm_mul_ps(yCoordsB, yCoordsB);
			zCoordsB = _mm_mul_ps(zCoordsB, zCoordsB);

			xCoordsB = _mm_add_ps(xCoordsB, yCoordsB);
			xCoordsB = _mm_add_ps(xCoordsB, zCoordsB);

			xCoordsB = _mm_cmplt_ps(xCoordsB, cutoffRegister);
			_mm_store_ps(distances[i] + j*4, xCoordsB);
		}
	}

}
