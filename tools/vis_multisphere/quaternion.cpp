/*
 * quaternion.cpp
 *
 *  Created on: 04.04.2011
 *      Author: eck
 */

using namespace std;

#include "quaternion.h"

Quaternion::Quaternion(){
	//Drehung um 0° um die x-Achse
	double input[4] = {0, 1, 0, 0};
	sCoordinates(input);
}

Quaternion::Quaternion(double q0, double q1, double q2, double q3){
	double input[4] = {q0, q1, q2, q3};
	sCoordinates(input);
}

double* Quaternion::gCoordinates(){
	return coordinates;
}

int Quaternion::sCoordinates(double input[4]){
	for (int i=0; i<=3; i++)
		coordinates[i] = input[i];
	return 0;
}

//Elementfunktion: überladener *Operator
Quaternion Quaternion::operator*(Quaternion quat1){
	double* c1 = quat1.gCoordinates();
	double* c2 = coordinates;
	double k[4];
	k[0] = c1[0]*c2[0]-c1[1]*c2[1]-c1[2]*c2[2]-c1[3]*c2[3];
	k[1] = c1[0]*c2[1]+c1[1]*c2[0]+c1[2]*c2[3]-c1[3]*c2[2];
	k[2] = c1[0]*c2[2]+c1[2]*c2[0]-c1[1]*c2[3]+c1[3]*c2[1];
	k[3] = c1[0]*c2[3]+c1[3]*c2[0]+c1[1]*c2[2]-c1[2]*c2[1];
	sCoordinates(k);
	return *this;
}

int Quaternion::rotatePoint(double point[3], double (&image)[3]){
	double* c = coordinates;
	//Rotationsmatrix D=[D11 D12 D13; D21 D22 D23; D31 D32 D33]
	double D11, D12, D13, D21, D22, D23, D31, D32, D33;
	D11 = 1-2*(c[2]*c[2]+c[3]*c[3]);
	D12 = -2*c[0]*c[3]+2*c[1]*c[2];
	D13 = 2*c[0]*c[2]+2*c[1]*c[3];
	D21 = 2*c[0]*c[3]+2*c[1]*c[2];
	D22 = 1-2*(c[1]*c[1]+c[3]*c[3]);
	D23 = -2*c[0]*c[1]+2*c[2]*c[3];
	D31 = -2*c[0]*c[2]+2*c[1]*c[3];
	D32 = 2*c[0]*c[1]+2*c[2]*c[3];
	D33 = 1-2*(c[1]*c[1]+c[2]*c[2]);
	image[0] = D11*point[0]+D12*point[1]+D13*point[2];
	image[1] = D21*point[0]+D22*point[1]+D23*point[2];
	image[2] = D31*point[0]+D32*point[1]+D33*point[2];
	return 0;
}
