/*
 * molecule.cpp
 *
 *  Created on: 10.11.2010
 *      Author: eck
 */

using namespace std;

#include "molecule.h"

int Molecule::gIdNumber() const{
	return idNumber;
}

int Molecule::sIdNumber(int input){
	idNumber = input;
	return 0;
}

int Molecule::gComponentType() const{
	return componentType;
}

int Molecule::sComponentType(int input){
	componentType = input;
	return 0;
}

double* Molecule::gPosition(){
	return position;
}

int Molecule::sPosition(double input[3]){
	for (int i=0; i<=2; i++)
		position[i] = input[i];
	return 0;
}

Quaternion Molecule::gOrientationQuaternion() const{
	return orientationQuaternion;
}

int Molecule::sOrientationQuaternion(Quaternion quat){
	orientationQuaternion = quat;
	return 0;
}

int Molecule::sOrientationQuaternion(double input[4]){
	orientationQuaternion.sCoordinates(input);
	return 0;
}
