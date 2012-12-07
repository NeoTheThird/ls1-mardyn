/*
 * molecule.h
 *
 *  Created on: 10.11.2010
 *      Author: eck
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <ostream>
#include "quaternion.h"

class Molecule{
private:
	int idNumber;
	int componentType;
	double position[3];
	Quaternion orientationQuaternion; //Orientierung in Quaternionen (q0, q1, q2, q3) mit q = q0+q1*i+q2*j+q3*k
public:
	int gIdNumber() const;
	int sIdNumber(int input);
	int gComponentType() const;
	int sComponentType(int input);
	double* gPosition();
	int sPosition(double input[3]);
	Quaternion gOrientationQuaternion() const;
	int sOrientationQuaternion(Quaternion quat);
	int sOrientationQuaternion(double input[4]);
};

#endif /* MOLECULE_H_ */
