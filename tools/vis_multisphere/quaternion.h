/*
 * quaternion.h
 *
 *  Created on: 04.04.2011
 *      Author: eck
 */

#ifndef QUATERNION_H_
#define QUATERNION_H_

class Quaternion{
private:
	double coordinates[4];
public:
	Quaternion();
	Quaternion(double q0, double q1, double q2, double q3);
	double* gCoordinates();
	int sCoordinates(double input[4]);
	//Elementfunktion: Multiplizieren mit einem anderen Quaternion
	Quaternion operator*(Quaternion quat1); //quat1 * *this; drehen
	int rotatePoint(double point[3], double (&image)[3]);
};

#endif /* QUATERNION_H_ */
