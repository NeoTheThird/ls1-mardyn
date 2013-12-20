/*
 * Vector3d.h
 *
 *  Created on: 29.10.2013
 *      Author: mheinen
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <iostream>

class Vector3d;
class Vector3d
{
	friend std::ostream& operator<<(std::ostream&, Vector3d &);

public:
	Vector3d();
	Vector3d(double dX, double dY, double dZ);
	Vector3d(Vector3d & rhs);
	Vector3d(Vector3d& vecOrigin, Vector3d& vecEnd);
	~Vector3d();

	double GetX(){return _dX;}
	double GetY(){return _dY;}
	double GetZ(){return _dZ;}

	Vector3d operator+(Vector3d & rhs);
	Vector3d operator-(Vector3d & rhs);
	double operator*(Vector3d & rhs);
	Vector3d operator*(double dFac);

	Vector3d CrossProduct(Vector3d & vec);
	Vector3d Normalize();
	double GetLength();
	void CalcRadialTangentialComponentLength( Vector3d &vec3dMidpoint, Vector3d vec3dPosition,
			                                  double &dVelocity_r, double &dVelocity_t1, double &dVelocity_t2);

private:
	double _dX;
	double _dY;
	double _dZ;

};  // class Vector3d


#endif /* VECTOR3D_H_ */
