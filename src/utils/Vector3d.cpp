/*
 * Vector3d.cpp
 *
 *  Created on: 29.10.2013
 *      Author: mheinen
 */

#include "Vector3d.h"
#include <iostream>
#include <cmath>

using namespace std;

Vector3d::Vector3d()
{
	_dX = 0.0;
	_dY = 0.0;
	_dZ = 0.0;
}

Vector3d::Vector3d(double dX, double dY, double dZ)
{
	_dX = dX;
	_dY = dY;
	_dZ = dZ;
}

Vector3d::Vector3d(Vector3d & rhs)
{
	_dX = rhs.GetX();
	_dY = rhs.GetY();
	_dZ = rhs.GetZ();
}

Vector3d::Vector3d(Vector3d& vecOrigin, Vector3d& vecEnd)
{
	_dX = vecEnd.GetX() - vecOrigin.GetX();
	_dY = vecEnd.GetY() - vecOrigin.GetY();
	_dZ = vecEnd.GetZ() - vecOrigin.GetZ();
}

Vector3d::~Vector3d()
{
}



Vector3d Vector3d::operator+(Vector3d & rhs)
{
	double dX, dY, dZ;

	dX = _dX + rhs.GetX();
	dY = _dY + rhs.GetY();
	dZ = _dZ + rhs.GetZ();

	Vector3d vec(dX, dY, dZ);

	return vec;
}


Vector3d Vector3d::operator- (Vector3d & rhs)
{
	double dX, dY, dZ;

	dX = _dX - rhs.GetX();
	dY = _dY - rhs.GetY();
	dZ = _dZ - rhs.GetZ();

	Vector3d vec(dX, dY, dZ);

	return vec;
}

double Vector3d::operator* (Vector3d & rhs)
{
	double dScalar;

	dScalar = _dX * rhs.GetX() + _dY * rhs.GetY() + _dZ * rhs.GetZ();

	return dScalar;
}


Vector3d Vector3d::operator*(double dFac)
{
	double dX, dY, dZ;

	dX = _dX * dFac;
	dY = _dY * dFac;
	dZ = _dZ * dFac;

	Vector3d res(dX, dY, dZ);

	return res;
}


Vector3d Vector3d::CrossProduct(Vector3d & vec)
{
	double dX, dY, dZ;

	dX = _dY*vec.GetZ() - _dZ*vec.GetY();
	dY = _dZ*vec.GetX() - _dX*vec.GetZ();
	dZ = _dX*vec.GetY() - _dY*vec.GetX();

	Vector3d res(dX, dY, dZ);

	return res;
}

Vector3d Vector3d::Normalize()
{
	double dX, dY, dZ;

	dX = _dX / this->GetLength();
	dY = _dY / this->GetLength();
	dZ = _dZ / this->GetLength();

	Vector3d res(dX, dY, dZ);

	return res;
}

double Vector3d::GetLength()
{
	double dLength;

	dLength = sqrt( _dX*_dX + _dY*_dY + _dZ*_dZ );

	return dLength;
}

void Vector3d::CalcRadialTangentialComponentLength( Vector3d &vec3dMidpoint, Vector3d vec3dPosition,
		                                            double &dVelocity_r, double &dVelocity_t1, double &dVelocity_t2)
{
	// cout << "calculating component length" << endl;

	Vector3d p, ey, aeq, t1, t2, r, et1, et2, er, hilf;
	Vector3d vt1, vt2, vr;

	p = vec3dPosition - vec3dMidpoint;

	// help vectors to calculate directions
	ey = Vector3d(0.0, 1.0, 0.0);
	hilf = Vector3d(0.000001, p.GetY(), 0);
	aeq = p - hilf;

//	cout << "ey, aeq:" << endl;
//	cout << ey << endl;
//	cout << aeq << endl;

	// calculate direction vectors
	t1 = ey.CrossProduct(aeq);
	t2 = t1.CrossProduct(p);
	r  = t1.CrossProduct(t2);

//	cout << "t1, t2, r:" << endl;
//	cout << t1 << endl;
//	cout << t2 << endl;
//	cout << r << endl;

	// normalize vectors
	er  = r.Normalize();
	et1 = t1.Normalize();
	et2 = t2.Normalize();

//	cout << "er, et1, et2:" << endl;
//	cout << er << endl;
//	cout << et1 << endl;
//	cout << et2 << endl;

	// calculate vector components
	vr  = er  * (er  * (*this) );
	vt1 = et1 * (et1 * (*this) );
	vt2 = et2 * (et2 * (*this) );

//	cout << "vr, vt1, vt2:" << endl;
//	cout << vr << endl;
//	cout << vt1 << endl;
//	cout << vt2 << endl;
//
//	cout << endl;
//	cout << "Length(vr, vt1, vt2):" << endl;

	// calculate and store length of vector components
	dVelocity_r  = vr.GetLength();
	dVelocity_t1 = vt1.GetLength();
	dVelocity_t2 = vt2.GetLength();
}

std::ostream& operator<<(std::ostream& os, Vector3d & vec)
{
    os << '(' << vec.GetX() << ", " << vec.GetY() << ", " << vec.GetZ() << ')' << std::endl;

    return os;
}
