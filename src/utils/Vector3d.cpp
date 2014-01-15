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

const long double pi = 3.14159265358979323846264338327950288419716939937510;

Vector3d::Vector3d()
{
	_dX = 1.0;
	_dY = 1.0;
	_dZ = 1.0;

	_dRadius = 1.0;
	_dTheta  = 45.0;
	_dPhi    = 45.0;

	_fptrX = &Vector3d::ReturnX;
	_fptrY = &Vector3d::ReturnY;
	_fptrZ = &Vector3d::ReturnZ;

	_fptrRadius = &Vector3d::ReturnRadius;
	_fptrTheta  = &Vector3d::ReturnTheta;
	_fptrPhi    = &Vector3d::ReturnPhi;

	// quadrant
	_fptrQuadrant = &Vector3d::DetermineQuadrant;
}

Vector3d::Vector3d(double dX, double dY, double dZ, int nOption)
{
	// store nOption
	_nOption = nOption;

	// quadrant
	_fptrQuadrant = &Vector3d::DetermineQuadrant;

	// cartesian coordinates
	if( nOption == VEC3DOPT_CARTESIAN_COORDS )
	{
		_dX = dX;
		_dY = dY;
		_dZ = dZ;

		_fptrX = &Vector3d::ReturnX;
		_fptrY = &Vector3d::ReturnY;
		_fptrZ = &Vector3d::ReturnZ;

		_fptrRadius = &Vector3d::CalcRadius;
		_fptrTheta  = &Vector3d::CalcTheta;
		_fptrPhi    = &Vector3d::CalcPhi;
	}

	// sphere coordinates
	else if( nOption == VEC3DOPT_SPHERE_COORDS_RAD )
	{
		_dRadius = dX;  // radius
		_dTheta  = dY;  // theta
		_dPhi    = dZ;  // phi

		// angles > 180° restricted
		while(_dTheta > pi)
			_dTheta -= pi;

		// angles > 360° restricted
		while(_dPhi > 2*pi)
			_dPhi -= 2*pi;

		_fptrX = &Vector3d::CalcX;
		_fptrY = &Vector3d::CalcY;
		_fptrZ = &Vector3d::CalcZ;

		_fptrRadius = &Vector3d::ReturnRadius;
		_fptrTheta  = &Vector3d::ReturnTheta;
		_fptrPhi    = &Vector3d::ReturnPhi;
	}

	else if( nOption == VEC3DOPT_SPHERE_COORDS_DEGREE )
	{
		_dRadius = dX;  // radius
		_dTheta  = dY * pi / 180;  // theta
		_dPhi    = dZ * pi / 180;  // phi

		// angles > 180° restricted
		while(_dTheta > pi)
			_dTheta -= pi;

		// angles > 360° restricted
		while(_dPhi > 2*pi)
			_dPhi -= 2*pi;

		_fptrX = &Vector3d::CalcX;
		_fptrY = &Vector3d::CalcY;
		_fptrZ = &Vector3d::CalcZ;

		_fptrRadius = &Vector3d::ReturnRadius;
		_fptrTheta  = &Vector3d::ReturnTheta;
		_fptrPhi    = &Vector3d::ReturnPhi;
	}

}

Vector3d::Vector3d(Vector3d & rhs)
{
	// quadrant
	_nQuadrant = rhs.GetQuadrant();
	// alternativ (performance): _fptrQuadrant = &Vector3d::DetermineQuadrant;

	// cartesian coordinates
	if( rhs.GetOption() == VEC3DOPT_CARTESIAN_COORDS )
	{
		_dX = rhs.GetX();
		_dY = rhs.GetY();
		_dZ = rhs.GetZ();

		_fptrX = &Vector3d::ReturnX;
		_fptrY = &Vector3d::ReturnY;
		_fptrZ = &Vector3d::ReturnZ;

		_fptrRadius = &Vector3d::CalcRadius;
		_fptrTheta  = &Vector3d::CalcTheta;
		_fptrPhi    = &Vector3d::CalcPhi;
	}
	// sphere coordinates
	else
	{
		_dRadius = rhs.GetRadius();
		_dTheta  = rhs.GetTheta();
		_dPhi    = rhs.GetPhi();

		_fptrX = &Vector3d::CalcX;
		_fptrY = &Vector3d::CalcY;
		_fptrZ = &Vector3d::CalcZ;

		_fptrRadius = &Vector3d::ReturnRadius;
		_fptrTheta  = &Vector3d::ReturnTheta;
		_fptrPhi    = &Vector3d::ReturnPhi;
	}

	// nOption
	_nOption = rhs.GetOption();

}

Vector3d::Vector3d(Vector3d& vecOrigin, Vector3d& vecEnd, int nOption)
{
	// nOption
	_nOption = nOption;

	// quadrant
	_fptrQuadrant = &Vector3d::DetermineQuadrant;

	Vector3d vec;
	vec = vecEnd - vecOrigin;

	// cartesian coordinates
	_dX = vec.GetX();
	_dY = vec.GetY();
	_dZ = vec.GetZ();

	_fptrX = &Vector3d::ReturnX;
	_fptrY = &Vector3d::ReturnY;
	_fptrZ = &Vector3d::ReturnZ;

	_fptrRadius = &Vector3d::CalcRadius;
	_fptrTheta  = &Vector3d::CalcTheta;
	_fptrPhi    = &Vector3d::CalcPhi;

	// sphere coordinates
	if( nOption == VEC3DOPT_SPHERE_COORDS_RAD )
	{
		_dRadius = vec.GetRadius();
		_dTheta  = vec.GetTheta();
		_dPhi    = vec.GetPhi();

		_fptrRadius = &Vector3d::ReturnRadius;
		_fptrTheta  = &Vector3d::ReturnTheta;
		_fptrPhi    = &Vector3d::ReturnPhi;
	}

	else if( nOption == VEC3DOPT_SPHERE_COORDS_DEGREE )
	{
		_dRadius = vec.GetRadius();
		_dTheta  = vec.GetTheta();
		_dPhi    = vec.GetPhi();

		_fptrRadius = &Vector3d::ReturnRadius;
		_fptrTheta  = &Vector3d::ReturnTheta;
		_fptrPhi    = &Vector3d::ReturnPhi;
	}
}

Vector3d::~Vector3d()
{
}



Vector3d Vector3d::operator+ (Vector3d & rhs)
{
	Vector3d vec;
	double dX, dY, dZ;

	dX = this->GetX() + rhs.GetX();
	dY = this->GetY() + rhs.GetY();
	dZ = this->GetZ() + rhs.GetZ();

	vec = Vector3d(dX, dY, dZ);

	return vec;
}


Vector3d Vector3d::operator- (Vector3d & rhs)
{
	Vector3d vec;
	double dX, dY, dZ;

	dX = this->GetX() - rhs.GetX();
	dY = this->GetY() - rhs.GetY();
	dZ = this->GetZ() - rhs.GetZ();

	vec = Vector3d(dX, dY, dZ);

	return vec;
}

double Vector3d::operator* (Vector3d & rhs)
{
	double dScalar;

	dScalar = this->GetX() * rhs.GetX() + this->GetY() * rhs.GetY() + this->GetZ() * rhs.GetZ();

	return dScalar;
}


Vector3d Vector3d::operator*(double dFac)
{
	double dX, dY, dZ;

	dX = this->GetX() * dFac;
	dY = this->GetY() * dFac;
	dZ = this->GetZ() * dFac;

	Vector3d res(dX, dY, dZ);

	return res;
}


Vector3d Vector3d::CrossProduct(Vector3d & vec)
{
	double dX, dY, dZ;

	dX = this->GetY()*vec.GetZ() - this->GetZ()*vec.GetY();
	dY = this->GetZ()*vec.GetX() - this->GetX()*vec.GetZ();
	dZ = this->GetX()*vec.GetY() - this->GetY()*vec.GetX();

	Vector3d res(dX, dY, dZ);

	return res;
}

Vector3d Vector3d::Normalize()
{
	double dX, dY, dZ;

	dX = this->GetX() / this->GetLength();
	dY = this->GetY() / this->GetLength();
	dZ = this->GetZ() / this->GetLength();

	Vector3d res(dX, dY, dZ);

	return res;
}

double Vector3d::GetLength()
{
	double dLength;

	dLength = sqrt( this->GetX()*this->GetX() + this->GetY()*this->GetY() + this->GetZ()*this->GetZ() );

	return dLength;
}

void Vector3d::CalcRadialTangentialComponentLength( Vector3d &vec3dMidpoint, Vector3d &vec3dPosition,
		                                            double &dVelocity_r, double &dVelocity_t1, double &dVelocity_t2)
{
	// cout << "calculating component length" << endl;

	Vector3d p, ey, aeq, t1, t2, r, et1, et2, er, hilf;
	Vector3d vt1, vt2, vr;

	p = vec3dPosition - vec3dMidpoint;

	// help vectors to calculate directions
	ey = Vector3d(0.0, 1.0, 0.0);
	hilf = Vector3d(0.00000001, p.GetY(), 0);
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


void Vector3d::CalcRadialTangentialComponentLength( Vector3d &M, Vector3d &P, double &dv_r, double &dv_t1, double &dv_t2,
		                                            Vector3d &vr, Vector3d &vt1, Vector3d &vt2 )
{
	// cout << "calculating component length" << endl;

	Vector3d r, ey, aeq, evt1, evt2, evr, hilf;

	r = P - M;

	// help vectors to calculate directions
	ey = Vector3d(0.0, 1.0, 0.0);
	hilf = Vector3d(0.00000001, r.GetY(), 0);
	aeq = r - hilf;  // vector on aequator layer

	// calculate direction vectors, normalize
	vt1 = ey.CrossProduct(aeq);
	evt1 = vt1.Normalize();

	vt2 = vt1.CrossProduct(r);
	evt2 = vt2.Normalize();

	evr = r.Normalize();

	// calculate vector components
	vr  = evr  * (evr  * (*this) );
	vt1 = evt1 * (evt1 * (*this) );
	vt2 = evt2 * (evt2 * (*this) );

	// calculate and store length of vector components
	dv_r  = vr.GetLength();
	dv_t1 = vt1.GetLength();
	dv_t2 = vt2.GetLength();
}


// cartesian coords

// X
double Vector3d::ReturnX()
{
	return _dX;
}

double Vector3d::CalcX()
{
	_dX = this->GetRadius() * sin( this->GetTheta() ) * sin( this->GetPhi() );
	_fptrX = &Vector3d::ReturnX;

	return _dX;
}

double Vector3d::GetX()
{
	return (this->*_fptrX)();
}

// Y
double Vector3d::ReturnY()
{
	return _dY;
}

double Vector3d::CalcY()
{
	_dY = this->GetRadius() * cos( this->GetTheta() );
	_fptrY = &Vector3d::ReturnY;

	return _dY;
}

double Vector3d::GetY()
{
	return (this->*_fptrY)();
}

// Z
double Vector3d::ReturnZ()
{
	return _dZ;
}

double Vector3d::CalcZ()
{
	_dZ = this->GetRadius() * sin( this->GetTheta() ) * cos( this->GetPhi() );
	_fptrZ = &Vector3d::ReturnZ;

	return _dZ;
}

double Vector3d::GetZ()
{
	return (this->*_fptrZ)();
}

// get sphere coordinates

// Radius
double Vector3d::ReturnRadius()
{
	//std::cout << "called ReturnRadius()" << std::endl;
	return _dRadius;
}

double Vector3d::CalcRadius()
{
	//std::cout << "called CalcRadius()" << std::endl;

	_dRadius = sqrt( this->GetX()*this->GetX() + this->GetY()*this->GetY() + this->GetZ()*this->GetZ() );

	// set function pointer
	_fptrRadius = &Vector3d::ReturnRadius;

	return _dRadius;
}

double Vector3d::GetRadius()
{
	return (this->*_fptrRadius)();
}

// Theta
double Vector3d::ReturnTheta()
{
	return _dTheta;
}

double Vector3d::CalcTheta()
{
	// calc theta
	double dY, dR;

	dY = this->GetY();
	dR = this->GetRadius();
/*
	double dQuadrant[8] = {0, 0, 2*pi, 2*pi, 0, 0, 2*pi, 2*pi};
	double dHelpFact[8] = {1, 1, -1, -1, 1, 1, -1, -1};
	int nQuadrant = this->GetQuadrant();

	_dTheta = ( acos(dY / dR ) - dQuadrant[nQuadrant] ) * dHelpFact[nQuadrant];
*/
	_dTheta = acos(dY / (dR + 0.00000000000001) );

	// set function pointer
	_fptrTheta = &Vector3d::ReturnTheta;

	return _dTheta;
}

double Vector3d::GetTheta(int nOption)
{
	double dRet = (this->*_fptrTheta)();

	if( nOption == VEC3DOPT_SPHERE_COORDS_DEGREE )
		dRet = dRet / pi * 180.0;

	return dRet;
}

// Phi
double Vector3d::ReturnPhi()
{
	return _dPhi;
}

double Vector3d::CalcPhi()
{
	// calc Phi
	double dX, dZ;
	dX = this->GetX();
	dZ = this->GetZ();
	double nQuadrant[8] = {0, pi, pi, 2*pi, 0, pi, pi, 2*pi};

	_dPhi = atan( dX / (dZ + 0.00000000000001) ) + nQuadrant[ this->GetQuadrant() ];

	// set function pointer
	_fptrPhi = &Vector3d::ReturnPhi;

	return _dPhi;
}

double Vector3d::GetPhi(int nOption)
{
	double dRet = (this->*_fptrPhi)();

	if( nOption == VEC3DOPT_SPHERE_COORDS_DEGREE )
		dRet = dRet / pi * 180.0;

	return dRet;
}

/*
void Vector3d::CalcCartesianCoords()
{
	if(_bCartesianCoordsCalculated)
		return;

	_dX = _dRadius * sin(_dTheta) * sin(_dPhi);
	_dY = _dRadius * cos(_dTheta);
	_dZ = _dRadius * sin(_dTheta) * cos(_dPhi);

	_bCartesianCoordsCalculated = true;
}

void Vector3d::CalcSphereCoords()
{
	if(_bSphereCoordsCalculated)
		return;

	_dRadius = this->GetLength();
	_dTheta = acos(_dY / _dRadius);
	_dPhi = atan(_dX / _dZ);

	_bSphereCoordsCalculated = true;
}
*/

bool Vector3d::IsEqualTo(Vector3d &vec, double dDiffMaxCartesianCoords)
{
	bool bIsEqual = true;
	double dDeltaX, dDeltaY, dDeltaZ;

	// calc differences
	dDeltaX = this->GetX() - vec.GetX();
	dDeltaY = this->GetY() - vec.GetY();
	dDeltaZ = this->GetZ() - vec.GetZ();

	// calc absolute value
	dDeltaX -= 2.0*dDeltaX * (dDeltaX < 0.0);
	dDeltaY -= 2.0*dDeltaY * (dDeltaY < 0.0);
	dDeltaZ -= 2.0*dDeltaZ * (dDeltaZ < 0.0);

	// test acceptance
	bIsEqual = bIsEqual && dDeltaX < dDiffMaxCartesianCoords;
	bIsEqual = bIsEqual && dDeltaY < dDiffMaxCartesianCoords;
	bIsEqual = bIsEqual && dDeltaZ < dDiffMaxCartesianCoords;

	return bIsEqual;
}


// quadrant

int Vector3d::ReturnQuadrant()
{
	return _nQuadrant;
}

int Vector3d::DetermineQuadrant()
{
	double dX, dY, dZ;

	dX = this->GetX();
	dY = this->GetY();
	dZ = this->GetZ();

	_nQuadrant = 0;

	if(dY >= 0)
	{
		if(dX >= 0)
		{
			if(dZ >= 0)
				_nQuadrant = 0;
			else
				_nQuadrant = 1;
		}
		else
		{
			if(dZ >= 0)
				_nQuadrant = 3;
			else
				_nQuadrant = 2;
		}
	}
	else
	{
		if(dX >= 0)
		{
			if(dZ >= 0)
				_nQuadrant = 4;
			else
				_nQuadrant = 5;
		}
		else
		{
			if(dZ >= 0)
				_nQuadrant = 7;
			else
				_nQuadrant = 6;
		}
	}

	_fptrQuadrant = &Vector3d::ReturnQuadrant;

	return _nQuadrant;
}

int Vector3d::GetQuadrant()
{
	return (this->*_fptrQuadrant)();
}


// absolute values
double Vector3d::GetXAbsoluteValue()
{
	double dX;

	dX = this->GetX();
	dX -= 2.0*dX * (dX < 0.0);

	return dX;
}

double Vector3d::GetYAbsoluteValue()
{
	double dY;

	dY = this->GetY();
	dY -= 2.0*dY * (dY < 0.0);

	return dY;
}

double Vector3d::GetZAbsoluteValue()
{
	double dZ;

	dZ = this->GetZ();
	dZ -= 2.0*dZ * (dZ < 0.0);

	return dZ;
}



std::ostream& operator<<(std::ostream& os, Vector3d & vec)
{
    os << "cartesian:   " << '(' << vec.GetX() << ", " << vec.GetY() << ", " << vec.GetZ() << ')' << std::endl;
    os << "sphere_rad: " << '(' << vec.GetRadius() << ", " << vec.GetTheta() << ", " << vec.GetPhi() << ')' << std::endl;
    os << "sphere_deg: " << '(' << vec.GetRadius() << ", " << vec.GetTheta(VEC3DOPT_SPHERE_COORDS_DEGREE) << ", " << vec.GetPhi(VEC3DOPT_SPHERE_COORDS_DEGREE) << ')' << std::endl;

    return os;
}
