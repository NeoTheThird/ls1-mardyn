/*
 * Vector3d.h
 *
 *  Created on: 29.10.2013
 *      Author: mheinen
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

// sphere coordinates defined like in literature:
// title:  Matemathik fuer Ingenieure
// author: Lothar Papula
// pages:  123 ff.

// only axes x, y, z are switched:

// Papula      Vector3d
// ------      --------
//    x    -->    z
//    y    -->    x
//    z    -->    y

//  theta  -->   theta
//  phi    -->   phi


#include <iostream>

enum Vector3dOptions
{
	VEC3DOPT_CARTESIAN_COORDS = 1,
	VEC3DOPT_SPHERE_COORDS_RAD,
	VEC3DOPT_SPHERE_COORDS_DEGREE,
};

class Vector3d;

typedef double (Vector3d::*func_ptr_double)(void);
typedef int (Vector3d::*func_ptr_int)(void);

class Vector3d
{
	friend std::ostream& operator<<(std::ostream&, Vector3d &);

public:
	Vector3d();
	Vector3d(double dX, double dY, double dZ, int nOption = VEC3DOPT_CARTESIAN_COORDS);
	// Vector3d(double dRadius, double dTheta, double dPhi, unsigned int nOption);  // sphere coordinates
	Vector3d(Vector3d & rhs);
	Vector3d(Vector3d& vecOrigin, Vector3d& vecEnd, int nOption = VEC3DOPT_CARTESIAN_COORDS);
	~Vector3d();

	// cartesian coords
	double GetX();
	double GetY();
	double GetZ();

	// get sphere coordinates
	double GetRadius();
	double GetTheta(int nOption = VEC3DOPT_SPHERE_COORDS_RAD);
	double GetPhi(int nOption = VEC3DOPT_SPHERE_COORDS_RAD);

	Vector3d operator+( Vector3d & rhs );
	Vector3d operator-( Vector3d & rhs );
	double operator*(Vector3d & rhs);
	Vector3d operator*(double dFac);

	Vector3d CrossProduct(Vector3d & vec);
	Vector3d Normalize();
	double GetLength();
	void CalcRadialTangentialComponentLength( Vector3d &vec3dMidpoint, Vector3d &vec3dPosition,
			                                  double &dVelocity_r, double &dVelocity_t1, double &dVelocity_t2);

	void CalcRadialTangentialComponentLength( Vector3d &M, Vector3d &P, double &dv_r, double &dv_t1, double &dv_t2,
			                                  Vector3d &vr, Vector3d &vt1, Vector3d &vt2 );

	// translate coordinates
	void CalcCartesianCoords();
	void CalcSphereCoords();
/*
	// check coords calculation
	bool SphereCoordsCalculated() { return _bSphereCoordsCalculated; }
	bool CartesianCoordsCalculated() { return _bCartesianCoordsCalculated; }
*/
	unsigned int GetOption(){return _nOption;}
	bool IsEqualTo(Vector3d &vec, double dDiffMaxCartesianCoords = 0.0001);

	// quadrants
	int GetQuadrant();

	// absolute values
	double GetXAbsoluteValue();
	double GetYAbsoluteValue();
	double GetZAbsoluteValue();

private:

	// cartesian coordinates
	double ReturnX();
	double CalcX();
	double ReturnY();
	double CalcY();
	double ReturnZ();
	double CalcZ();

	// sphere coordinates
	double ReturnRadius();
	double CalcRadius();
	double ReturnTheta();
	double CalcTheta();
	double ReturnPhi();
	double CalcPhi();

	// quadrants
	int DetermineQuadrant();
	int ReturnQuadrant();

private:
	// cartesian coordinates
	double _dX;
	double _dY;
	double _dZ;

	// sphere coordinates
	double _dRadius;
	double _dTheta;
	double _dPhi;

	func_ptr_double _fptrX;
	func_ptr_double _fptrY;
	func_ptr_double _fptrZ;

	func_ptr_double _fptrRadius;
	func_ptr_double _fptrTheta;
	func_ptr_double _fptrPhi;

	int _nOption;

	// quadrants
	int _nQuadrant;
	func_ptr_int _fptrQuadrant;

};  // class Vector3d


#endif /* VECTOR3D_H_ */
