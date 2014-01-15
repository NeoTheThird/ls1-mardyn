//============================================================================
// Name        : boost_Vector3d_test.cpp
// Author      : mheinen
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#define BOOST_TEST_MODULE Vector3d
#include <boost/test/included/unit_test.hpp>
#include "Vector3d.h"

#define VECTOR_SIZE_MAX 1000
const long double pi = 3.14159265358979323846264338327950288419716939937510;

BOOST_AUTO_TEST_CASE( constructor_test )
{
	double dx[VECTOR_SIZE_MAX];
	double dy[VECTOR_SIZE_MAX];
	double dz[VECTOR_SIZE_MAX];

	Vector3d vec[VECTOR_SIZE_MAX];

	// create random radius vectors
	srand(time(NULL)); //initialisiert den pseudo-random number generator

	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		dx[i] = ( rand()%1000 + 1) / 100.0;
		dy[i] = ( rand()%1000 + 1) / 100.0;
		dz[i] = ( rand()%1000 + 1) / 100.0;

		BOOST_REQUIRE_MESSAGE( dx[i] > 0.0, "dx[" << i << "] = " << dx[i] );
		BOOST_REQUIRE_MESSAGE( dy[i] > 0.0, "dy[" << i << "] = " << dy[i] );
		BOOST_REQUIRE_MESSAGE( dz[i] > 0.0, "dz[" << i << "] = " << dz[i] );

		vec[i] = Vector3d(dx[i], dy[i], dz[i]);
	}

	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		// value stored correctly
		BOOST_REQUIRE( vec[i].GetX() == dx[i] );
		BOOST_REQUIRE( vec[i].GetY() == dy[i] );
		BOOST_REQUIRE( vec[i].GetZ() == dz[i] );


		// copy constructor
		Vector3d vec2(vec[i]);

		BOOST_REQUIRE( vec2.GetX() == vec[i].GetX() );
		BOOST_REQUIRE( vec2.GetY() == vec[i].GetY() );
		BOOST_REQUIRE( vec2.GetZ() == vec[i].GetZ() );
	}


	// test constructor options

	Vector3d vec1;
	double pi = 3.141592654;
	double dRadius, dTheta, dPhi;
	double dX, dY, dZ;

	vec1 = Vector3d(1, sqrt(2), 1, VEC3DOPT_CARTESIAN_COORDS);
	dRadius = vec1.GetRadius();
	dTheta  = vec1.GetTheta(VEC3DOPT_SPHERE_COORDS_DEGREE);
	dPhi    = vec1.GetPhi(VEC3DOPT_SPHERE_COORDS_DEGREE);

	BOOST_CHECK_CLOSE( dTheta, 45, 0.0001);
	BOOST_CHECK_CLOSE( dPhi,   45, 0.0001);

    if( !( abs(dRadius - 2) < 0.00001 && abs(dTheta - 45) < 0.00001 && abs(dPhi - 45) < 0.00001 ) )
    {
    	std::cout << vec1;
    	BOOST_FAIL( "Fehler!" );
    }

	vec1 = Vector3d(2, 45, 45, VEC3DOPT_SPHERE_COORDS_DEGREE);
	dRadius = vec1.GetRadius();
	dTheta = vec1.GetTheta(VEC3DOPT_SPHERE_COORDS_DEGREE);
	dPhi   = vec1.GetPhi(VEC3DOPT_SPHERE_COORDS_DEGREE);

	// angle theta and phi should equal 45°
	BOOST_CHECK_CLOSE( dTheta, 45, 0.0001);
	BOOST_CHECK_CLOSE( dPhi,   45, 0.0001);

	if( !( abs(dRadius - 2) < 0.00001 && abs(dTheta - 45) < 0.00001 && abs(dPhi - 45) < 0.00001 ) )
    {
    	std::cout << vec1;
    	BOOST_FAIL( "Fehler!" );
    }
}


BOOST_AUTO_TEST_CASE( get_absolute_value_test )
{
	Vector3d vec;
	double dX, dY, dZ;
	double dXabs, dYabs, dZabs;

	dX = -3.3;
	dY =  5.4;
	dZ = -2.8;

	dXabs = 3.3;
	dYabs = 5.4;
	dZabs = 2.8;

	vec = Vector3d(dX, dY, dZ);

	BOOST_REQUIRE_MESSAGE( vec.GetXAbsoluteValue() == dXabs,
			"vec.GetXAbsoluteValue() = " << vec.GetXAbsoluteValue() << ", dXabs = " << dXabs << ", dX = " << dX);
	BOOST_REQUIRE_MESSAGE( vec.GetYAbsoluteValue() == dYabs,
			"vec.GetYAbsoluteValue() = " << vec.GetYAbsoluteValue() << ", dYabs = " << dYabs << ", dY = " << dY);
	BOOST_REQUIRE_MESSAGE( vec.GetZAbsoluteValue() == dZabs,
			"vec.GetZAbsoluteValue() = " << vec.GetZAbsoluteValue() << ", dZabs = " << dZabs << ", dZ = " << dZ);
}

BOOST_AUTO_TEST_CASE( get_quadrant_test )
{
	Vector3d vec[8];
	int nQuadrant;

	vec[0] = Vector3d( 1,  1,  1);
	vec[1] = Vector3d( 1,  1, -1);
	vec[2] = Vector3d(-1,  1, -1);
	vec[3] = Vector3d(-1,  1,  1);

	vec[4] = Vector3d( 1, -1,  1);
	vec[5] = Vector3d( 1, -1, -1);
	vec[6] = Vector3d(-1, -1, -1);
	vec[7] = Vector3d(-1, -1,  1);

	for(int i=0; i<8; i++)
	{
		nQuadrant = vec[i].GetQuadrant();
		BOOST_REQUIRE_MESSAGE( i == nQuadrant, "i = " << i << ", nQuadrant = " << nQuadrant );
	}
}

BOOST_AUTO_TEST_CASE( is_equal_test )
{
	Vector3d vec1, vec2, vec3;
	double dx, dy, dz;
	double dRadius, dTheta, dPhi;

	dx = 0;
	dy = 7;
	dz = 0;

	dRadius = 7;
	dTheta  = 0;
	dPhi    = 0;

	vec1 = Vector3d(dx, dy, dz);
	vec2 = Vector3d(vec1);
	vec3 = Vector3d(dRadius, dTheta, dPhi, VEC3DOPT_SPHERE_COORDS_RAD);

	if( !vec1.IsEqualTo(vec3) )
	{
		std::cout << vec1 << std::endl;
		// std::cout << vec3 << std::endl;
		BOOST_FAIL( "vectors not equal!" );
	}


	// main directions, constructor: cartesian/sphere coordinates

	Vector3d vecCartes[26];
	Vector3d vecSphere[26];
	dRadius = 1.0;

	// create cartesian coordinates vectors
	vecCartes[0] = Vector3d( 1,  0,  0 );  //  x
	vecCartes[1] = Vector3d(-1,  0,  0 );  // -x
	vecCartes[2] = Vector3d( 0,  1,  0 );  //  y
	vecCartes[3] = Vector3d( 0, -1,  0 );  // -y
	vecCartes[4] = Vector3d( 0,  0,  1 );  //  z
	vecCartes[5] = Vector3d( 0,  0, -1 );  // -z

	// 45° x,y layer
	vecCartes[6] = Vector3d( 1,  1,  0 );  //  x,  y
	vecCartes[7] = Vector3d(-1,  1,  0 );  // -x,  y
	vecCartes[8] = Vector3d(-1, -1,  0 );  // -x, -y
	vecCartes[9] = Vector3d( 1, -1,  0 );  //  x, -Y

	// 45° x,z layer
	vecCartes[10] = Vector3d( 1,  0,  1 );  //  x,  z
	vecCartes[11] = Vector3d(-1,  0,  1 );  // -x,  z
	vecCartes[12] = Vector3d(-1,  0, -1 );  // -x, -z
	vecCartes[13] = Vector3d( 1,  0, -1 );  //  x, -z

	// 45° y,z layer
	vecCartes[14] = Vector3d( 0,  1,  1 );  //  y,  z
	vecCartes[15] = Vector3d( 0,  1, -1 );  //  y, -z
	vecCartes[16] = Vector3d( 0, -1, -1 );  // -y, -z
	vecCartes[17] = Vector3d( 0, -1,  1 );  // -y,  z

	{
		double y = sqrt(2);

		// 45° x, y, z
		vecCartes[18] = Vector3d( 1,  y,  1 );  //  x,  y,  z
		vecCartes[19] = Vector3d( 1,  y, -1 );  //  x,  y, -z
		vecCartes[20] = Vector3d(-1,  y, -1 );  // -x,  y, -z
		vecCartes[21] = Vector3d(-1,  y,  1 );  // -x,  y,  z

		vecCartes[22] = Vector3d( 1, -y,  1 );  //  x, -y,  z
		vecCartes[23] = Vector3d( 1, -y, -1 );  //  x, -y, -z
		vecCartes[24] = Vector3d(-1, -y, -1 );  // -x, -y, -z
		vecCartes[25] = Vector3d(-1, -y,  1 );  // -x, -y,  z
	}

	// create sphere coordinates vectors
	vecSphere[0] = Vector3d( dRadius, 90,  90,  VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x
	vecSphere[1] = Vector3d( dRadius, 90,  270, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x
	vecSphere[2] = Vector3d( dRadius, 0,   0,   VEC3DOPT_SPHERE_COORDS_DEGREE );  //  y
	vecSphere[3] = Vector3d( dRadius, 180, 0,   VEC3DOPT_SPHERE_COORDS_DEGREE );  // -y
	vecSphere[4] = Vector3d( dRadius, 90,  0,   VEC3DOPT_SPHERE_COORDS_DEGREE );  //  z
	vecSphere[5] = Vector3d( dRadius, 90,  180, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -z


	dRadius = sqrt(2);

	// 45° x,y layer
	vecSphere[6] = Vector3d( dRadius, 45,  90,  VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x,  y
	vecSphere[7] = Vector3d( dRadius, 45,  270, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x,  y
	vecSphere[8] = Vector3d( dRadius, 135, 270, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x, -y
	vecSphere[9] = Vector3d( dRadius, 135, 90,  VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x, -Y

	// 45° x,z layer
	vecSphere[10] = Vector3d( dRadius, 90, 45,  VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x, z
	vecSphere[11] = Vector3d( dRadius, 90, 315, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x, z
	vecSphere[12] = Vector3d( dRadius, 90, 225, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x,-z
	vecSphere[13] = Vector3d( dRadius, 90, 135, VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x,-z

	// 45° y,z layer
	vecSphere[14] = Vector3d( dRadius, 45,  0,   VEC3DOPT_SPHERE_COORDS_DEGREE );  //  y,  z
	vecSphere[15] = Vector3d( dRadius, 45,  180, VEC3DOPT_SPHERE_COORDS_DEGREE );  //  y, -z
	vecSphere[16] = Vector3d( dRadius, 135, 180, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -y, -z
	vecSphere[17] = Vector3d( dRadius, 135, 0,   VEC3DOPT_SPHERE_COORDS_DEGREE );  // -y,  z


	dRadius = 2;

	// 45° x, y, z
	vecSphere[18] = Vector3d( dRadius, 45, 45,  VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x,  y,  z
	vecSphere[19] = Vector3d( dRadius, 45, 135, VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x,  y, -z
	vecSphere[20] = Vector3d( dRadius, 45, 225, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x,  y, -z
	vecSphere[21] = Vector3d( dRadius, 45, 315, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x,  y,  z

	vecSphere[22] = Vector3d( dRadius, 135, 45,  VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x, -y,  z
	vecSphere[23] = Vector3d( dRadius, 135, 135, VEC3DOPT_SPHERE_COORDS_DEGREE );  //  x, -y, -z
	vecSphere[24] = Vector3d( dRadius, 135, 225, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x, -y, -z
	vecSphere[25] = Vector3d( dRadius, 135, 315, VEC3DOPT_SPHERE_COORDS_DEGREE );  // -x, -y,  z

	// in different ways constructed vectors equal?
	for(int i=0; i<26; i++)
	{
		if( !vecCartes[i].IsEqualTo( vecSphere[i] ) )
		{
			std::cout << "vecCartes[" << i << "]:" << std::endl << vecCartes[i] << std::endl;
			std::cout << "vecSphere[" << i << "]:" << std::endl << vecSphere[i] << std::endl;
			BOOST_FAIL( "vectors not equal!" );
		}
	}


	// different cases

	// vectors not equal, coords(vec1) > coords(vec2)

	vec1 = Vector3d(3, 3, 3);
	vec2 = Vector3d(2, 3, 3);

	if( vec1.IsEqualTo( vec2 ) )
	{
		std::cout << "vec1" << std::endl << vec1 << std::endl;
		std::cout << "vec2" << std::endl << vec2 << std::endl;
		BOOST_FAIL( "vectors equal, but shouldnt!" );
	}

	vec1 = Vector3d(3, 3, 3);
	vec2 = Vector3d(3, 2, 3);

	if( vec1.IsEqualTo( vec2 ) )
	{
		std::cout << "vec1" << std::endl << vec1 << std::endl;
		std::cout << "vec2" << std::endl << vec2 << std::endl;
		BOOST_FAIL( "vectors equal, but shouldnt!" );
	}

	vec1 = Vector3d(3, 3, 3);
	vec2 = Vector3d(3, 3, 2);

	if( vec1.IsEqualTo( vec2 ) )
	{
		std::cout << "vec1" << std::endl << vec1 << std::endl;
		std::cout << "vec2" << std::endl << vec2 << std::endl;
		BOOST_FAIL( "vectors equal, but shouldnt!" );
	}


	// vectors not equal, coords(vec1) < coords(vec2)

	vec1 = Vector3d(2, 3, 3);
	vec2 = Vector3d(3, 3, 3);

	if( vec1.IsEqualTo( vec2 ) )
	{
		std::cout << "vec1" << std::endl << vec1 << std::endl;
		std::cout << "vec2" << std::endl << vec2 << std::endl;
		BOOST_FAIL( "vectors equal, but shouldnt!" );
	}

	vec1 = Vector3d(3, 2, 3);
	vec2 = Vector3d(3, 3, 3);

	if( vec1.IsEqualTo( vec2 ) )
	{
		std::cout << "vec1" << std::endl << vec1 << std::endl;
		std::cout << "vec2" << std::endl << vec2 << std::endl;
		BOOST_FAIL( "vectors equal, but shouldnt!" );
	}

	vec1 = Vector3d(3, 3, 2);
	vec2 = Vector3d(3, 3, 3);

	if( vec1.IsEqualTo( vec2 ) )
	{
		std::cout << "vec1" << std::endl << vec1 << std::endl;
		std::cout << "vec2" << std::endl << vec2 << std::endl;
		BOOST_FAIL( "vectors equal, but shouldnt!" );
	}

}


BOOST_AUTO_TEST_CASE( math_operations_test )
{
	// add
	double x1, y1, z1;
	double dRadius, dTheta, dPhi;

	x1 = 3.0;
	y1 = 4.0;
	z1 = 5.0;

	dRadius = 7;
	dTheta  = 0;
	dPhi    = 0;

	Vector3d vec1(x1, y1, z1);
	Vector3d vec2(vec1);
	Vector3d vecSphere(dRadius, dTheta, dPhi, VEC3DOPT_SPHERE_COORDS_DEGREE);
	Vector3d vecRes;

	vecRes = vec1 + vec2;

	BOOST_CHECK_MESSAGE( vecRes.GetX() == 2 * vec1.GetX(), "x = " << vecRes.GetX() );
	BOOST_CHECK_MESSAGE( vecRes.GetY() == 2 * vec1.GetY(), "y = " << vecRes.GetY() );
	BOOST_CHECK_MESSAGE( vecRes.GetZ() == 2 * vec1.GetZ(), "z = " << vecRes.GetZ() );


	// sphere coords
	vecRes = vec1 + vecSphere;

	BOOST_CHECK_MESSAGE( vecRes.GetX() == vec1.GetX() + 0, "x = " << vecRes.GetX() << ", radius = " << vecSphere.GetRadius() );
	BOOST_CHECK_MESSAGE( vecRes.GetY() == vec1.GetY() + dRadius, "y = " << vecRes.GetY() << ", theta  = " << vecSphere.GetTheta() );
	BOOST_CHECK_MESSAGE( vecRes.GetZ() == vec1.GetZ() + 0, "z = " << vecRes.GetZ() << ", phi    = " << vecSphere.GetPhi() );


	// sub
	vecRes = vec1 - vec2;

	BOOST_CHECK_MESSAGE( vecRes.GetX() == 0, "x = " << vecRes.GetX() );
	BOOST_CHECK_MESSAGE( vecRes.GetY() == 0, "y = " << vecRes.GetY() );
	BOOST_CHECK_MESSAGE( vecRes.GetZ() == 0, "z = " << vecRes.GetZ() );


	// multiply by factor
	vecRes = vec1 * 2.0;

	BOOST_CHECK_MESSAGE( vecRes.GetX() == 2 * vec1.GetX(), "x = " << vecRes.GetX() );
	BOOST_CHECK_MESSAGE( vecRes.GetY() == 2 * vec1.GetY(), "y = " << vecRes.GetY() );
	BOOST_CHECK_MESSAGE( vecRes.GetZ() == 2 * vec1.GetZ(), "z = " << vecRes.GetZ() );


	// scalar product
	double dScalar = vec1 * vec2;
	double dRes = x1*x1 + y1*y1 + z1*z1;

	BOOST_CHECK_MESSAGE( dScalar == dRes, "IST = " << dScalar << " SOLL = " << dRes );


	// cross product
	Vector3d vecCross;
	vecCross = vec1.CrossProduct(vec2);
	dRes = vec1 * vecCross;

	// test if scalar product of result vector with origin vector == null vector
	BOOST_CHECK_MESSAGE( dRes == 0, "IST = " << dRes << ", Soll = 0");

}


BOOST_AUTO_TEST_CASE( normalize_test )
{
	// create random vectors
	srand(time(NULL)); //initialisiert den pseudo-random number generator

	Vector3d vec;
	Vector3d vecRes;
	double dx, dy, dz;

	for(int i=0; i<100; i++)
	{
		dx = (double)(rand()%1000) / 100.0;
		dy = (double)(rand()%1000) / 100.0;
		dz = (double)(rand()%1000) / 100.0;
		vec = Vector3d(dx, dy, dz);

		vecRes = vec.Normalize();

		BOOST_CHECK_CLOSE( vecRes.GetLength(), 1.0, 0.000001 );
	}
}


BOOST_AUTO_TEST_CASE( vector_breakdown_test )
{
	// velocity vectors in radial direction only

	Vector3d m(11, 12, 13);  // midpoint of sphere
	Vector3d p[VECTOR_SIZE_MAX];
	Vector3d r[VECTOR_SIZE_MAX];
	Vector3d v[VECTOR_SIZE_MAX];
	Vector3d v_t[VECTOR_SIZE_MAX];  // tangential direction only

	// velocity vector components for tangential velocity vectors only
	Vector3d vr[VECTOR_SIZE_MAX];
	Vector3d vt1[VECTOR_SIZE_MAX];
	Vector3d vt2[VECTOR_SIZE_MAX];

	// velocity vector components for tangential velocity vectors only
	Vector3d vr_t[VECTOR_SIZE_MAX];
	Vector3d vt1_t[VECTOR_SIZE_MAX];
	Vector3d vt2_t[VECTOR_SIZE_MAX];

	double dRadius = 7.0;

	// create specific radius vectors

	r[0] = Vector3d( 1,  0,  0 );  //  x
	r[1] = Vector3d(-1,  0,  0 );  // -x
	r[2] = Vector3d( 0,  1,  0 );  //  y
	r[3] = Vector3d( 0, -1,  0 );  // -y
	r[4] = Vector3d( 0,  0,  1 );  //  z
	r[5] = Vector3d( 0,  0, -1 );  // -z

	// 45° x,y layer
	r[6] = Vector3d( 1,  1,  0 );
	r[7] = Vector3d(-1,  1,  0 );
	r[8] = Vector3d(-1, -1,  0 );
	r[9] = Vector3d( 1, -1,  0 );

	// 45° x,z layer
	r[10] = Vector3d( 1,  0,  1 );
	r[11] = Vector3d(-1,  0,  1 );
	r[12] = Vector3d(-1,  0, -1 );
	r[13] = Vector3d( 1,  0, -1 );

	// 45° y,z layer
	r[14] = Vector3d( 0,  1,  1 );
	r[15] = Vector3d( 0,  1, -1 );
	r[16] = Vector3d( 0, -1, -1 );
	r[17] = Vector3d( 0, -1,  1 );


	// create random radius vectors
	srand(time(NULL)); //initialisiert den pseudo-random number generator

	double dx[VECTOR_SIZE_MAX];
	double dy[VECTOR_SIZE_MAX];
	double dz[VECTOR_SIZE_MAX];

	for(int i=18; i < VECTOR_SIZE_MAX; i++)
	{
		dx[i] = ( rand()%1000 + 10) / 100.0;
		dy[i] = ( rand()%1000 + 10) / 100.0;
		dz[i] = ( rand()%1000 + 10) / 100.0;
	}

	// check random values
	for(int i=18; i < VECTOR_SIZE_MAX; i++)
	{
		BOOST_REQUIRE_MESSAGE( dx[i] >= 0.1, "dx[" << i << "] = " << dx[i] );
		BOOST_REQUIRE_MESSAGE( dy[i] >= 0.1, "dy[" << i << "] = " << dy[i] );
		BOOST_REQUIRE_MESSAGE( dz[i] >= 0.1, "dz[" << i << "] = " << dz[i] );

		r[i] = Vector3d(dx[i], dy[i], dz[i]);
	}

	// check random vectors
	for(int i=18; i < VECTOR_SIZE_MAX; i++)
	{
		BOOST_CHECK_MESSAGE( r[i].GetX() >= 0.1, "_x = " << r[i].GetX() << ", dx = " << dx[i] );
		BOOST_CHECK_MESSAGE( r[i].GetY() >= 0.1, "_y = " << r[i].GetY() << ", dy = " << dy[i] );
		BOOST_CHECK_MESSAGE( r[i].GetZ() >= 0.1, "_z = " << r[i].GetZ() << ", dz = " << dz[i] );
	}

	// radius vector not null vector
	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		BOOST_REQUIRE_MESSAGE( r[i].GetLength() > 0.0, "x = " << r[i].GetX() <<
				                                     ", y = " << r[i].GetY() <<
				                                     ", z = " << r[i].GetZ() <<
				                                     ", length[" << i << "] = " << r[i].GetLength() );
	}

	// normalize, scale to radius
	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		r[i] = r[i].Normalize();
		r[i] = r[i] * dRadius;

		// length of radius vectors have to equal radius value
		BOOST_CHECK_CLOSE( r[i].GetLength(), dRadius, 0.000001 );
	}

	// velocity vectors in radial direction
	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		v[i] = r[i];

		BOOST_REQUIRE( v[i].GetX() == r[i].GetX() );
		BOOST_REQUIRE( v[i].GetY() == r[i].GetY() );
		BOOST_REQUIRE( v[i].GetZ() == r[i].GetZ() );
	}

	// velocity vectors in tangential direction
	Vector3d vecHelp;
	double dDeltaTheta_deg, dDeltaTheta_rad, dHelpRadius, dScalar;

	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		dDeltaTheta_deg = (rand()%5001 + 1000) / 100.0;
		BOOST_REQUIRE(dDeltaTheta_deg >= 10 && dDeltaTheta_deg <= 60);  // angle between: 10.00° .. 60.00°

		dDeltaTheta_rad = dDeltaTheta_deg * pi/180;
		BOOST_REQUIRE_MESSAGE(dDeltaTheta_rad > 0.17 && dDeltaTheta_rad < 1.05,
				"dDeltaTheta_deg = " << dDeltaTheta_deg << ", dDeltaTheta_rad = " << dDeltaTheta_rad );  // radian between: 0.17 .. 1.05

		dHelpRadius = r[i].GetRadius() / cos(dDeltaTheta_rad);
		BOOST_REQUIRE(dHelpRadius > 0);

		// theta > 180° restricted
		if( (r[i].GetTheta(VEC3DOPT_SPHERE_COORDS_DEGREE) + dDeltaTheta_deg) < 180 )
			vecHelp = Vector3d( dHelpRadius, r[i].GetTheta() + dDeltaTheta_rad, r[i].GetPhi(), VEC3DOPT_SPHERE_COORDS_RAD );
		else
			vecHelp = Vector3d( dHelpRadius, r[i].GetTheta() - dDeltaTheta_rad, r[i].GetPhi(), VEC3DOPT_SPHERE_COORDS_RAD );

		v_t[i] = Vector3d(r[i], vecHelp);

		// scalar product should be 0
		dScalar = v_t[i] * r[i];

		if(dScalar > 0.00001 || dScalar < -0.00001)
		{
			std::cout << "dScalar = " << dScalar << " SOLL = 0" << std::endl;
			std::cout << "dDeltaTheta_deg = " << dDeltaTheta_deg << std::endl;
			std::cout << "dDeltaTheta_rad = " << dDeltaTheta_rad << std::endl;
			std::cout << "dHelpRadius = " << dHelpRadius << std::endl;

			std::cout << "r[i]:" << std::endl << r[i];
			std::cout << "vecHelp:" << std::endl << vecHelp;
			std::cout << "v_t[i]:" << std::endl << v_t[i];

			BOOST_FAIL("Fehler!");
		}


		// vectors valid and not null vector
		double abs_x, abs_y, abs_z;

		abs_x = v_t[i].GetXAbsoluteValue();
		abs_y = v_t[i].GetYAbsoluteValue();
		abs_z = v_t[i].GetZAbsoluteValue();

		bool bIsNullVector = (abs_x < 0.000001) && (abs_y < 0.000001) && (abs_z < 0.000001);

		if( bIsNullVector == true )
		{
			std::cout << "x = " << v_t[i].GetX() << ", y = " << v_t[i].GetY() << ", z = " << v_t[i].GetZ() << std::endl;
			std::cout << "abs_x = " << abs_x << ", abs_y = " << abs_y << ", abs_z = " << abs_z << std::endl;

			std::cout << "r[" << i << "]:" << std::endl << r[i] << std::endl;
			std::cout << "vecHelp:" << std::endl << vecHelp;
			std::cout << "v_t[" << i << "]:" << std::endl << v_t[i] << std::endl;

			BOOST_FAIL("Fehler!");
		}
	}

	// calculate position vectors
	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		p[i] = m + r[i];

		// distance from origin of each point at least equals dist(midpoint) - radius
		BOOST_CHECK( p[i].GetLength() >= m.GetLength() - dRadius );
	}

	// breakdown velocity vectors
	double dvr, dvt1, dvt2;

	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		// radial velocity vectors
		v[i].CalcRadialTangentialComponentLength(p[i], m, dvr, dvt1, dvt2, vr[i], vt1[i], vt2[i]);

		// length of velocity components t1, t2 should equal 0, that of r should equal length of v
		BOOST_CHECK( dvt1 < 0.000001 );
		BOOST_CHECK( dvt2 < 0.000001 );
		BOOST_CHECK_CLOSE( dvr, v[i].GetLength(), 0.000001);

		BOOST_CHECK_EQUAL( dvr,  vr[i].GetLength() );
		BOOST_CHECK_EQUAL( dvt1, vt1[i].GetLength() );
		BOOST_CHECK_EQUAL( dvt2, vt2[i].GetLength() );


		// tangential velocity vectors
		v_t[i].CalcRadialTangentialComponentLength(p[i], m, dvr, dvt1, dvt2, vr_t[i], vt1_t[i], vt2_t[i]);

		// length of velocity components t1, t2 should be > 0, that of r should equal 0
		BOOST_CHECK( (dvt1 + dvt2) > 0.000001 );
		BOOST_CHECK( dvr < 0.000001 );

		if( abs(dvt1 + dvt2) < 0.000001 || abs(dvr) > 0.000001 )
		{
			std::cout << "dvt1 = " << dvt1 << ", dvt2 = " << dvt2 << ", dvr = " << dvr << std::endl;

			std::cout << "r[" << i << "]:" << std::endl;
			std::cout << r[i] << std::endl;
			std::cout << "v_t[" << i << "]:" << std::endl;
			std::cout << v_t[i] << std::endl;
			std::cout << "vr_t[" << i << "]:" << std::endl;
			std::cout << vr_t[i] << std::endl;
			std::cout << "vt1_t[" << i << "]:" << std::endl;
			std::cout << vt1_t[i] << std::endl;
			std::cout << "vt2_t[" << i << "]:" << std::endl;
			std::cout << vt2_t[i] << std::endl;
		}

		BOOST_CHECK_EQUAL( dvr,  vr_t[i].GetLength() );
		BOOST_CHECK_EQUAL( dvt1, vt1_t[i].GetLength() );
		BOOST_CHECK_EQUAL( dvt2, vt2_t[i].GetLength() );
	}


	// scalar products of components should equal 0
	dScalar = 0;

	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		// radial velocity vectors
		dScalar = vr[i] * vt1[i];
		BOOST_CHECK_MESSAGE( dScalar < 0.000001, "dScalar = " << dScalar );

		dScalar = vr[i] * vt2[i];
		BOOST_CHECK_MESSAGE( dScalar < 0.000001, "dScalar = " << dScalar );

		dScalar = vt1[i] * vt2[i];
		BOOST_CHECK_MESSAGE( dScalar < 0.000001, "dScalar = " << dScalar );


		// tangential velocity vectors only
		dScalar = vr_t[i] * vt1_t[i];
		BOOST_CHECK_MESSAGE( dScalar < 0.000001, "dScalar = " << dScalar );

		dScalar = vr_t[i] * vt2_t[i];
		BOOST_CHECK_MESSAGE( dScalar < 0.000001, "dScalar = " << dScalar );

		dScalar = vt1_t[i] * vt2_t[i];
		BOOST_CHECK_MESSAGE( dScalar < 0.000001, "dScalar = " << dScalar );
	}


	// sum of vector components should reproduce origin vector

	Vector3d vSum;

	for(int i=0; i < VECTOR_SIZE_MAX; i++)
	{
		vSum = vr[i] + vt1[i] + vt2[i];

		BOOST_CHECK( v[i].GetX() - vSum.GetX() < 0.0000001 );
		BOOST_CHECK( v[i].GetY() - vSum.GetY() < 0.0000001 );
		BOOST_CHECK( v[i].GetZ() - vSum.GetZ() < 0.0000001 );
	}

}



/*
    // seven ways to detect and report the same error:
    BOOST_CHECK( add( a,b ) == 4 );        // #1 continues on error

    BOOST_REQUIRE( add( a,b ) == 4 );      // #2 throws on error

    if( add( 2,2 ) != 4 )
      BOOST_ERROR( "Ouch..." );            // #3 continues on error

    if( add( 2,2 ) != 4 )
      BOOST_FAIL( "Ouch..." );             // #4 throws on error

    if( add( 2,2 ) != 4 ) throw "Ouch..."; // #5 throws on error

    BOOST_CHECK_MESSAGE( add( 2,2 ) == 4,  // #6 continues on error
                         "add(..) result: " << add( 2,2 ) );

    BOOST_CHECK_EQUAL( add( a,b ), 4 );	  // #7 continues on error
*/

