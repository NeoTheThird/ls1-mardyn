/*
 * Region.cpp
 *
 *  Created on: 30.08.2013
 *      Author: mheinen
 */

#include "io/Region.h"

Region::Region(double* dLowerCorner, double* dUpperCorner, unsigned short nType)
{
	// init values
	_dLowerCorner[0] = dLowerCorner[0];
	_dLowerCorner[1] = dLowerCorner[1];
	_dLowerCorner[2] = dLowerCorner[2];

	_dUpperCorner[0] = dUpperCorner[0];
	_dUpperCorner[1] = dUpperCorner[1];
	_dUpperCorner[2] = dUpperCorner[2];

	_nType = nType;
}

Region::~Region()
{
}

bool Region::MoleculeIsInside(Molecule* mol)
{
	double dX, dY, dZ;

	dX = mol->r(0);
	dY = mol->r(1);
	dZ = mol->r(2);

	return dX > _dLowerCorner[0] && dX < _dUpperCorner[0] &&
           dY > _dLowerCorner[1] && dY < _dUpperCorner[1] &&
           dZ > _dLowerCorner[2] && dZ < _dUpperCorner[2];
}



