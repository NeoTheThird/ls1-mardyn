/*
 * Region.h
 *
 *  Created on: 30.08.2013
 *      Author: mheinen
 */

#include "molecules/Molecule.h"

#ifndef REGION_H_
#define REGION_H_

class Region
{
public:
	Region(double* dLowerCorner, double* dUpperCorner, unsigned short nType);
	~Region();

	unsigned short GetType(void) {return _nType;}
	double GetLeftBoundary(void) {return _dLowerCorner[0];}
	double GetMidpoint(unsigned int d) {return (_dUpperCorner[d] - _dLowerCorner[d]) / 2.0;}
	double GetRightBoundary(void) {return _dUpperCorner[0];}
	bool MoleculeIsInside(Molecule* mol);

private:
	double _dLowerCorner[3];
	double _dUpperCorner[3];
	unsigned short int _nType;
};

#endif /* REGION_H_ */
