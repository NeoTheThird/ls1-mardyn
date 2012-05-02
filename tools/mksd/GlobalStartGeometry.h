/*
 * GlobalStartGeometry.h
 *
 *  Created on: 12.01.2012
 *      Author:Stefan Becker <stefan.becker@mv.uni-kl.de>
 *
 * The instance of this class defines the global geometry of the start configuration, e.g. the dimensions of the simulation box
 * or the length of the liquid cuboid in relation to the box length, amongst others.
 */

#ifndef GLOBALSTARTGEOMETRY_H_
#define GLOBALSTARTGEOMETRY_H_

#include<iostream>
#include<string>
#include<math.h>
#include<map>
#include"RandomNumber.h"

using namespace std;

class GlobalStartGeometry
{
public:
	GlobalStartGeometry(unsigned in_nFluid, double in_rhoLiq, double in_rhoVap, double in_alpha, double in_beta, double in_gamma);
	~GlobalStartGeometry();

	void calculateBoxFluidOffset(double hWall, double shielding);
	void calculateFillProbabilityArray();

	double gBoxLength(unsigned direction);
	unsigned gFluidUnits(unsigned direction);
	double gOffset(unsigned direction);
	double gFluidUnit(unsigned direction);
	bool gFillArray(unsigned fluidUnits0, unsigned fluidUnits1, unsigned flunidUnits2, unsigned particleInElementaryBox);
	unsigned gNFilledSlots();

private:
	unsigned _nFluid;
	//unsigned _wallThick; // thickness of the wall as a multiple of the lattice constant in y-direction
	double _alpha;
	double _beta;
	double _gamma;

	//@brief: the next three variables needed to compute the density of the fluid (so far: 1CLJ)
	double _rhoLiq;
	double _rhoVap;
	double _sigmaFluid1CLJ;

	double _nLiq;
	double _nVap;
	unsigned _nFilledSlots;

	//@brief: lattice contant of the solid phase (wall) in the three dimensions
	//double _latticeConstSolid[3];
	//@brief: dimensions of the simulation box
	double _box[3];
	//@brief: offset of the fluid particles posistion with respect to the simulation box (x,z direction)
	// and with respect to the wall (y-direction)
	double _off[3];
	//@brief: number of elementary liquid lattices in three directions
	unsigned _fluidUnits[3];
	//@brief: lengths of a single elementary liquid lattice in three directions
	double _fluidUnit[3];
	//@brief: gross probability of a liquid elementary box to be filled
	double _fluidFillProbability;
	//@brief: 4-dimensional array addressing single slots that may be filled with a particle
	// => 3 diemsions for a fluid elementary box (due to 3 directions in space)
	// and one dimension addressing one of three slots within an elementary box
	map< unsigned, map<unsigned, map<unsigned, map<unsigned, bool> > > > _fill;
	//bool _fill[_fluidUnits[0]][_fluidUnits[1]][_fluidUnits[2]][3];
	// vector< vector< vector < vector <bool> > > > _fill;
};
#endif /* GLOBALSTARTGEOMETRY_H_ */
