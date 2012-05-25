/*
 * GlobalStartGeometry.cpp
 *
 *  Created on: 12.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 */


#include"GlobalStartGeometry.h"



const unsigned PRECISION = 5; // the number of times the loops are run determining the fill array _fill[i][j][k][l]

extern double LATTICE_CONST_WALL_LJTS;
//extern const double LATTICE_CONST_CU;

//@brief: zero temperature lattice constant of the LJ model used for copper, needed to build up the wall
//@todo: more robust way necessary: example of a problem: if sigma of the LJ wall is changed => Lattice Constant will change, too!
//@todo: CU_LC to be shifted to the class PhaseSpaceWriter
const double CU_LC = 6.83325;

GlobalStartGeometry::GlobalStartGeometry(unsigned in_nFluid, double in_rhoLiq, double in_rhoVap, double in_alpha, double in_beta, double in_gamma)
{
	_nFluid = in_nFluid;
	_rhoLiq = in_rhoLiq;
	_rhoVap = in_rhoVap;
	_alpha = in_alpha;
	_beta = in_beta;
	_gamma = in_gamma;
	_nLiq = _nFluid / (1+ (_alpha*_beta*_gamma-1.0) *_rhoVap/_rhoLiq );
	_nVap = _nFluid - _nLiq;
	cout << "rhoVap = "<< _rhoVap << "\t rhoLiq = " << _rhoLiq << "\n";
	cout << "N liquid: " << _nLiq << "\n";
	cout << "N vapor: " << _nVap << "\n";
}

GlobalStartGeometry::~GlobalStartGeometry(){

}

void GlobalStartGeometry::calculateBoxFluidOffset(double hWall, double shielding)
{
	double aApprox;
	unsigned nEkX;	// number of elementary lattice boxes of the solid wall in x-direction
	unsigned nEkZ;	// "		"			"			"		"		"	 "	z-direction
	double effFluid[3];
	double effLiq[3];

	// @brief: 2nd step: determining the simulation box, offset and fluid cuboid dimensions
	aApprox = pow(_nLiq/_rhoLiq, 1.0/3.0);
	nEkX = round(_alpha*aApprox / LATTICE_CONST_WALL_LJTS); 	// the number of lattices in the x-direction (Ek == elementary cristal)
	nEkZ = ceil(_gamma*aApprox /  LATTICE_CONST_WALL_LJTS);	// the number of lattices in the z-direction (Ek == elementary cristal)
	//cout << "nEKx = "<< nEkX << " \t nEkz = "<< nEkZ << "\n";


	_box[0] = nEkX * LATTICE_CONST_WALL_LJTS;
	_box[2] = nEkZ * LATTICE_CONST_WALL_LJTS;

	effLiq[0] = _box[0] / _alpha;
	effLiq[2] = _box[2] / _gamma;
	effLiq[1] = (_nLiq / _rhoLiq) / (effLiq[0]*effLiq[2]);
/*	cout << "effLiq[0]: " << effLiq[0] << "\n";
	cout << "effLiq[1]: " << effLiq[1] << "\n";
	cout << "effLiq[2]: " << effLiq[2] << "\n";
*/

	_box[1] = _nVap/_rhoVap/_box[0]/_box[2] + (effLiq[0]*effLiq[1]*effLiq[2])/_box[0]/_box[2] + hWall + shielding;
	cout << "hWall: " << hWall << "\n";
	cout << "shielding: " << shielding <<"\n";
	cout << "box[0]: " << _box[0] << "\n";
	cout << "box[1]: " << _box[1] << "\n";
	cout << "box[2]: " << _box[2] << "\n";

	effFluid[0] = pow( _nFluid/_rhoLiq , 1.0/3.0 );
	effFluid[1] = sqrt(_nFluid/_rhoLiq/effFluid[0]);
	effFluid[2] = _nFluid / (_rhoLiq*effFluid[0]*effFluid[1]);
/*	cout << "effFluid[0]: " << effFluid[0]<< "\n";
	cout << "effFluid[1]: " << effFluid[1]<< "\n";
	cout << "effFluid[2]: " << effFluid[2]<< "\n";
*/
	_off[0] = 0.5*(_box[0] - effFluid[0]);
	_off[1] = hWall + shielding + 2*LATTICE_CONST_WALL_LJTS;
	_off[2] = 0.5*(_box[2]-effFluid[2]);

	// @brief: 3rd step: setting up a cuboid containing all the fluid particles at liquid phase density
	double nFluidBoxes;
	nFluidBoxes = _nFluid/3;
	_fluidUnits[0] = round( pow( effFluid[0]*effFluid[0]*nFluidBoxes/effFluid[1]/effFluid[2], 1.0/3.0 ) );
	//double fluidUnitsYz = (double)nFluidBoxes/_fluidUnits[0];
	_fluidUnits[1] = round( sqrt(effFluid[1]*nFluidBoxes/_fluidUnits[0]/ effFluid[2]) );
	_fluidUnits[2] = ceil( nFluidBoxes/_fluidUnits[0]/_fluidUnits[1] );

	_fluidUnit[0] = (double)effFluid[0]/_fluidUnits[0];
	_fluidUnit[1] = (double)effFluid[1]/_fluidUnits[1];
	_fluidUnit[2] = (double)effFluid[2]/_fluidUnits[2];


	_fluidFillProbability = nFluidBoxes / ( _fluidUnits[0] * _fluidUnits[1] * _fluidUnits[2] );
}

//@todo: the tensor accounting for the filling of a single slot, i.e. _fill[][][][], should be part of an extra class that also determines the
// posistion of the particles(both fluid and wall)
void GlobalStartGeometry::calculateFillProbabilityArray()
{
	//@brief: initializing each element of _fill[][][][] by true
	_fill[_fluidUnits[0]][_fluidUnits[1]][_fluidUnits[2]][3] = new bool;
	for(unsigned i = 0; i< _fluidUnits[0]; i++){
		for(unsigned j = 0; j< _fluidUnits[1]; j++){
			for(unsigned k = 0; k< _fluidUnits[2]; k++){
				for(unsigned l = 0; l < 3; l++){
					_fill[i][j][k][l] = true;
				}
			}
		}
	}

	// number of filled slots, initialized by the total number of slots (both filled and not filled)
	//unsigned nFilledSlots;
	// total number of slots
	double totalNSlots;
	// a measure of the probability by which an element of _fill is swapped
	double pSwap;
	// number of slots beeing ideally filled
	double nIdeallyFilled;
	// bool valued quantity checking if particles have to be added or removed respectively
	bool tSwap;

	_nFilledSlots = 3*_fluidUnits[0]*_fluidUnits[1]*_fluidUnits[2];
//	cout << "number of filled slots at the beginning of Gloablstartgeometry:" << _nFilledSlots <<"\n";
	totalNSlots = _nFilledSlots; // slots is and "always" will be the total number of slots
	nIdeallyFilled = _fluidFillProbability * totalNSlots;
	RandomNumber rdm;

	for(unsigned m = 0; m < PRECISION; m++){
		tSwap = (_nFilledSlots < nIdeallyFilled);
		pSwap = (nIdeallyFilled - (double)_nFilledSlots)/( (tSwap ? totalNSlots : 0.0) - _nFilledSlots );
		for(unsigned i=0; i < _fluidUnits[0]; i++){
			for(unsigned j=0; j < _fluidUnits[1]; j++){
				for(unsigned k=0; k < _fluidUnits[2]; k++){
					for(unsigned l=0; l < 3; l++){
						if(pSwap >= rdm.randNum()){
							if(_fill[i][j][k][l]) _nFilledSlots--;
							_fill[i][j][k][l] = tSwap;
							if(tSwap) _nFilledSlots++;
						}
					}
				}
			}
		}
	}
/*	cout << "Filling" << _nFilledSlots << " out of a total number of " << totalNSlots << "fluid slots. Ideally a number of "
		<< nIdeallyFilled << "fluid slots was to be filled.\n" << "fluidFillProbability = " << _fluidFillProbability << "\n";*/
}

bool GlobalStartGeometry::gFillArray(unsigned fluidUnits0, unsigned fluidUnits1, unsigned flunidUnits2, unsigned particleInElementaryBox){
	return _fill[fluidUnits0][fluidUnits1][flunidUnits2][particleInElementaryBox];
}

double GlobalStartGeometry::gBoxLength(unsigned direction){
	return _box[direction];
}

unsigned GlobalStartGeometry::gFluidUnits(unsigned direction){
	return _fluidUnits[direction];
}

double GlobalStartGeometry::gOffset(unsigned direction){
	return _off[direction];
}

double GlobalStartGeometry::gFluidUnit(unsigned direction){
	return _fluidUnit[direction];
}

unsigned GlobalStartGeometry::gNFilledSlots(){
	return _nFilledSlots;
}
