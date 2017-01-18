
#ifndef HOMOGENEOUS_H__
#define HOMOGENEOUS_H__

#include "LongRangeCorrection.h"

#include <cmath>

class Domain;
class LongRangeCorrection;

//class Homogeneous:public LongRangeCorrection{
class Homogeneous: public LongRangeCorrection{

public:
//	Homogeneous();
	Homogeneous(double cutoffRadius, double cutoffRadiusLJ,  Domain* domain, Simulation _simulation);
  
//	void initializeLongRange();
	void calculateLongRange();

private:
	double _UpotCorr;
	double _VirialCorr;
	/* TODO: Comments on all the functions */
	// Long range correction for the Lennard-Jones interactions based on Lustig (1988)
	double _TICCu(int n,double rc,double sigma2);
	double _TICSu(int n,double rc,double sigma2,double tau);
	double _TISSu(int n,double rc,double sigma2,double tau1,double tau2);
	double _TICCv(int n,double rc,double sigma2);
	double _TICSv(int n,double rc,double sigma2,double tau);
	double _TISSv(int n,double rc,double sigma2,double tau1,double tau2);
	
	//! Components resp. molecule types
	std::vector<Component> _components;
	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;
	
	Domain* _domain;
};

#endif /* __HOMOGENEOUS_H__ */
