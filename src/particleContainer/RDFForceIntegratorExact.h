/*
 * RDFForceIntegratorExact.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATOREXACT_H_
#define RDFFORCEINTEGRATOREXACT_H_

#include "RDFForceIntegrator.h"
#include "molecules/potforce.h"
#include "RDF.h"

class RDFForceIntegratorExact: public RDFForceIntegrator {
public:
	RDFForceIntegratorExact(ParticleContainer* moleculeContainer, double rc,
			std::vector<std::vector<double> >* globalADist, std::vector<
					std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorExact();

	double traverseMolecules();

	void getScalingFactor(double* mol_r, double* site_r,
			double x, double y, double z, int site_i, double* scale);

	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false);


	void prepareUnitTest(Molecule* m) {
		_dx = _dy = _dz = _d_level = _dn = _dr = 10;
		_d_alpha = 180;
		_rho = 1;
		_extension = m->ljcenter_disp(0);

	}
	static double* precomputeScalingFactorsX(bool unit_test = false);
private:
	static double _dx, _dy, _dz, _extension, _dn, _dr, *_scaling_factors_x,  _d_alpha, _d_level, _rho, _g_start;
	static int _n_r, _n_n, _n_levels, _n_alpha;
	static bool called_x;
	static int timestep;
	static double* rhos;

	double integrateRDFSiteCartesian(double xlim[2], double ylim[2],
			double zlim[2], Molecule* mol, int plane, unsigned int site,
			int boundary[3], bool add_influence, double* return_force, bool unit_test = false);

	double checkScalingFactor(int idx_level,
			int idx_n, int idx_r);

};

#endif /* RDFFORCEINTEGRATOREXACT_H_ */
