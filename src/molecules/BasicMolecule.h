/*************************************************************************
 * Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or (at *
 * your option) any later version.                                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            * 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the Free Software           *
 * Foundation, 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.   *
 *************************************************************************/

#ifndef BASICMOLECULE_H_
#define BASICMOLECULE_H_

/*
 * maximal size of the Tersoff neighbour list
 */
#define MAX_TERSOFF_NEIGHBOURS 10

#include "molecules/Quaternion.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"
#include "molecules/Component.h"
#include "integrators/Integrator.h"
class Domain;

#include <vector>
#include <iostream>
#include <cassert>

//! @brief Molecule modeled as LJ sphere with point polarities + Tersoff potential
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
class BasicMolecule {

public:
	// TODO Correct this constructor: the components vector is optional,
	// but if it is left away, all pointer data is not initialized (which is not
	// neccessarily bad), but then assertions fail (e.g. in the destructor) and we can't
	// use it's instances.
	BasicMolecule(unsigned long id = 0, int componentid = 0,
	         double rx = 0., double ry = 0., double rz = 0.,
	         double vx = 0., double vy = 0., double vz = 0.,
	         double q0 = 0., double q1 = 0., double q2 = 0., double q3 = 0.,
	         double Dx = 0., double Dy = 0., double Dz = 0.,
	         const std::vector<Component>* components = NULL
	);
	BasicMolecule(const BasicMolecule& m);

	~BasicMolecule() {
		assert(_sites_d); delete[] _sites_d;
		assert(_osites_e); delete[] _osites_e;
		assert(_sites_F); delete[] _sites_F;
	}

	BasicMolecule& operator=(const BasicMolecule& rhs);


	/** get the ID */
	unsigned long id() const { return _id; }
	void setid(unsigned long id) { this->_id = id; }
	/** get the Component */
	int componentid() const { return _componentid; }
	/** get the position */
	double r(unsigned short d) const { return _r[d]; }

	/** get the velocity */
	double v(unsigned short d) const { return _v[d]; }
	/** get the Orientation */
	const Quaternion& q() const { return _q; }

	inline void move(int d, double dr) { _r[d] += dr; }

	/** get the rotatational speed */
	double D(unsigned short d) const { return _D[d]; }

	/** get F */
	double F(unsigned short d) const {return _F[d]; }
	/** get M */
	double M(unsigned short d) const {return _M[d]; }

	double Utrans() const { return .5*_m*(_v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]); }
	double Urot();

	/** get number of sites */
	unsigned int numSites() const { return _numsites; }
	unsigned int numLJcenters() const { return _ljcenters->size(); }
	unsigned int numCharges() const { return _charges->size(); }
	unsigned int numDipoles() const { return _dipoles->size(); }
	unsigned int numQuadrupoles() const { return _quadrupoles->size(); }
	unsigned int numTersoff() const { assert(_tersoff); return _tersoff->size(); }

	const double* site_d(unsigned int i) const { return &(_sites_d[3*i]); }
	const double* site_F(unsigned int i) const { return &(_sites_F[3*i]); }
	const double* ljcenter_d(unsigned int i) const { return &(_ljcenters_d[3*i]); }
	const double* charge_d(unsigned int i) const { return &(_charges_d[3*i]); }
	const double* dipole_d(unsigned int i) const { return &(_dipoles_d[3*i]); }
	const double* dipole_e(unsigned int i) const { return &(_dipoles_e[3*i]); }
	const double* quadrupole_d(unsigned int i) const { return &(_quadrupoles_d[3*i]); }
	const double* quadrupole_e(unsigned int i) const { return &(_quadrupoles_e[3*i]); }

	/**
	 * get the total object memory size, together with all its members
	 * \Note You can retrieve the size of the molecule class itself simply
	 *       with the sizeof()-operator.
	 */
	unsigned long totalMemsize() const;

	/** set the position */
	void setr(unsigned short d, double r) { _r[d]=r; }

	/** Calculate the difference vector and return the square (euclidean) distance.
	 *
	 *  \param molecule2 molecule to which the distance shall be calculated
	 */
	double dist2(const BasicMolecule& molecule2, double dr[3]) const {
		double d2 = 0.;
		for (unsigned short d = 0; d < 3; d++) {
			dr[d] = molecule2._r[d] - _r[d];
			d2 += dr[d] * dr[d];
		}
		return d2;
	}
	
	/** calculate and return the square velocity */
	double v2() const {return _v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]; }

	/** set force acting on molecule
	 * @param F force vector (x,y,z)
	 */
	void setF(double F[3]) { for(int d = 0; d < 3; d++ ) { _F[d] = F[d]; } }

	/** set momentum acting on molecule 
	 * @param M force vector (x,y,z)
	 */
	void setM(double M[3]) { for(int d = 0; d < 3; d++ ) { _M[d] = M[d]; } }

	void scale_v(double s) { for(unsigned short d=0;d<3;++d) _v[d]*=s; }
	void scale_v(double s, double offx, double offy, double offz);
	void scale_F(double s) { for(unsigned short d=0;d<3;++d) _F[d]*=s; }
	void scale_D(double s) { for(unsigned short d=0;d<3;++d) _D[d]*=s; }
	void scale_M(double s) { for(unsigned short d=0;d<3;++d) _M[d]*=s; }

	void Fadd(const double a[]) { for(unsigned short d=0;d<3;++d) _F[d]+=a[d]; }

	void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) _M[d]+=a[d]; }

	void vadd(const double ax, const double ay, const double az) {
		_v[0] += ax; _v[1] += ay; _v[2] += az;
	}
	void vsub(const double ax, const double ay, const double az) {
		_v[0] -= ax; _v[1] -= ay; _v[2] -= az;
	}
	void setXY() { fixedx = _r[0]; fixedy = _r[1]; }
	void resetXY()
	{
		_v[0] = 0.0;
		_v[1] = 0.0;
		_F[1] = 0.0;
		_F[0] = 0.0;
		_r[0] = fixedx;
		_r[1] = fixedy;
	}

	void Fljcenteradd(unsigned int i, double a[])
	{ double* Fsite=&(_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fljcentersub(unsigned int i, double a[])
	{ double* Fsite=&(_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fchargeadd(unsigned int i, double a[])
	{ double* Fsite=&(_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fchargesub(unsigned int i, double a[])
	{ double* Fsite=&(_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fdipoleadd(unsigned int i, double a[])
	{ double* Fsite=&(_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fdipolesub(unsigned int i, double a[])
	{ double* Fsite=&(_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fquadrupoleadd(unsigned int i, double a[])
	{ double* Fsite=&(_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fquadrupolesub(unsigned int i, double a[])
	{ double* Fsite=&(_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Ftersoffadd(unsigned int i, double a[])
	{ double* Fsite=&(_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }

	void upd_preF(double dt, double vcorr=1., double Dcorr=1.);
	void upd_cache();
	void upd_postF(double dt_halve, double& summv2, double& sumIw2);

	//! calculate summv2 and sumIw2
	//! @todo what is sumIw2?
	//! @todo comment
	void calculate_mv2_Iw2(double& summv2, double& sumIw2);
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz);

	/** write information to stream */
	void write(std::ostream& ostrm) const;

	inline unsigned getCurTN() { return this->_numTersoffNeighbours; }
	inline BasicMolecule* getTersoffNeighbour(unsigned i) { return this->_Tersoff_neighbours_first[i]; }
	inline void clearTersoffNeighbourList() { this->_numTersoffNeighbours = 0; }
	void addTersoffNeighbour(BasicMolecule* m, bool pairType);
	double tersoffParameters(double params[15]); //returns delta_r

	// clear forces and moments
	void clearFM();
	// calculate forces and moments for already given site forces
	void calcFM();
	void check(unsigned long id);

	//! @brief find out whether m1 is before m2 (in some global ordering)
	//!
	//! Compares this molecule to m2 based on their coordinates.
	//!
	//! @return true if this molecule is smaller than m2 (according to the order
	//!         induced by the coordinates (z,y,x)
	//!
	//! At the boundary between two processes (if used in parallel mode), the forces
	//! for pairs which cross the boundary are calculated twice (once by each proc who
	//! owns one of the particles). But the contribution to macroscopic value must be
	//! counted only once, which is done by the process who owns the "first" particle.
	//! As order criterion, the spacial position is used int this method. The particles
	//! with lower x-coordinate is first (if equal, then y- or z-coordinate).
	//! For pairs which are completely on one process, the first particle can be
	//! determined from the cell structure. But for pairs on different procs, the
	//! corresponding cell discretisations might be different as well, and therefore
	//! the cell structure must not be used to determine the order.
	bool isLessThan(const BasicMolecule& m2) const;

private:

	unsigned long _id; // IDentification number of that molecule
	int _componentid;  // IDentification number of its component type
	double _r[3];  // position coordinates
	double _F[3];  // forces
	double _v[3];  // velocity
	Quaternion _q; // angular orientation
	double _M[3];  // torsional moment
    // TODO: We should rename _D to _L with respect to the literature.
	double _D[3];  // angular momentum 

	const std::vector<LJcenter>* _ljcenters;
	const std::vector<Charge>* _charges;
	const std::vector<Dipole>* _dipoles;
	const std::vector<Quadrupole>* _quadrupoles;
	const std::vector<Tersoff>* _tersoff;

	double _m; // total mass
	double _I[3],_invI[3];  // moment of inertia for principal axes and it's inverse
	std::size_t _numsites; // number of sites
	std::size_t _numorientedsites; // number of oriented sites (subset of sites)
	// global site coordinates relative to site origin
	// row order: dx1,dy1,dz1,dx2,dy2,dz2,...
	double *_sites_d;
	double *_ljcenters_d, *_charges_d, *_dipoles_d,
	       *_quadrupoles_d, *_tersoff_d;
	// site orientation
	double *_osites_e;
	double *_dipoles_e, *_quadrupoles_e;
	// site Forces
	// row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
	double* _sites_F;
	double *_ljcenters_F, *_charges_F, *_dipoles_F,
	       *_quadrupoles_F, *_tersoff_F;

	BasicMolecule* _Tersoff_neighbours_first[MAX_TERSOFF_NEIGHBOURS];
	bool _Tersoff_neighbours_second[MAX_TERSOFF_NEIGHBOURS];
	int _numTersoffNeighbours;
	double fixedx, fixedy;

	// setup cache values/properties
	void setupCache(const std::vector<Component>* components);

};

#endif /*BASICMOLECULE_H_*/
