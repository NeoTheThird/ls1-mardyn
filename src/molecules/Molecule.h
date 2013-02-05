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

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <vector>
#include <iostream>
#include <cassert>

#include "molecules/Quaternion.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"
#include "molecules/Component.h"
#include "integrators/Integrator.h"

/*
 * maximal size of the Tersoff neighbour list
 */
#define MAX_TERSOFF_NEIGHBOURS 10

class Domain;

//! @brief Molecule modeled as LJ sphere with point polarities + Tersoff potential
class Molecule {

public:
	// TODO Correct this constructor: the components vector is optional,
	// but if it is left away, all pointer data is not initialized (which is not
	// neccessarily bad), but then assertions fail (e.g. in the destructor) and we can't
	// use it's instances.
	Molecule(unsigned long id = 0, unsigned int componentid = 0,
	         double rx = 0., double ry = 0., double rz = 0.,
	         double vx = 0., double vy = 0., double vz = 0.,
	         double q0 = 0., double q1 = 0., double q2 = 0., double q3 = 0.,
	         double Dx = 0., double Dy = 0., double Dz = 0.,
	         const std::vector<Component>* components = NULL
	);
	Molecule(const Molecule& m);

private:
	Molecule& operator=(const Molecule& m);

public:
	~Molecule() {
		assert(_sites_d);
		delete[] _sites_d;
		assert(_osites_e);
		delete[] _osites_e;
		assert(_sites_F);
		delete[] _sites_F;
	}

	/** get molecule ID */
	unsigned long id() const { return _id; }
	/** set molecule ID */
	void setid(unsigned long id) { _id = id; }
	/** get the molecule's component ID */
	unsigned int componentid() const { return _componentid; }
	/** set the molecule's component ID */
	void setComponentid( unsigned int id ) { _componentid = id; }
	/** get position coordinate */
	double r(unsigned short d) const { return _r[d]; }
	/** set position coordinate */
	void setr(unsigned short d, double r) { _r[d] = r; }
	/** get velocity coordinate */
	double v(unsigned short d) const { return _v[d]; }
	/** set the velocity coordiate */
	void setv(unsigned short d, double vd) { _v[d] = vd; }
	/** get coordinate of current force onto molecule */
	double F(unsigned short d) const {return _F[d]; }
	/** get the orientation */
	const Quaternion& q() const { return _q; }
	/** get coordinate of the angular momentum */
	double D(unsigned short d) const { return _D[d]; } /* TODO: we should rename D to L with respect to literature. */
	/** set the coordinate of the angular momentum */
	void setD(unsigned short d, double Dd){ _D[d] = Dd; }
	/** get coordinate of the current angular momentum  onto molecule */ 
	double M(unsigned short d) const {return _M[d]; }


	inline void move(int d, double dr) { _r[d] += dr; } /* TODO: is this realy needed? */
	
	// by Stefab Becker <stefan.becker@mv.uni-kl.de> 
	// method returns the total mass of a particle
	double gMass(){return _m;}


	/** calculate and return the square velocity */
	double v2() const {return _v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]; }
	
	/** return the translational energy of the molecule */
	double U_trans() const { return 0.5 * _m * v2(); }
	/** return the rotational energy of the molecule */
	double U_rot();
	/** return total kinetic energy of the molecule */
	double U_kin() { return U_trans() + U_rot(); }
	
	/* TODO: Maybe we should better do this using the component directly? 
	 * In the GNU STL vector.size() causes two memory accesses and one subtraction!
	 */
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


	/* TODO: Is this realy necessary? Better use a function like dist2(m1.r(), m2.r()). */
	/** Calculate the difference vector and return the square (euclidean) distance.
	 *
	 *  \param molecule2 molecule to which the distance shall be calculated
	 */
	double dist2(const Molecule& molecule2, double dr[3]) const {
		double d2 = 0.;
		for (unsigned short d = 0; d < 3; d++) {
			dr[d] = molecule2._r[d] - _r[d];
			d2 += dr[d] * dr[d];
		}
		return d2;
	}
	

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

	
	
	/** First step of the leap frog integrator */
	void upd_preF(double dt, double vcorr=1., double Dcorr=1.);
	/** update the molecules site position caches (rotate sites and save relative positions) */
	void upd_cache();
	/** second step of the leap frog integrator */
	void upd_postF(double dt_halve, double& summv2, double& sumIw2);

	//! calculate summv2 and sumIw2
	//! @todo what is sumIw2?
	//! @todo comment
	void calculate_mv2_Iw2(double& summv2, double& sumIw2);
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz);

	/** write information to stream */
	void write(std::ostream& ostrm) const;

	inline unsigned getCurTN() { return this->_numTersoffNeighbours; }
	inline Molecule* getTersoffNeighbour(unsigned i) { return this->_Tersoff_neighbours_first[i]; }
	inline void clearTersoffNeighbourList() { this->_numTersoffNeighbours = 0; }
	void addTersoffNeighbour(Molecule* m, bool pairType);
	double tersoffParameters(double params[15]); //returns delta_r

	/** clear forces and moments */
	void clearFM();
	/** calculate forces and moments for already given site forces */
	void calcFM();
	
	/** perform data consistency check for the molecule (only debug mode) */
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
	bool isLessThan(const Molecule& m2) const;

private:

	unsigned long _id; 	/**< IDentification number of that molecule */
	unsigned int _componentid;  /**< IDentification number of its component type */
	double _r[3];  /**< position coordinates */
	double _F[3];  /**< forces */
	double _v[3];  /**< velocity */
	Quaternion _q; /**< angular orientation */
	double _M[3];  /**< torsional moment */
    // TODO: We should rename _D to _L with respect to the literature.
	double _D[3];  /**< angular momentum */

	/* Caches for component data */
	const std::vector<LJcenter>* _ljcenters;
	const std::vector<Charge>* _charges;
	const std::vector<Dipole>* _dipoles;
	const std::vector<Quadrupole>* _quadrupoles;
	const std::vector<Tersoff>* _tersoff;

	double _m; /**< total mass */
	double _I[3],_invI[3];  // moment of inertia for principal axes and it's inverse
	std::size_t _numsites; // number of sites
	std::size_t _numorientedsites; // number of oriented sites (subset of sites)
	// global site coordinates relative to site origin
	// row order: dx1,dy1,dz1,dx2,dy2,dz2,...
	/* TODO: Maybe change to absolute positions for many center molecules. */
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

	Molecule* _Tersoff_neighbours_first[MAX_TERSOFF_NEIGHBOURS];
	bool _Tersoff_neighbours_second[MAX_TERSOFF_NEIGHBOURS]; /* TODO: Comment */
	int _numTersoffNeighbours;
	double fixedx, fixedy;

	// setup cache values/properties
	void setupCache(const std::vector<Component>* components);

};


std::ostream& operator<<( std::ostream& os, const Molecule& m );

#endif /*MOLECULE_H_*/
