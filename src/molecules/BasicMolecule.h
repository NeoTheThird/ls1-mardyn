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

//#define TERSOFF_SUPPORT

/*
 * maximal size of the Tersoff neighbour list
 */
#define MAX_TERSOFF_NEIGHBOURS 10

#include "molecules/Quaternion.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"
#include "molecules/Component.h"
#include "integrators/Integrator.h"
#include "utils/Logger.h"
#include "utils/DynamicArray.h"
class Domain;

#include <vector>
#include <iostream>
#include <cassert>

//! @brief Molecule modeled as LJ sphere with point polarities + Tersoff potential
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
class BasicMolecule {

public:

	typedef float fp_type;

	friend class ParticleData;

	// TODO Correct this constructor: the components vector is optional,
	// but if it is left away, all pointer data is not initialized (which is not
	// neccessarily bad), but then assertions fail (e.g. in the destructor) and we can't
	// use it's instances.
	BasicMolecule(unsigned long id = 0, unsigned int componentid = 0,
	         double rx = 0., double ry = 0., double rz = 0.,
	         double vx = 0., double vy = 0., double vz = 0.,
	         double q0 = 0., double q1 = 0., double q2 = 0., double q3 = 0.,
	         double Dx = 0., double Dy = 0., double Dz = 0.,
	         const std::vector<Component>* components = NULL
	);
	BasicMolecule(const BasicMolecule& m);

	~BasicMolecule() {
	}

	BasicMolecule& operator=(const BasicMolecule& rhs);


	/** get the ID */
	unsigned long id() const { return _id; }
	void setid(unsigned long id) { this->_id = id; }
	/** get the Component */
	unsigned int componentid() const { return _componentid; }
	/** get the position */
	double r(unsigned short d) const { return _r[d]; }

	/** get the velocity */
	double v(unsigned short d) const { return _v[d]; }
	/** get the Orientation */
	const Quaternion<fp_type>& q() const { return _q; }

	inline void move(int d, double dr) { _r[d] += dr; }

	/** get the rotatational speed */
	double D(unsigned short d) const { return _D[d]; }

	/** get F */
	double F(unsigned short d) const {return _F[d]; }
	/** get M */
	double M(unsigned short d) const {return _M[d]; }

	double U_trans() const {
		double m = (*components)[_componentid].m();
		return 0.5 * m * v2();
	}
	/** return the rotational energy of the molecule */
	double U_rot();
	/** return total kinetic energy of the molecule */
	double U_kin() { return U_trans() + U_rot(); }

	/** get number of sites */
	unsigned int numSites() const { return (*components)[_componentid].numSites(); }
	unsigned int numLJcenters() const { return (*components)[_componentid].numLJcenters(); }
	unsigned int numCharges() const { return (*components)[_componentid].numCharges(); }
	unsigned int numDipoles() const { return (*components)[_componentid].numDipoles(); }
	unsigned int numQuadrupoles() const { return (*components)[_componentid].numQuadrupoles(); }
	unsigned int numTersoff() const {
		//std::cout << "NumTersoff cid=" << _componentid << " ntersoff="<<(*components)[_componentid].numTersoff() << std::endl;
		return (*components)[_componentid].numTersoff(); }

	void ljcenter_d(unsigned int i, double d[3]) const {
		assert(i < numLJcenters());
		_q.rotateinv((*components)[_componentid].ljcenter(i).r(), d);
	}

	void charge_d(unsigned int i, double d[3]) const {
		assert(i < numCharges());
		_q.rotateinv((*components)[_componentid].charge(i).r(), d);
	}

	void dipole_d(unsigned int i, double d[3]) const {
		assert(i < numDipoles());
		_q.rotateinv((*components)[_componentid].dipole(i).r(), d);
	}

	void dipole_e(unsigned int i, double d[3]) const {
		assert(i < numDipoles());
		_q.rotateinv((*components)[_componentid].dipole(i).e(), d);
	}

	void quadrupole_d(unsigned int i, double d[3]) const {
		assert(i < numQuadrupoles());
		_q.rotateinv((*components)[_componentid].quadrupole(i).r(), d);
	}

	void quadrupole_e(unsigned int i, double d[3]) const {
		assert(i < numQuadrupoles());
		_q.rotateinv((*components)[_componentid].quadrupole(i).e(), d);
	}

	void tersoff_d(unsigned int i, double d[3]) const {
		assert(i < numTersoff());
		_q.rotateinv((*components)[_componentid].tersoff(i).r(), d);
	}

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

	void Fadd(const double a[]) {
		for(unsigned short d=0;d<3;++d) {
			_F[d]+=a[d];
		}
	}

	void Fsub(const double a[]) {
		for(unsigned short d=0;d<3;++d) {
			_F[d]-=a[d];
		}
	}

	void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) _M[d]+=a[d]; }

	void vadd(const double ax, const double ay, const double az) {
		_v[0] += ax; _v[1] += ay; _v[2] += az;
	}
	void vsub(const double ax, const double ay, const double az) {
		_v[0] -= ax; _v[1] -= ay; _v[2] -= az;
	}

	/*	void setXY() { fixedx = _r[0]; fixedy = _r[1]; }
	void resetXY()
	{
		_v[0] = 0.0;
		_v[1] = 0.0;
		_F[1] = 0.0;
		_F[0] = 0.0;
		_r[0] = fixedx;
		_r[1] = fixedy;
	}
*/

	// add force f effective at distance r
	void Fadd(const double f[3], const double r[3]) {
		Fadd(f);
		_M[0] += r[1] * f[2] - r[2] * f[1];
		_M[1] += r[2] * f[0] - r[0] * f[2];
		_M[2] += r[0] * f[1] - r[1] * f[0];
	}

	// substrace force f effective at distance r
	void Fsub(const double f[3], const double r[3]) {
		Fsub(f);
		_M[0] -= r[1] * f[2] - r[2] * f[1];
		_M[1] -= r[2] * f[0] - r[0] * f[2];
		_M[2] -= r[0] * f[1] - r[1] * f[0];
	}

	void Fljcenteradd(unsigned int i, double a[]) {
		double r[3];
		ljcenter_d(i, r);
		Fadd(a, r);
	}

	void Fljcentersub(unsigned int i, double a[]) {
		double r[3];
		ljcenter_d(i, r);
		Fsub(a, r);
	}

	void Fchargeadd(unsigned int i, double a[]){
		double r[3];
		charge_d(i, r);
		Fadd(a, r);
	}

	void Fchargesub(unsigned int i, double a[]) {
		double r[3];
		charge_d(i, r);
		Fsub(a, r);
	}

	void Fdipoleadd(unsigned int i, double a[]) {
		double r[3];
		dipole_d(i, r);
		Fadd(a, r);
	}
	void Fdipolesub(unsigned int i, double a[]) {
		double r[3];
		dipole_d(i, r);
		Fsub(a, r);
	}

	void Fquadrupoleadd(unsigned int i, double a[]) {
		double r[3];
		quadrupole_d(i, r);
		Fadd(a, r);
	}

	void Fquadrupolesub(unsigned int i, double a[]) {
		double r[3];
		quadrupole_d(i, r);
		Fsub(a, r);
	}

	void Ftersoffadd(unsigned int i, double a[]) {
		double r[3];
		tersoff_d(i, r);
		Fadd(a, r);
	}

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

	inline unsigned getCurTN() {
		#ifdef TERSOFF_SUPPORT
			assert(_Tersoff_neighbours_first.size() == _Tersoff_neighbours_second.size());
			return _Tersoff_neighbours_first.size();
		#else
			warn_tersoff("getCurTN");
			return 0;
		#endif
	}

	inline BasicMolecule* getTersoffNeighbour(unsigned i) {
		#ifdef TERSOFF_SUPPORT
			assert(i < _Tersoff_neighbours_first.size());
			return this->_Tersoff_neighbours_first[i];
		#else
			warn_tersoff("getTersoffNeighbour"); return NULL;
		#endif
	}

	inline void clearTersoffNeighbourList() {
		#ifdef TERSOFF_SUPPORT
			_Tersoff_neighbours_first.clear();
			_Tersoff_neighbours_second.clear();
		#else
			warn_tersoff("clearTersoffNeighbourList");
		#endif
	}

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

	void toString(std::ostream& ostrm) const {
		ostrm << "ID=" << _id << " cid=" << _componentid << std::endl;
		ostrm << "r=" << _r[0] << "," << _r[1] << "," << _r[2] << " v=" << _v[0] << "," << _v[1] << "," << _v[2] << std::endl;
		ostrm << "F=" << _F[0] << "," << _F[1] << "," << _F[2] << " M=" << _M[0] << "," << _M[1] << "," << _M[2] << std::endl;
		ostrm << "D=" << _D[0] << "," << _D[1] << "," << _D[2] << "  q=" << _q.qw() << "," << _q.qx() << "," << _q.qy() << "," << _q.qz() << std::endl;
	}

	// public for testing purpose...
	static void setComponentsVector(std::vector<Component>* components) {
		BasicMolecule::components = components;
	}

private:
	static const std::vector<Component>* components;

	unsigned long int _id; // IDentification number of that molecule
	unsigned int _componentid;  // IDentification number of its component type
	fp_type _r[3];  // position coordinates
	fp_type _v[3];  // velocity
	Quaternion<fp_type> _q; // angular orientation
    // TODO: We should rename _D to _L with respect to the literature.
	fp_type _D[3];  // angular momentum; the angular velocity;
	fp_type _F[3];  // forces
	fp_type _M[3];  // torsional moment

	// memory: 5 x 3 x 4 + 1 x 4 x 4 + 1 x 4 + 1 x 8 = 76 Bytes


#ifdef TERSOFF_SUPPORT
	utils::DynamicArray<BasicMolecule*, true, false> _Tersoff_neighbours_first;
	utils::DynamicArray<bool, true, false> _Tersoff_neighbours_second;
#endif

//	double fixedx, fixedy;

	// setup cache values/properties
	void setupCache(const std::vector<Component>* components);

	void warn_tersoff(const std::string& method) {
#ifndef NDEBUG
		static std::vector<std::string> warnings;
		for (size_t i = 0; i < warnings.size(); i++) {
			if (warnings[i] == method) {
				//Log::global_log->info() << "Suppressing Warning " << method << std::endl;
				return;
			}
		}
		warnings.push_back(method);
		Log::global_log->warning() << "Method BasicMolecule::" << method << "!"<< std::endl;
		Log::global_log->warning() << "Compile with -DTERSOFF_SUPPORT to enable tersoff support in BasicMolecule!" << std::endl;
#endif
	}

};

std::ostream& operator<<( std::ostream& os, const BasicMolecule& m );

#endif /*BASICMOLECULE_H_*/
