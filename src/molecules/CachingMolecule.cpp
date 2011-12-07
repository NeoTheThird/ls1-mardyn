/***************************************************************************
 *   Copyright (C) 2009 by Martin Bernreuther and colleagues               *
 *   Martin.Bernreuther@informatik.uni-stuttgart.de                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   Due to copyleft all future versions of this program must be           *
 *   distributed as Free Software (e.g., using a BSD-like license).        *
 ***************************************************************************/
#include "molecules/CachingMolecule.h"


#include <cmath>
#include <fstream>
#include <cassert>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;


CachingMolecule::CachingMolecule(unsigned long id, unsigned int componentid,
	                 double rx, double ry, double rz,
	                 double vx, double vy, double vz,
	                 double q0, double q1, double q2, double q3,
	                 double Dx, double Dy, double Dz,
	                 const vector<Component>* components)
		: _q(q0, q1, q2, q3), _ljcenters(NULL), _tersoff(NULL) {
	_id = id;
	_componentid = componentid;
	_r[0] = rx;
	_r[1] = ry;
	_r[2] = rz;
	_v[0] = vx;
	_v[1] = vy;
	_v[2] = vz;
	_D[0] = Dx;
	_D[1] = Dy;
	_D[2] = Dz;
	_sites_d = _sites_F =_osites_e = NULL;
#ifdef TERSOFF_SUPPORT
	_numTersoffNeighbours = 0;
#endif
	fixedx = rx;
	fixedy = ry;
	if (components)
		setupCache(components);
}

CachingMolecule::CachingMolecule(const CachingMolecule& m) {
	_id = m._id;
	_componentid = m._componentid;
	_r[0] = m._r[0];
	_r[1] = m._r[1];
	_r[2] = m._r[2];
	_v[0] = m._v[0];
	_v[1] = m._v[1];
	_v[2] = m._v[2];
	_q = m._q;
	_D[0] = m._D[0];
	_D[1] = m._D[1];
	_D[2] = m._D[2];
	_F[0] = m._F[0];
	_F[1] = m._F[1];
	_F[2] = m._F[2];
	_M[0] = m._M[0];
	_M[1] = m._M[1];
	_M[2] = m._M[2];

	_ljcenters = m._ljcenters;
	_charges = m._charges;
	_dipoles = m._dipoles;
	_quadrupoles = m._quadrupoles;
	if (!m._tersoff) {
		global_log->warning() << "Tersoff vector null pointer detected for CachingMolecule " << _id << endl;
	}
	_tersoff = m._tersoff;
	_m = m._m;
	_I[0] = m._I[0];
	_I[1] = m._I[1];
	_I[2] = m._I[2];
	_invI[0] = m._invI[0];
	_invI[1] = m._invI[1];
	_invI[2] = m._invI[2];

	_numsites = m._numsites;
	_numorientedsites = m._numorientedsites;
	assert(_numsites);
	_sites_d = new double[_numsites*3];
	assert(_sites_d);
	//for(unsigned int i=0;i<_numsites*3;++i) _sites_d[i]=m._sites_d[i]; // not necessary -> cache only
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[numLJcenters()*3]);
	_dipoles_d = &(_charges_d[numCharges()*3]);
	_quadrupoles_d = &(_dipoles_d[numDipoles()*3]);
	_tersoff_d = &(_quadrupoles_d[numQuadrupoles()*3]);

	_osites_e = new double[_numorientedsites*3];
	assert(_osites_e);
	//for(unsigned int i=0;i<_numorientedsites*3;++i) _osites_e[i]=m._osites_e[i]; // not necessary -> cache only
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[numDipoles()*3]);

	_sites_F = new double[_numsites*3];

	assert(_sites_F);
	//for(unsigned int i=0;i<_numsites*3;++i) _sites_F[i]=m._sites_F[i]; // not necessary -> cache only
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[numLJcenters()*3]);
	_dipoles_F = &(_charges_F[numCharges()*3]);
	_quadrupoles_F = &(_dipoles_F[numDipoles()*3]);
	_tersoff_F = &(_quadrupoles_F[numQuadrupoles()*3]);
#ifdef TERSOFF_SUPPORT
	_numTersoffNeighbours = 0;
#endif
	fixedx = m.fixedx;
	fixedy = m.fixedy;
}


CachingMolecule& CachingMolecule::operator=(const CachingMolecule& rhs) {
	_id = rhs._id;
	_componentid = rhs._componentid;
	_r[0] = rhs._r[0];
	_r[1] = rhs._r[1];
	_r[2] = rhs._r[2];
	_v[0] = rhs._v[0];
	_v[1] = rhs._v[1];
	_v[2] = rhs._v[2];
	_q = rhs._q;
	_D[0] = rhs._D[0];
	_D[1] = rhs._D[1];
	_D[2] = rhs._D[2];
	_F[0] = rhs._F[0];
	_F[1] = rhs._F[1];
	_F[2] = rhs._F[2];
	_M[0] = rhs._M[0];
	_M[1] = rhs._M[1];
	_M[2] = rhs._M[2];

	_ljcenters = rhs._ljcenters;
	_charges = rhs._charges;
	_dipoles = rhs._dipoles;
	_quadrupoles = rhs._quadrupoles;
	if (!rhs._tersoff) {
		global_log->warning() << "Tersoff vector null pointer detected for CachingMolecule " << _id << endl;
	}
	_tersoff = rhs._tersoff;
	_m = rhs._m;
	_I[0] = rhs._I[0];
	_I[1] = rhs._I[1];
	_I[2] = rhs._I[2];
	_invI[0] = rhs._invI[0];
	_invI[1] = rhs._invI[1];
	_invI[2] = rhs._invI[2];

	_numsites = rhs._numsites;
	_numorientedsites = rhs._numorientedsites;
	assert(_numsites);

	if (_sites_d) {
		delete[] _sites_d;
	}

	_sites_d = new double[_numsites*3];
	assert(_sites_d);
	//for(unsigned int i=0;i<_numsites*3;++i) _sites_d[i]=m._sites_d[i]; // not necessary -> cache only
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[numLJcenters()*3]);
	_dipoles_d = &(_charges_d[numCharges()*3]);
	_quadrupoles_d = &(_dipoles_d[numDipoles()*3]);
	_tersoff_d = &(_quadrupoles_d[numQuadrupoles()*3]);

	if (_osites_e) {
		delete[] _osites_e;
	}

	_osites_e = new double[_numorientedsites*3];
	assert(_osites_e);
	//for(unsigned int i=0;i<_numorientedsites*3;++i) _osites_e[i]=m._osites_e[i]; // not necessary -> cache only
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[numDipoles()*3]);

	if (_sites_F) {
		delete[] _sites_F;
	}

	_sites_F = new double[_numsites*3];
	assert(_sites_F);
	//for(unsigned int i=0;i<_numsites*3;++i) _sites_F[i]=m._sites_F[i]; // not necessary -> cache only
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[numLJcenters()*3]);
	_dipoles_F = &(_charges_F[numCharges()*3]);
	_quadrupoles_F = &(_dipoles_F[numDipoles()*3]);
	_tersoff_F = &(_quadrupoles_F[numQuadrupoles()*3]);
#ifdef TERSOFF_SUPPORT
	_numTersoffNeighbours = 0;
#endif
	fixedx = rhs.fixedx;
	fixedy = rhs.fixedy;
	return *this;
}



void CachingMolecule::upd_preF(double dt, double vcorr, double Dcorr) {
	assert(_m);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / _m;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] = vcorr * _v[d] + dtInv2m * _F[d];
		_r[d] += dt * _v[d];
	}

	double w[3];
	_q.rotate(_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qhalfstep;
	_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		_D[d] = Dcorr * _D[d] + dt_halve * _M[d];
	qhalfstep.rotate(_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	_q.add(qincr);
	qcorr = 1. / sqrt(_q.magnitude2());
	_q.scale(qcorr);

}

void CachingMolecule::upd_cache() {
	// update Cache (rotate sites and save relative positions)
	unsigned int i;
	unsigned int ns = numLJcenters();
	for (i = 0; i < ns; ++i)
		_q.rotateinv((*_ljcenters)[i].r(), &(_ljcenters_d[i*3]));
	ns = numCharges();
	for (i = 0; i < ns; ++i)
		_q.rotateinv((*_charges)[i].r(), &(_charges_d[i*3]));
	ns = numDipoles();

	for (i = 0; i < ns; ++i) {
		const Dipole& di = (*_dipoles)[i];
		_q.rotateinv(di.r(), &(_dipoles_d[i*3]));
		_q.rotateinv(di.e(), &(_dipoles_e[i*3]));
	}
	ns = numQuadrupoles();

	for (i = 0; i < ns; ++i) {
		const Quadrupole& qi = (*_quadrupoles)[i];
		_q.rotateinv(qi.r(), &(_quadrupoles_d[i*3]));
		_q.rotateinv(qi.e(), &(_quadrupoles_e[i*3]));
	}
	ns = this->numTersoff();
	for (i = 0; i < ns; i++)
		_q.rotateinv((*_tersoff)[i].r(), &(_tersoff_d[i*3]));

	clearFM();
	//_Upot=0.;
}


void CachingMolecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {

	calcFM();

	double dtInv2m = dt_halve / _m;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		v2 += _v[d] * _v[d];
		_D[d] += dt_halve * _M[d];
	}
    assert(!isnan(v2)); // catches NaN
    summv2 += _m * v2;

	double w[3];
	_q.rotate(_D, w); // L = D = Iw
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
    assert(!isnan(Iw2)); // catches NaN
	sumIw2 += Iw2;
}

double CachingMolecule::U_rot() {
	double w[3];
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	return 0.5 * Iw2;
}

void CachingMolecule::calculate_mv2_Iw2(double& summv2, double& sumIw2) {
	summv2 += _m * v2();
	double w[3];
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void CachingMolecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {
	double vcx = _v[0] - offx;
	double vcy = _v[1] - offy;
	double vcz = _v[2] - offz;
	summv2 += _m * (vcx*vcx + vcy*vcy + vcz*vcz);

	double w[3];
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void CachingMolecule::scale_v(double s, double offx, double offy, double offz) {
	this->vsub(offx, offy, offz);
	this->scale_v(s);
	this->vadd(offx, offy, offz);
}

void CachingMolecule::write(ostream& ostrm) const {
	ostrm << _id << "\t" << (_componentid+1) << "\t"
	      << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
	      << _v[0] << " " << _v[1] << " " << _v[2] << "\t"
	      << _q.qw() << " " << _q.qx() << " " << _q.qy() << " " << _q.qz() << "\t"
	      << _D[0] << " " << _D[1] << " " << _D[2] << "\t"
	      << endl;
}

void CachingMolecule::addTersoffNeighbour(CachingMolecule* m, bool pairType) {
#ifdef TERSOFF_SUPPORT
	// this->_Tersoff_neighbours.insert(pair<CachingMolecule*, bool>(m, (pairType > 0)));
	for (int j = 0; j < _numTersoffNeighbours; j++) {
		if (m->_id == _Tersoff_neighbours_first[j]->id()) {
			this->_Tersoff_neighbours_first[j] = m;
			this->_Tersoff_neighbours_second[j] = pairType;
			return;
		}
	}

	this->_Tersoff_neighbours_first[_numTersoffNeighbours] = m;
	this->_Tersoff_neighbours_second[_numTersoffNeighbours] = pairType;
	this->_numTersoffNeighbours++;
	if (_numTersoffNeighbours > MAX_TERSOFF_NEIGHBOURS) {
		global_log->error() << "Tersoff neighbour list overflow: CachingMolecule " << m->_id << " has more than " << MAX_TERSOFF_NEIGHBOURS << " Tersoff neighbours." << endl;
		exit(1);
	}
#else
	std::cout << "CachingMolecule.cpp:382 Tersoff not supported!" << std::endl;
	exit(-1);
#endif
}

double CachingMolecule::tersoffParameters(double params[15]) //returns delta_r
{
	const Tersoff* t = &((*(this->_tersoff))[0]);
	params[ 0] = t->R();
	params[ 1] = t->S();
	params[ 2] = t->h();
	params[ 3] = t->cSquare();
	params[ 4] = t->dSquare();
	params[ 5] = t->A();
	params[ 6] = t->minusLambda();
	params[ 7] = t->minusMu();
	params[ 8] = t->beta();
	params[ 9] = t->n();
	params[10] = M_PI / (t->S() - t->R());
	params[11] = 1.0 + t->cSquare() / t->dSquare();
	params[12] = t->S() * t->S();
	params[13] = -(t->B());
	params[14] = -0.5 / t->n();

	return 0.000001 * (t->S() - t->R());
}

// private functions
// these are only used when compiling CachingMolecule.cpp and therefore might be inlined without any problems

inline void CachingMolecule::setupCache(const vector<Component>* components) {
	assert(components);
	if (components->size() == 0)
		return;
	_numsites = _numorientedsites = 0;
	_ljcenters = &(*components)[_componentid].ljcenters();
	_numsites += _ljcenters->size();
	_charges = &(*components)[_componentid].charges();
	_numsites += _charges->size();
	_dipoles = &(*components)[_componentid].dipoles();
	_numsites += _dipoles->size();
	_numorientedsites += _dipoles->size();
	_quadrupoles = &(*components)[_componentid].quadrupoles();
	_numsites += _quadrupoles->size();
	_numorientedsites += _quadrupoles->size();
	_tersoff = &(*components)[_componentid].tersoff();
#ifndef NDEBUG
	if (!_tersoff) {
		global_log->error() << "Tersoff vector null pointer detected for CachingMolecule " << _id << endl;
		exit(1);
	}
#endif
	_numsites += _tersoff->size();

	_m = (*components)[_componentid].m();
	_I[0] = (*components)[_componentid].I11();
	_I[1] = (*components)[_componentid].I22();
	_I[2] = (*components)[_componentid].I33();
	for (unsigned short d = 0; d < 3; ++d) {
		if (_I[d] != 0.)
			_invI[d] = 1. / _I[d];
		else
			_invI[d] = 0.;
	}

	assert(_numsites);

	_sites_d = new double[_numsites*3];
	assert(_sites_d);
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[numLJcenters()*3]);
	_dipoles_d = &(_charges_d[numCharges()*3]);
	_quadrupoles_d = &(_dipoles_d[numDipoles()*3]);
	_tersoff_d = &(_quadrupoles_d[numQuadrupoles()*3]);

	_osites_e = new double[_numorientedsites*3];
	assert(_osites_e);
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[numDipoles()*3]);

	_sites_F = new double[_numsites*3];

	assert(_sites_F);
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[numLJcenters()*3]);
	_dipoles_F = &(_charges_F[numCharges()*3]);
	_quadrupoles_F = &(_dipoles_F[numDipoles()*3]);
	_tersoff_F = &(_quadrupoles_F[numQuadrupoles()*3]);

	this->clearFM();
}

void CachingMolecule::clearFM() {
	for (unsigned int i = 0; i < _numsites * 3; ++i)
		_sites_F[i] = 0.;
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
}

void CachingMolecule::calcFM() {
	unsigned int ns = numSites();
	for (unsigned int si = 0; si < ns; ++si) {
		const double* Fsite = site_F(si);
		const double* dsite = site_d(si);
#ifndef NDEBUG
		/*
		 * catches NaN assignments
		 */
		for (int d = 0; d < 3; d++) {
			if (isnan(dsite[d])) {
				global_log->error() << "Severe dsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
			if (isnan(Fsite[d])) {
				global_log->error() << "Severe Fsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
		}
#endif
		Fadd(Fsite);
		_M[0] += dsite[1] * Fsite[2] - dsite[2] * Fsite[1];
		_M[1] += dsite[2] * Fsite[0] - dsite[0] * Fsite[2];
		_M[2] += dsite[0] * Fsite[1] - dsite[1] * Fsite[0];
	}
}


/*
 * catches NaN values and missing data
 *
 * @note Use isnan from cmath to check for nan.
 * If that's not available (C99), compare the value with itself. If the value
 * is NaN, the comparison will evaluate to false (according to IEEE754 spec.)
 */
void CachingMolecule::check(unsigned long id) {
	assert(_id == id);
	assert(_m > 0.0);
	assert(_numsites > 0);
	assert(_numorientedsites >= 0);
	for (int d = 0; d < 3; d++) {
		assert(!isnan(_r[d]));
		assert(!isnan(_v[d]));
		assert(!isnan(_D[d]));
		assert(!isnan(_F[d]));
		assert(!isnan(_M[d]));
		assert(!isnan(_I[d]));
		assert(!isnan(_invI[d]));
	}
}

bool CachingMolecule::isLessThan(const CachingMolecule& m2) const {
	if (_r[2] < m2.r(2))
		return true;
	else if (_r[2] > m2.r(2))
		return false;
	else {
		if (_r[1] < m2.r(1))
			return true;
		else if (_r[1] > m2.r(1))
			return false;
		else {
			if (_r[0] < m2.r(0))
				return true;
			else if (_r[0] > m2.r(0))
				return false;
			else {
				std::stringstream msg;
				msg << "LinkedCells::isFirstParticle: both Particles have the same position" << endl;
				msg << "m1:" << endl;
				this->write(msg);
				msg << "m2:" << endl;
				m2.write(msg);
				global_log->error() << msg.str();
				exit(1);
			}
		}
	}
}




unsigned long CachingMolecule::totalMemsize() const {
	unsigned long size = sizeof (*this);

	//_sites_d
	size += sizeof(double) * _numsites * 3;
	// site orientation _osites_e
	size += sizeof(double) * _numorientedsites * 3;
	// site Forces _sites_F
	// row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
	size += sizeof(double) * _numsites * 3;

	return size;
}


std::ostream& operator<<( std::ostream& os, const CachingMolecule& m ) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n" ;
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy() << ", " << m.q().qz()<< "]]\n" ;
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")" ;
	return os;
}
