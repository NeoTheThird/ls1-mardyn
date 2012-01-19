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
#include "molecules/BasicMolecule.h"


#include <cmath>
#include <fstream>
#include <cassert>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;

const std::vector<Component>* BasicMolecule::components = NULL;


BasicMolecule::BasicMolecule(unsigned long id, unsigned int componentid,
	                 double rx, double ry, double rz,
	                 double vx, double vy, double vz,
	                 double q0, double q1, double q2, double q3,
	                 double Dx, double Dy, double Dz,
	                 const vector<Component>* components)
		: _q(q0, q1, q2, q3) /* _ljcenters(NULL), _tersoff(NULL), _components(*components) */ {
	if (BasicMolecule::components == NULL) {
		BasicMolecule::components = components;
	}
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
//	_numTersoffNeighbours = 0;
//	fixedx = rx;
//	fixedy = ry;
	if (components)
		setupCache(components);
}

BasicMolecule::BasicMolecule(const BasicMolecule& m) {
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

/*	_numTersoffNeighbours = 0;
	fixedx = m.fixedx;
	fixedy = m.fixedy;
*/
}


BasicMolecule& BasicMolecule::operator=(const BasicMolecule& rhs) {
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

/*	_numTersoffNeighbours = 0;
	fixedx = rhs.fixedx;
	fixedy = rhs.fixedy;
*/	return *this;
}



void BasicMolecule::upd_preF(double dt, double vcorr, double Dcorr) {
	double m = (*components)[_componentid].m();
	assert(m);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / m;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] = vcorr * _v[d] + dtInv2m * _F[d];
		_r[d] += dt * _v[d];
	}

	fp_type w[3];
	_q.rotate(_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= (*components)[_componentid].invI(d);
	Quaternion<fp_type> qhalfstep;
	_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		_D[d] = Dcorr * _D[d] + dt_halve * _M[d];
	qhalfstep.rotate(_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= (*components)[_componentid].invI(d);
	Quaternion<fp_type> qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	_q.add(qincr);
	qcorr = 1. / sqrt(_q.magnitude2());
	_q.scale(qcorr);

}

void BasicMolecule::upd_cache() {
	clearFM();
}


void BasicMolecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {

	calcFM();

	double m = (*components)[_componentid].m();
	double dtInv2m = dt_halve / m;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		v2 += _v[d] * _v[d];
		_D[d] += dt_halve * _M[d];
	}
    assert(!isnan(v2)); // catches NaN
    summv2 += m * v2;

    fp_type w[3];
	_q.rotate(_D, w); // L = D = Iw
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= (*components)[_componentid].invI(d);
		Iw2 += (*components)[_componentid].I(d) * w[d] * w[d];
	}
    assert(!isnan(Iw2)); // catches NaN
	sumIw2 += Iw2;
}

double BasicMolecule::U_rot() {
	fp_type w[3];
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= (*components)[_componentid].invI(d);
		Iw2 += (*components)[_componentid].I(d) * w[d] * w[d];
	}
	return 0.5 * Iw2;
}

void BasicMolecule::calculate_mv2_Iw2(double& summv2, double& sumIw2) {
	double m = (*components)[_componentid].m();
	summv2 += m * v2();
	fp_type w[3];
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= (*components)[_componentid].invI(d);
		Iw2 += (*components)[_componentid].I(d) * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void BasicMolecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {
	double vcx = _v[0] - offx;
	double vcy = _v[1] - offy;
	double vcz = _v[2] - offz;
	double m = (*components)[_componentid].m();
	summv2 += m * (vcx*vcx + vcy*vcy + vcz*vcz);

	fp_type w[3];
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= (*components)[_componentid].invI(d);
		Iw2 += (*components)[_componentid].I(d) * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void BasicMolecule::scale_v(double s, double offx, double offy, double offz) {
	this->vsub(offx, offy, offz);
	this->scale_v(s);
	this->vadd(offx, offy, offz);
}

void BasicMolecule::write(ostream& ostrm) const {
	ostrm << _id << "\t" << (_componentid+1) << "\t"
	      << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
	      << _v[0] << " " << _v[1] << " " << _v[2] << "\t"
	      << _q.qw() << " " << _q.qx() << " " << _q.qy() << " " << _q.qz() << "\t"
	      << _D[0] << " " << _D[1] << " " << _D[2] << "\t"
	      << endl;
}

void BasicMolecule::addTersoffNeighbour(BasicMolecule* m, bool pairType) {
#ifdef TERSOFF_SUPPORT
	for (int j = 0; j < _Tersoff_neighbours_first.size(); j++) {
		if (m->_id == _Tersoff_neighbours_first[j]->id()) {
			this->_Tersoff_neighbours_first[j] = m;
			this->_Tersoff_neighbours_second[j] = pairType;
			return;
		}
	}

	this->_Tersoff_neighbours_first.push_back(m);
	this->_Tersoff_neighbours_second.push_back(pairType);

#else
	warn_tersoff("addTersoffNeighbour");
#endif
}

double BasicMolecule::tersoffParameters(double params[15]) //returns delta_r
{
	const std::vector<Tersoff>* tersoff = &((*components)[_componentid].tersoff());
	const Tersoff* t = &((*(tersoff))[0]);
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
// these are only used when compiling BasicMolecule.cpp and therefore might be inlined without any problems

inline void BasicMolecule::setupCache(const vector<Component>* components) {
	assert(components);
	if (components->size() == 0)
		return;

	this->clearFM();
}

void BasicMolecule::clearFM() {
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
}

void BasicMolecule::calcFM() {
}


/*
 * catches NaN values and missing data
 *
 * @note Use isnan from cmath to check for nan.
 * If that's not available (C99), compare the value with itself. If the value
 * is NaN, the comparison will evaluate to false (according to IEEE754 spec.)
 */
void BasicMolecule::check(unsigned long id) {
	assert(_id == id);
	assert((*components)[_componentid].m() > 0.0);
	for (int d = 0; d < 3; d++) {
		assert(!isnan(_r[d]));
		assert(!isnan(_v[d]));
		assert(!isnan(_D[d]));
		assert(!isnan(_F[d]));
		assert(!isnan(_M[d]));
		assert(!isnan((*components)[_componentid].I(d)));
		assert(!isnan((*components)[_componentid].invI(d)));
	}
}

bool BasicMolecule::isLessThan(const BasicMolecule& m2) const {
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

unsigned long BasicMolecule::totalMemsize() const {
#ifdef TERSOFF_SUPPORT
	int size = sizeof (*this);
	size += _Tersoff_neighbours_first.capacity() * sizeof (BasicMolecule*);
	size += _Tersoff_neighbours_second.capacity() * sizeof (bool);
	return size;
#else
	return sizeof (*this);
#endif
}

std::ostream& operator<<( std::ostream& os, const BasicMolecule& m ) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n" ;
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy() << ", " << m.q().qz()<< "]]\n" ;
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")" ;
	return os;
}

