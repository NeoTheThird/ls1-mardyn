/***************************************************************************
 *   Copyright (C) 2005 by Martin Bernreuther   *
 *   Martin.Bernreuther@informatik.uni-stuttgart.de   *
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
 ***************************************************************************/
#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <cmath>

/**
 @author Martin Bernreuther, Wolfgang Eckhardt

 @param fp_type must be either float or double
 */
template <typename fp_type>
class Quaternion {

public:
	Quaternion(fp_type qw = 0., fp_type qx = 0., fp_type qy = 0., fp_type qz = 0.)
			: m_qw(qw), m_qx(qx), m_qy(qy), m_qz(qz) {
	}

	fp_type qw() const { return m_qw; }
	fp_type qx() const { return m_qx; }
	fp_type qy() const { return m_qy; }
	fp_type qz() const { return m_qz; }
	fp_type magnitude2() {
		return m_qw * m_qw + m_qx * m_qx + m_qy * m_qy + m_qz * m_qz;
	}

	void scale(fp_type s) {
		m_qw *= s;
		m_qx *= s;
		m_qy *= s;
		m_qz *= s;
	}

	void scaleinv(fp_type s) {
		m_qw /= s;
		m_qx /= s;
		m_qy /= s;
		m_qz /= s;
	}

	void normalize() {
		scaleinv(sqrt(magnitude2()));
	}

	void conjugate() {
		m_qx = -m_qx;
		m_qy = -m_qy;
		m_qz = -m_qz;
	}

	void inverse() {
		conjugate();
		scaleinv(magnitude2());
	}

	void add(const Quaternion<fp_type>& q) {
		m_qw += q.m_qw;
		m_qx += q.m_qx;
		m_qy += q.m_qy;
		m_qz += q.m_qz;
	}

	void operator *=(const Quaternion<fp_type>& q);
	/**
	 * apply the rotation represented by tis quaternion to d
	 * @param d the vector to be rotated
	 * @param result vector of the rotation
	 */
	void rotate(const fp_type d[3], fp_type drot[3]) const;

	/**
	 * apply the rotation represented by tis quaternion to d. The result vector
	 * is stored to d.
	 * @param d the vector to be rotated
	 */
	void rotate(fp_type d[3]) const {
		fp_type drot[3];
		rotate(d, drot);
		d[0] = drot[0];
		d[1] = drot[1];
		d[2] = drot[2];
	}

	void rotateinv(const double d[3], double drot[3]) const;
	void rotateinv(double d[3]) const {
		double drot[3];
		rotateinv(d, drot);
		d[0] = drot[0];
		d[1] = drot[1];
		d[2] = drot[2];
	}

	//void differentiate_dbl(const double w[3], Quaternion& dqdt) const;
	void differentiate(const fp_type w[3], Quaternion<fp_type>& dqdt) const;
	//  { differentiate_dbl(w,dqdt); dqdt.scale(.5); }
	void getRotMatrix(fp_type R[3][3]) const;
	void getRotinvMatrix(double R[3][3]) const;

private:
	fp_type m_qw, m_qx, m_qy, m_qz; // components

};

template <typename fp_type>
void Quaternion<fp_type>::operator *=(const Quaternion<fp_type>& q) {
	fp_type qw=m_qw*q.m_qw-m_qx*q.m_qx-m_qy*q.m_qy-m_qz*q.m_qz;
	fp_type qx=m_qw*q.m_qx+m_qx*q.m_qw+m_qy*q.m_qz-m_qz*q.m_qy;
	fp_type qy=m_qw*q.m_qy+m_qy*q.m_qw+m_qz*q.m_qx-m_qx*q.m_qz;
	m_qz=m_qw*q.m_qz+m_qz*q.m_qw+m_qx*q.m_qy-m_qy*q.m_qx;
	m_qy=qy;
	m_qx=qx;
	m_qw=qw;
}

template <typename fp_type>
void Quaternion<fp_type>::rotate(const fp_type d[3], fp_type drot[3]) const {
	fp_type ww=m_qw*m_qw;
	fp_type xx=m_qx*m_qx;
	fp_type yy=m_qy*m_qy;
	fp_type zz=m_qz*m_qz;
	fp_type xy=m_qx*m_qy;
	fp_type zw=m_qz*m_qw;
	fp_type xz=m_qx*m_qz;
	fp_type yw=m_qy*m_qw;
	//       1-2*(yy+zz)
	drot[0]=(ww+xx-yy-zz)*d[0]+2.*(xy+zw)*d[1]+2.*(xz-yw)*d[2];
	fp_type yz=m_qy*m_qz;
	fp_type xw=m_qx*m_qw;
	//                       1-2*(xx+zz)
	drot[1]=2.*(xy-zw)*d[0]+(ww-xx+yy-zz)*d[1]+2.*(yz+xw)*d[2];
	//                                       1-2*(xx+yy)
	drot[2]=2.*(xz+yw)*d[0]+2.*(yz-xw)*d[1]+(ww-xx-yy+zz)*d[2];
}

template <typename fp_type>
void Quaternion<fp_type>::rotateinv(const double d[3], double drot[3]) const {
	double ww=m_qw*m_qw;
	double xx=m_qx*m_qx;
	double yy=m_qy*m_qy;
	double zz=m_qz*m_qz;
	double xy=m_qx*m_qy;
	double zw=m_qz*m_qw;
	double xz=m_qx*m_qz;
	double yw=m_qy*m_qw;
	//       1-2*(yy+zz)
	drot[0]=(ww+xx-yy-zz)*d[0]+2.*(xy-zw)*d[1]+2.*(xz+yw)*d[2];
	double yz=m_qy*m_qz;
	double xw=m_qx*m_qw;
	//                       1-2*(xx+zz)
	drot[1]=2.*(xy+zw)*d[0]+(ww-xx+yy-zz)*d[1]+2.*(yz-xw)*d[2];
	//                                       1-2*(xx+yy)
	drot[2]=2.*(xz-yw)*d[0]+2.*(yz+xw)*d[1]+(ww-xx-yy+zz)*d[2];
}


template <typename fp_type>
void Quaternion<fp_type>::differentiate(const fp_type w[3], Quaternion<fp_type>& dqdt) const {
	dqdt.m_qw=.5*(-m_qx*w[0]-m_qy*w[1]-m_qz*w[2]);
	dqdt.m_qx=.5*( m_qw*w[0]-m_qz*w[1]+m_qy*w[2]);
	dqdt.m_qy=.5*( m_qz*w[0]+m_qw*w[1]-m_qx*w[2]);
	dqdt.m_qz=.5*(-m_qy*w[0]+m_qx*w[1]+m_qw*w[2]);
}

template <typename fp_type>
void Quaternion<fp_type>::getRotMatrix(fp_type R[3][3]) const {
	fp_type ww=m_qw*m_qw;
	fp_type xx=m_qx*m_qx;
	fp_type yy=m_qy*m_qy;
	fp_type zz=m_qz*m_qz;
	fp_type xy=m_qx*m_qy;
	fp_type zw=m_qz*m_qw;
	fp_type xz=m_qx*m_qz;
	fp_type yw=m_qy*m_qw;
	R[0][0]=ww+xx-yy-zz;
	R[0][1]=2.*(xy+zw);
	R[0][2]=2.*(xz-yw);
	fp_type yz=m_qy*m_qz;
	fp_type xw=m_qx*m_qw;
	R[1][0]=2.*(xy-zw);
	R[1][1]=ww-xx+yy-zz;
	R[1][2]=2.*(yz+xw);
	R[2][0]=2.*(xz+yw);
	R[2][1]=2.*(yz-xw);
	R[2][2]=ww-xx-yy+zz;
}

template <typename fp_type>
void Quaternion<fp_type>::getRotinvMatrix(double R[3][3]) const {
	fp_type ww=m_qw*m_qw;
	fp_type xx=m_qx*m_qx;
	fp_type yy=m_qy*m_qy;
	fp_type zz=m_qz*m_qz;
	fp_type xy=m_qx*m_qy;
	fp_type zw=m_qz*m_qw;
	fp_type xz=m_qx*m_qz;
	fp_type yw=m_qy*m_qw;
	R[0][0]=ww+xx-yy-zz;
	R[0][1]=2.*(xy-zw);
	R[0][2]=2.*(xz+yw);
	fp_type yz=m_qy*m_qz;
	fp_type xw=m_qx*m_qw;
	R[1][0]=2.*(xy+zw);
	R[1][1]=ww-xx+yy-zz;
	R[1][2]=2.*(yz-xw);
	R[2][0]=2.*(xz-yw);
	R[2][1]=2.*(yz+xw);
	R[2][2]=ww-xx-yy+zz;
}


#endif /*QUATERNION_H_*/
