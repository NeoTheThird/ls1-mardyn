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
#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <vector>
#include <string>
using namespace std;

#include "molecules/Site.h"

/**
Komponente (Molekltyp)

@author Martin Bernreuther
*/
class Component{
public:
  Component(unsigned int id=0);
  //Component(std::istream& istrm);

  void setID(unsigned int id) { m_id=id; };
  unsigned int ID() const { return m_id; }
  unsigned int numSites() const
  {
#ifdef COMPLEX_POTENTIAL_SET
     return numLJcenters() + numCharges() + numQuadrupoles() + numDipoles() + numTersoff();
#else
     return numLJcenters() + numCharges() + numQuadrupoles();
#endif
  }
  unsigned int numLJcenters() const { return m_ljcenters.size(); }
  unsigned numCharges() const { return m_charges.size(); }
#ifdef COMPLEX_POTENTIAL_SET
  unsigned int numDipoles() const { return m_dipoles.size(); }
#endif
  unsigned int numQuadrupoles() const { return m_quadrupoles.size(); }
#ifdef COMPLEX_POTENTIAL_SET
  unsigned int numTersoff() const { return m_tersoff.size(); }
#endif
  double m() const { return m_m; }
  const double* I() const { return m_I; }
  const double* Ipa() const { return m_Ipa; }
  double Ixx() const { return m_I[0]; }
  double Iyy() const { return m_I[1]; }
  double Izz() const { return m_I[2]; }
  double Ixy() const { return m_I[3]; }
  double Ixz() const { return m_I[4]; }
  double Iyz() const { return m_I[5]; }
  double I11() const { return m_Ipa[0]; }
  double I22() const { return m_Ipa[1]; }
  double I33() const { return m_Ipa[2]; }
  void setIxx(double I) { m_I[0]=I; }
  void setIyy(double I) { m_I[1]=I; }
  void setIzz(double I) { m_I[2]=I; }
  void setIxy(double I) { m_I[3]=I; }
  void setIxz(double I) { m_I[4]=I; }
  void setIyz(double I) { m_I[5]=I; }
  void setI11(double I) { m_Ipa[0]=I; }
  void setI22(double I) { m_Ipa[1]=I; }
  void setI33(double I) { m_Ipa[2]=I; }
  void setI(double Ixx=0.,double Iyy=0.,double Izz=0.
           ,double Ixy=0.,double Ixz=0.,double Iyz=0.);
  void setI(unsigned short d,double I) { m_I[d]=I; }
  void addI(double Ixx=0.,double Iyy=0.,double Izz=0.
           ,double Ixy=0.,double Ixz=0.,double Iyz=0.);
  unsigned int rot_dof() const { return m_rot_dof; }

  //const std::vector<Site>& sites() const { return m_sites; }  // see below
  //const Site& site(unsigned int i) const { return m_sites[i]; }
  const std::vector<LJcenter>& ljcenters() const { return m_ljcenters; }
  const LJcenter& ljcenter(unsigned int i) const { return m_ljcenters[i]; }
  const std::vector<Charge>& charges() const { return m_charges; }
  const Charge& charge(unsigned i) const { return m_charges[i]; }
#ifdef COMPLEX_POTENTIAL_SET
  const std::vector<Dipole>& dipoles() const { return m_dipoles; }
  const Dipole& dipole(unsigned int i) const { return m_dipoles[i]; }
#endif
  const std::vector<Quadrupole>& quadrupoles() const { return m_quadrupoles; }
  const Quadrupole& quadrupole(unsigned int i) const { return m_quadrupoles[i]; }
#ifdef COMPLEX_POTENTIAL_SET
  const std::vector<Tersoff>& tersoff() const { return m_tersoff; }
  const Tersoff& tersoff(unsigned int i) const { return m_tersoff[i]; }
#endif
  void incrnumMolecules() { ++m_numMolecules; }
  unsigned long numMolecules() const { return m_numMolecules; }

  void addLJcenter(double x, double y, double z, double m, double eps, double sigma
#ifdef TRUNCATED_SHIFTED
     ,
     double rc, int ownrank
#endif
  );
  void addCharge(double x, double y, double z, double m, double q);
#ifdef COMPLEX_POTENTIAL_SET
  void addDipole( double x, double y, double z,
                  double eMyx, double eMyy, double eMyz, double eMyabs );
#endif
  void addQuadrupole( double x, double y, double z,
                      double eQx, double eQy, double eQz, double eQabs );
#ifdef COMPLEX_POTENTIAL_SET
  void addTersoff( double x, double y, double z,
                   double m, double A, double B, double lambda, double mu, double R,
                   double S, double c, double d, double h, double n, double beta  );
#endif

  /** write information to stream */
  void write(std::ostream& ostrm) const;

  /** write POVray object definition to stream */
  void writePOVobjs(std::ostream& ostrm, std::string para=std::string("pigment {color rgb<1,0,0>}")) const;

#ifdef COMPLEX_POTENTIAL_SET
   double getTersoffRadius() { return this->maximalTersoffExternalRadius; }
   void setTersoffRadius(double mTER) { this->maximalTersoffExternalRadius = mTER; }
#endif

   void writeVIM(std::ostream& ostrm);

private:
  unsigned int m_id;  // IDentification number
  // LJcenter,Dipole,Quadrupole have different size -> not suitable to store in a _Site_-array
  //std::vector<Site> m_sites;
  // use separate vectors instead...
  vector<LJcenter> m_ljcenters;
  vector<Charge> m_charges;
  vector<Quadrupole> m_quadrupoles;
#ifdef COMPLEX_POTENTIAL_SET
  vector<Dipole> m_dipoles;
  vector<Tersoff> m_tersoff;
#endif

  // for performance reasons better(?) omit Site-class indirection and use cached values
  double m_m; // total mass
  //       Ixx,Iyy,Izz,Ixy,Ixz,Iyz
  double m_I[6];  // moments of inertia tensor
  unsigned long m_rot_dof;  // number of rotational degrees of freedom
  double m_Ipa[3];  // moments of inertia for principal axes
  unsigned long m_numMolecules; // number of molecules for this molecule type
#ifdef COMPLEX_POTENTIAL_SET
  double maximalTersoffExternalRadius;
#endif
};

#endif /*COMPONENT_H_*/
