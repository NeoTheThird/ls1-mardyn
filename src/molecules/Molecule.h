/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther et al.                       *
 *   bernreuther@hlrs.de                                                   *
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
#ifndef MOLECULE_H_
#define MOLECULE_H_

#ifdef COMPLEX_POTENTIAL_SET
#define MAXTN 12
#endif

/**
Molecule modeled as LJ and/or Tersoff centre with point polarities

@author Martin Bernreuther and Martin Horsch
*/

#include "molecules/Quaternion.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"
#include "molecules/Component.h"
#include "integrators/Integrator.h"
class Domain;

#include <set>
#include <vector>
#include <map>
#include <iostream>
#include <string>
using namespace std;

#include <cassert>

class Molecule{
//friend class integrators::Integrator;
//friend class integrators::Leapfrog;
public:
  enum streamtype { RESTART };

  Molecule(unsigned long id=0, int componentid=0
          , double rx=0., double ry=0., double rz=0.
          , double vx=0., double vy=0., double vz=0.
          , double q0=0., double q1=0., double q2=0., double q3=0.
          , double Dx=0., double Dy=0., double Dz=0.
          , const std::vector<Component>* components=NULL
          );
  Molecule(const Molecule& m);
  Molecule(std::istream& istrm, streamtype type, const std::vector<Component>* components=NULL);

  ~Molecule()
  {
     /*
#ifdef COMPLEX_POTENTIAL_SET
     if(m_tersoff->size() > 0)
     {
        map<Molecule*, bool>::iterator itn;
        for(itn=m_Tersoff_neighbours.begin(); itn!=m_Tersoff_neighbours.end(); itn++)
        {
           itn->first->m_Tersoff_neighbours.erase(this);
#ifndef NDEBUG
           // cout << "deleting Tersoff connection " << itn->first->m_id << "->" << m_id << ".\n";
	   // cout.flush();
#endif
        }
     }
#endif  
     */

     delete[] m_sites_d; delete[] m_osites_e; delete[] m_sites_F;
  }
  /** get the ID */
  unsigned long id() const { return m_id; }
#ifdef GRANDCANONICAL
  void setid(unsigned long id) { this->m_id = id; }
#endif
  /** get the Component */
  int componentid() const { return m_componentid; }
  /** get the position */
  double r(unsigned short d) const { return m_r[d]; }
  /** get the velocity */
  double v(unsigned short d) const { return m_v[d]; }
  /** get the Orientation */
  const Quaternion& q() const { return m_q; }

  inline void move(int d, double dr) { m_r[d] += dr; }

  /** get the rotatational speed */
  double D(unsigned short d) const { return m_D[d]; }

  /** get F */
  double F(unsigned short d) const {return m_F[d]; }
  /** get M */
  double M(unsigned short d) const {return m_M[d]; }

  //double Upot() const { return m_Upot; }
  double Utrans() const { return .5*m_m*(m_v[0]*m_v[0]+m_v[1]*m_v[1]+m_v[2]*m_v[2]); }
  double Urot();
  double m() { return m_m; }

  /** get number of sites */
  unsigned int numSites() const { return m_numsites; }
  unsigned int numLJcenters() const { return m_ljcenters->size(); }
  unsigned int numCharges() const { return m_charges->size(); }
  unsigned int numQuadrupoles() const { return m_quadrupoles->size(); }
#ifdef COMPLEX_POTENTIAL_SET
  unsigned int numDipoles() const { return m_dipoles->size(); }
  unsigned int numTersoff() const { assert(m_tersoff); return m_tersoff->size(); }
#endif

  const double* site_d(unsigned int i) const { return &(m_sites_d[3*i]); }
  const double* osite_e(unsigned int i) const { return &(m_osites_e[3*i]); }
  const double* site_F(unsigned int i) const { return &(m_sites_F[3*i]); }
  const double* ljcenter_d(unsigned int i) const { return &(m_ljcenters_d[3*i]); }
  const double* ljcenter_F(unsigned int i) const { return &(m_ljcenters_F[3*i]); }
  const double* charge_d(unsigned int i) const { return &(m_charges_d[3*i]); }
  const double* charge_F(unsigned int i) const { return &(m_charges_F[3*i]); }
  const double* quadrupole_d(unsigned int i) const { return &(m_quadrupoles_d[3*i]); }
  const double* quadrupole_e(unsigned int i) const { return &(m_quadrupoles_e[3*i]); }
  const double* quadrupole_F(unsigned int i) const { return &(m_quadrupoles_F[3*i]); }
#ifdef COMPLEX_POTENTIAL_SET
  const double* dipole_d(unsigned int i) const { return &(m_dipoles_d[3*i]); }
  const double* dipole_e(unsigned int i) const { return &(m_dipoles_e[3*i]); }
  const double* dipole_F(unsigned int i) const { return &(m_dipoles_F[3*i]); }
  const double* tersoff_d(unsigned int i) const { return &(m_tersoff_d[3*i]); }
  const double* tersoff_F(unsigned int i) const { return &(m_tersoff_F[3*i]); }
#endif

  /** get object memory size */
  static unsigned long memsize() { return sizeof(Molecule); }

  /** set the position */
  void setr(unsigned short d, double r) { m_r[d]=r; }

  /** calculate the difference vector and return the square (euclidean) distance */
  double dist2(const Molecule& a, double dr[]) const
    { double d2=0.; for(unsigned short d=0;d<3;++d) { dr[d]=a.m_r[d]-m_r[d]; d2+=dr[d]*dr[d]; } return d2; }
  double dist2(const Molecule& a, double L[3], double dr[]) const;
  /** calculate and return the square velocity */
  double v2() const {return m_v[0]*m_v[0]+m_v[1]*m_v[1]+m_v[2]*m_v[2]; }

  void setFM(double Fx, double Fy, double Fz, double Mx, double My, double Mz)
    { m_F[0]=Fx; m_F[1]=Fy; m_F[2]=Fz; m_M[0]=Mx; m_M[1]=My; m_M[2]=Mz; }
  void scale_v(double s) { for(unsigned short d=0;d<3;++d) m_v[d]*=s; }
#ifdef COMPLEX_POTENTIAL_SET
  void scale_v(double s, double offx, double offy, double offz);
#endif
  void scale_F(double s) { for(unsigned short d=0;d<3;++d) m_F[d]*=s; }
  void scale_D(double s) { for(unsigned short d=0;d<3;++d) m_D[d]*=s; }

  void Fadd(const double a[]) { for(unsigned short d=0;d<3;++d) m_F[d]+=a[d]; }
  void Fsub(const double a[]) { for(unsigned short d=0;d<3;++d) m_F[d]-=a[d]; }

  void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) m_M[d]+=a[d]; }
  void Msub(const double a[]) { for(unsigned short d=0;d<3;++d) m_M[d]-=a[d]; }

  inline void move(const double* dr)
  {
     m_r[0] += dr[0];
     m_r[1] += dr[1];
     m_r[2] += dr[2];
  }
  /*
  void move(const double* dr, const double* globalLength)
  {
     for(unsigned d=0; d < 3; d++)
     {
        m_r[d] += dr[d];
        if(m_r[d] > globalLength[d]) m_r[d] -= globalLength[d];
        else if(m_r[d] < 0.0) m_r[d] += globalLength[d];
     }
  }
  */
  void vsub(const double ax, const double ay, const double az)
  {
     m_v[0] -= ax; m_v[1] -= ay; m_v[2] -= az;
  }
#ifdef COMPLEX_POTENTIAL_SET
  void vadd(const double ax, const double ay, const double az)
  {
     m_v[0] += ax; m_v[1] += ay; m_v[2] += az;
  }

  void setXY() { fixedx = m_r[0]; fixedy = m_r[1]; }
  void setXYZ() { fixedx = m_r[0]; fixedy = m_r[1]; fixedz = m_r[2]; }

  void resetXY()
  {
     m_v[0] = 0.0;
     m_v[1] = 0.0;
     m_F[0] = 0.0; 
     m_F[1] = 0.0;
     m_r[0] = fixedx;
     m_r[1] = fixedy;
  }
  void resetXYZ()
  {
     m_v[0] = 0.0;
     m_v[1] = 0.0;
     m_v[2] = 0.0;
     m_F[0] = 0.0; 
     m_F[1] = 0.0;
     m_F[2] = 0.0;
     m_r[0] = fixedx;
     m_r[1] = fixedy;
     m_r[2] = fixedz;
  }
#endif

  void Fsiteadd(unsigned int i, double a[])
    { double* Fsite=&(m_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fsitesub(unsigned int i, double a[])
    { double* Fsite=&(m_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fljcenteradd(unsigned int i, double a[])
    { double* Fsite=&(m_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fljcentersub(unsigned int i, double a[])
    { double* Fsite=&(m_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fchargeadd(unsigned int i, double a[])
    { double* Fsite=&(m_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fchargesub(unsigned int i, double a[])
    { double* Fsite=&(m_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fquadrupoleadd(unsigned int i, double a[])
    { double* Fsite=&(m_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fquadrupolesub(unsigned int i, double a[])
    { double* Fsite=&(m_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
#ifdef COMPLEX_POTENTIAL_SET
  void Fdipoleadd(unsigned int i, double a[])
    { double* Fsite=&(m_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fdipolesub(unsigned int i, double a[])
    { double* Fsite=&(m_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Ftersoffadd(unsigned int i, double a[])
    { double* Fsite=&(m_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Ftersoffsub(unsigned int i, double a[])
    { double* Fsite=&(m_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
#endif

  //void Upotadd(double u) { m_Upot+=u; }
  //void Upotsub(double u) { m_Upot-=u; }

  void upd_preF(double dt, double vcorr=1., double Dcorr=1.);
  void upd_cache();
  void upd_postF(double dt_halve, double& summv2, double& sumIw2);

  //! calculate summv2 and sumIw2
  void calculate_mv2_Iw2(double& summv2, double& sumIw2);
#ifdef COMPLEX_POTENTIAL_SET
  void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz);
#endif

  /** write information to stream */
  void write(std::ostream& ostrm) const;
  /** write binary information to stream */
  void save_restart(std::ostream& ostrm) const;

  /*
   * veraltet, aber aus Versehen von Martin Buchholz im Code gelassen (M.H. 17. Februar 2009)
   * nach dem neuen Schema von Martin Buchholz ist eine Cell als STL-Liste organisiert,
   * deshalb ist das hier unnötig.
   *
  //  Linked cells
  Molecule* nextinCell() const { return m_nextincell; }
   */
  
  static void setDomain(Domain* domain);
  static void setblubb(double a);

#ifdef COMPLEX_POTENTIAL_SET
  inline unsigned getCurTN() { return this->m_curTN; }
  inline Molecule* getTersoffNeighbour(unsigned i) { return this->m_Tersoff_neighbours_first[i]; }
  inline bool getPairCode(unsigned i) { return this->m_Tersoff_neighbours_second[i]; }
  // map<Molecule*, bool>* getTersoffNeighbours() { return &(this->m_Tersoff_neighbours); }
  inline void clearTersoffNeighbourList() { this->m_curTN = 0; }
  // void clearTersoffNeighbourList() { this->m_Tersoff_neighbours.clear(); }
  void addTersoffNeighbour(Molecule* m, bool pairType);
  double tersoffParameters(double params[15]); //returns delta_r
#endif

  // clear forces and moments
  void clearFM();
  void check(unsigned long id);

private:
  
  static Domain* _domain;
  unsigned long m_id; // IDentification number of that molecule
  int m_componentid;  // IDentification number of its component type
  double m_r[3];  // position coordinates
  double m_v[3];  // velocity
  Quaternion m_q; // orientation
  double m_D[3];  // angular momentum

  double m_F[3];  // forces
  double m_M[3];  // moments
  //double m_Upot;  // potential energy

#ifdef COMPLEX_POTENTIAL_SET
  const vector<Dipole>* m_dipoles;
  const vector<Tersoff>* m_tersoff;
#endif
  const vector<LJcenter>* m_ljcenters;
  const vector<Charge>* m_charges;
  const vector<Quadrupole>* m_quadrupoles;

  double m_m; // total mass
  double m_I[3],m_invI[3];  // moment of inertia for principal axes and it's inverse
  std::size_t m_numsites; // number of sites
  std::size_t m_numorientedsites; // number of oriented sites (subset of sites)
  // global site coordinates relative to site origin
  // row order: dx1,dy1,dz1,dx2,dy2,dz2,...
  double* m_sites_d;
  double *m_ljcenters_d, *m_charges_d, *m_quadrupoles_d;
#ifdef COMPLEX_POTENTIAL_SET
  double *m_dipoles_d, *m_tersoff_d;
#endif
  // site orientation
  double* m_osites_e;
  double* m_quadrupoles_e;
#ifdef COMPLEX_POTENTIAL_SET
  double* m_dipoles_e;
#endif
  // site Forces
  // row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
  double* m_sites_F;
  double *m_ljcenters_F, *m_charges_F, *m_quadrupoles_F;
#ifdef COMPLEX_POTENTIAL_SET
  double *m_dipoles_F, *m_tersoff_F;
#endif
 
  /*
   * veraltet, aber aus Versehen von Martin Buchholz im Code gelassen (M.H. 17. Februar 2009)
   *
  // used by CELL data structure for intrusive single-linked "collision" lists
  Molecule* m_nextincell;
  // used by Cells::addMolecule to modify m_nextincell
  void setNextinCell(Molecule* next) { m_nextincell=next; }
  //  friend void Cells::addMolecule(Molecule* atom);
   *
   */

#ifdef COMPLEX_POTENTIAL_SET
  //! near molecules with Tersoff centres
  //! key: pointer to a molecule, value: whether interaction counts for U_pot
  Molecule* m_Tersoff_neighbours_first[MAXTN];
  bool m_Tersoff_neighbours_second[MAXTN];
  int m_curTN;
  // map<Molecule*, bool> m_Tersoff_neighbours;
  double fixedx, fixedy, fixedz;
#endif

  // setup cache values/properties
  void setupCache(const std::vector<Component>* components);
  // calculate forces and moments for already given site forces
  void calcFM();
};

#endif /*MOLECULE_H_*/
