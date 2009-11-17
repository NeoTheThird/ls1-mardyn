/***************************************************************************
 *   Copyright (C) 2009 by Martin Bernreuther and colleagues               *
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
 *                                                                         *
 *   Due to copyleft all future versions of this program must be           *
 *   published as Free Software.                                           *
 ***************************************************************************/

#ifndef POTFORCE_H_
#define POTFORCE_H_

#include "molecules/Comp2Param.h"

/// helper function to calculate the distance between 2 sites
inline void SiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
  dr2=0.;
  for(unsigned short d=0;d<3;++d)
  {
    drs[d]=drm[d]+ds1[d]-ds2[d];
    dr2+=drs[d]*drs[d];
  }
}
inline void minusSiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
  dr2=0.;
  for(unsigned short d=0;d<3;++d)
  {
    drs[d] = ds2[d] - drm[d] - ds1[d];
    dr2+=drs[d]*drs[d];
  }
}

/// calculate potential and force between 2 Lennard-Jones 12-6 centers
//inline void PotForceLJ(const double dr[3], const double& dr2, ParaStrm& params, double f[3], double& u)
inline void PotForceLJ(const double dr[3], const double& dr2
                      , const double& eps24, const double& sig2
                      ,double f[3], double& u6)
{
  //double eps24;
  //params >> eps24;
  //double sig2;
  //params >> sig2;
  double invdr2=1./dr2;
  double lj6=sig2*invdr2; lj6=lj6*lj6*lj6;
  double lj12=lj6*lj6;
  double lj12m6=lj12-lj6;
  u6=eps24*lj12m6;
  double fac=eps24*(lj12+lj12m6)*invdr2;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d];
}

inline void PotForceLJ(const double dr[3]
                      , const double& eps24, const double& sig2
                      , double f[3], double& u6)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForceLJ(dr,dr2,eps24,sig2,f,u6);
}

#ifdef COMPLEX_POTENTIAL_SET
/// calculate potential and force between 2 Dipoles (dr2 given)
inline void PotForce2Dipole(const double dr[3], const double& dr2, const double* eii, const double* ejj
                            , const double& my2, const double& rffac
                            , double f[3], double m1[3], double m2[3], double& u, double& MyRF)
{
  double invdr2=1./dr2;
  double invdr1=sqrt(invdr2);
  double myfac=my2*invdr2*invdr1;
  double costi=0.;
  double costj=0.;
  double cosgij=0.;
  for(unsigned short d=0;d<3;++d)
  {
    const double& drd=dr[d];
    const double& eiid=eii[d];
    const double& ejjd=ejj[d];
    costi+=eiid*drd;
    costj+=ejjd*drd;
    cosgij+=eiid*ejjd;
  }
  costi*=invdr1;
  costj*=invdr1;
  u=myfac*(cosgij-3.*costi*costj);
  MyRF-=rffac*cosgij;
  const double partialRijInvdr1=-3.*u*invdr2;
  const double partialTiInvdr1=-myfac*3.*costj*invdr1;
  const double partialTjInvdr1=-myfac*3.*costi*invdr1;
  const double& partialGij=myfac;
  const double fac=-partialRijInvdr1+(costi*partialTiInvdr1+costj*partialTjInvdr1)*invdr1;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d]-partialTiInvdr1*eii[d]-partialTjInvdr1*ejj[d];

  double eiXej[3], eXrij[3];
  eiXej[0]=eii[1]*ejj[2]-eii[2]*ejj[1];
  eiXej[1]=eii[2]*ejj[0]-eii[0]*ejj[2];
  eiXej[2]=eii[0]*ejj[1]-eii[1]*ejj[0];
  eXrij[0]=eii[1]*dr[2]-eii[2]*dr[1];
  eXrij[1]=eii[2]*dr[0]-eii[0]*dr[2];
  eXrij[2]=eii[0]*dr[1]-eii[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m1[d]=-partialTiInvdr1*eXrij[d]+(-partialGij+rffac)*eiXej[d];
  eXrij[0]=ejj[1]*dr[2]-ejj[2]*dr[1];
  eXrij[1]=ejj[2]*dr[0]-ejj[0]*dr[2];
  eXrij[2]=ejj[0]*dr[1]-ejj[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m2[d]=-partialTjInvdr1*eXrij[d]+(partialGij-rffac)*eiXej[d];
}

/// calculate potential and force between 2 Dipoles
inline void PotForce2Dipole(const double dr[3], const double* eii, const double* ejj
                            , const double& my2, const double& rffac
                            ,double f[3], double m1[3], double m2[3], double& u, double& MyRF)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForce2Dipole(dr,dr2,eii,ejj,my2,rffac,f,m1,m2,u,MyRF);
}
#endif

/// calculate potential and force between 2 Quadrupoles (dr2 given)
inline void PotForce2Quadrupole(const double dr[3], const double& dr2, const double* eii, const double* ejj
                               , const double& q2075
                               , double f[3], double m1[3], double m2[3], double& u)
{
  double invdr2=1./dr2;
  double invdr1=sqrt(invdr2);
  double qfac=q2075*invdr2*invdr2*invdr1;
  double costi=0.;
  double costj=0.;
  double cosgij=0.;
  for(unsigned short d=0;d<3;++d)
  {
    const double& drd=dr[d];
    const double& eiid=eii[d];
    const double& ejjd=ejj[d];
    costi+=eiid*drd;
    costj+=ejjd*drd;
    cosgij+=eiid*ejjd;
  }
  costi*=invdr1;
  costj*=invdr1;
  double cos2ti=costi*costi;
  double cos2tj=costj*costj;
  double term=(cosgij-5.*costi*costj);
  u=qfac*(1.-5.*(cos2ti+cos2tj)-15.*cos2ti*cos2tj+2.*term*term);
  const double partialRijInvdr1=-5.*u*invdr2;
  const double partialTiInvdr1=-qfac*10.*(costi+3.*costi*cos2tj+2.*costj*term)*invdr1;
  const double partialTjInvdr1=-qfac*10.*(costj+3.*cos2ti*costj+2.*costi*term)*invdr1;
  const double& partialGij=qfac*4.*term;
  const double fac=-partialRijInvdr1+(costi*partialTiInvdr1+costj*partialTjInvdr1)*invdr1;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d]-partialTiInvdr1*eii[d]-partialTjInvdr1*ejj[d];

  double eiXej[3], eXrij[3];
  eiXej[0]=eii[1]*ejj[2]-eii[2]*ejj[1];
  eiXej[1]=eii[2]*ejj[0]-eii[0]*ejj[2];
  eiXej[2]=eii[0]*ejj[1]-eii[1]*ejj[0];
  eXrij[0]=eii[1]*dr[2]-eii[2]*dr[1];
  eXrij[1]=eii[2]*dr[0]-eii[0]*dr[2];
  eXrij[2]=eii[0]*dr[1]-eii[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m1[d]=-partialTiInvdr1*eXrij[d]-partialGij*eiXej[d];

  eXrij[0]=ejj[1]*dr[2]-ejj[2]*dr[1];
  eXrij[1]=ejj[2]*dr[0]-ejj[0]*dr[2];
  eXrij[2]=ejj[0]*dr[1]-ejj[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m2[d]=-partialTjInvdr1*eXrij[d]+partialGij*eiXej[d];
}

/// calculate potential and force between 2 Quadrupoles
inline void PotForce2Quadrupole(const double dr[3], const double* eii, const double* ejj
                               , const double& q2075
                               , double f[3], double m1[3], double m2[3], double& u)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForce2Quadrupole(dr,dr2,eii,ejj,q2075,f,m1,m2,u);
}

#ifdef COMPLEX_POTENTIAL_SET
/// calculate potential and force between a Dipole and Quadrupole (dr2 given)
inline void PotForceDiQuadrupole(const double dr[3], const double& dr2, const double* eii, const double* ejj
                                , const double& myq15
                                , double f[3], double m1[3], double m2[3], double& u)
{
  double invdr2=1./dr2;
  double invdr1=sqrt(invdr2);
  double myqfac=myq15*invdr2*invdr2;
  double costi=0.;
  double costj=0.;
  double cosgij=0.;
  for(unsigned short d=0;d<3;++d)
  {
    const double& drd=dr[d];
    const double& eiid=eii[d];
    const double& ejjd=ejj[d];
    costi+=eiid*drd;
    costj+=ejjd*drd;
    cosgij+=eiid*ejjd;
  }
  costi*=invdr1;
  costj*=invdr1;
  double cos2tj=costj*costj;
  u=myqfac*(-costi*(5.*cos2tj-1.)+2.*cosgij*costj);
  const double partialRijInvdr1=-4.*u*invdr2;
  const double partialTiInvdr1=myqfac*(-5.*cos2tj+1.)*invdr1;
  const double partialTjInvdr1=myqfac*2.*(-5.*costi*costj+cosgij)*invdr1;
  const double& partialGij=myqfac*2.*costj;
  const double fac=-partialRijInvdr1+(costi*partialTiInvdr1+costj*partialTjInvdr1)*invdr1;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d]-partialTiInvdr1*eii[d]-partialTjInvdr1*ejj[d];

  double eiXej[3], eXrij[3];
  eiXej[0]=eii[1]*ejj[2]-eii[2]*ejj[1];
  eiXej[1]=eii[2]*ejj[0]-eii[0]*ejj[2];
  eiXej[2]=eii[0]*ejj[1]-eii[1]*ejj[0];
  eXrij[0]=eii[1]*dr[2]-eii[2]*dr[1];
  eXrij[1]=eii[2]*dr[0]-eii[0]*dr[2];
  eXrij[2]=eii[0]*dr[1]-eii[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m1[d]=-partialTiInvdr1*eXrij[d]-partialGij*eiXej[d];

  eXrij[0]=ejj[1]*dr[2]-ejj[2]*dr[1];
  eXrij[1]=ejj[2]*dr[0]-ejj[0]*dr[2];
  eXrij[2]=ejj[0]*dr[1]-ejj[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m2[d]=-partialTjInvdr1*eXrij[d]+partialGij*eiXej[d];
}

/// calculate potential and force between a Dipole and Quadrupole
inline void PotForceDiQuadrupole(const double dr[3], const double* eii, const double* ejj
                                 , const double& myq15
                                 , double f[3], double m1[3], double m2[3], double& u)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForce2Quadrupole(dr,dr2,eii,ejj,myq15,f,m1,m2,u);
}
#endif

/// calculate potential and force between two point charges (dr2 given)
inline void PotForce2Charge( const double dr[3], const double& dr2,
                             const double& q1q2per4pie0, double f[3], double& u )
{
   double invdr2 = 1.0 / dr2;
   double invdr = sqrt(invdr2);
   u = q1q2per4pie0 * invdr;
   const double fac = u*invdr2;
   for(unsigned short d = 0; d < 3; d++) f[d] = fac*dr[d];
}

/// calculate potential and force between an electric charge and a quadrupole
inline void PotForceChargeQuadrupole( const double dr[3], const double& dr2,
                                      const double* ejj, const double& qQ025per4pie0,
                                      double f[3], double m2[3], double& u )
{
   double invdr2 = 1.0 / dr2;
   double invdr = sqrt(invdr2);
   double costj = 0;
   for(unsigned short d = 0; d < 3; d++)
   {
      costj += ejj[d] * dr[d];
   }
   costj *= invdr;
   double qQinv4dr3 = qQ025per4pie0 * invdr*invdr2;
   u = qQinv4dr3 * (3.0*costj*costj - 1);

   const double partialRijInvdr1 = -3.0 * u * invdr2;
   const double partialTjInvdr1 = 6.0 * costj * qQinv4dr3*invdr;
   const double fac = costj*partialTjInvdr1*invdr - partialRijInvdr1;
   for(unsigned short d=0;d<3;++d) f[d] = fac*dr[d] - partialTjInvdr1*ejj[d];

   double minuseXrij[3];
   minuseXrij[0] = ejj[2]*dr[1] - ejj[1]*dr[2];
   minuseXrij[1] = ejj[0]*dr[2] - ejj[2]*dr[0];
   minuseXrij[2] = ejj[1]*dr[0] - ejj[0]*dr[1];
   for(unsigned short d=0;d<3;++d) m2[d] = partialTjInvdr1 * minuseXrij[d];
}

#ifdef COMPLEX_POTENTIAL_SET
/// calculate potential and force between an electric charge and a dipole
inline void PotForceChargeDipole( const double dr[3], const double& dr2,
                                  const double* ejj, const double& minusqmyper4pie0,
                                  double f[3], double m2[3], double& u )
{
   double invdr2 = 1.0 / dr2;
   double invdr = sqrt(invdr2);
   double costj = 0;
   for(unsigned short d = 0; d < 3; d++)
   {
      costj += ejj[d] * dr[d];
   }
   costj *= invdr;
   double uInvcostj = minusqmyper4pie0 * invdr2;
   u = uInvcostj * costj;

   // const double partialRijInvdr1 = -2.0 * u * invdr2;
   const double partialTjInvdr1 = uInvcostj*invdr;
   // const double fac = costj*partialTjInvdr1*invdr - partialRijInvdr1;
   const double fac = 3.0 * u * invdr*invdr2;
   for(unsigned short d=0;d<3;++d) f[d] = fac*dr[d] - partialTjInvdr1*ejj[d];

   double minuseXrij[3];
   minuseXrij[0] = ejj[2]*dr[1] - ejj[1]*dr[2];
   minuseXrij[1] = ejj[0]*dr[2] - ejj[2]*dr[0];
   minuseXrij[2] = ejj[1]*dr[0] - ejj[0]*dr[1];
   for(unsigned short d=0;d<3;++d) m2[d] = partialTjInvdr1 * minuseXrij[d];
}
#endif

/** calculate Potential and Force between molecules (all site-site interactions)
    paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
    (cmp. comp2param.h/.cpp)

   drm == distance FROM j TO i ... !!!
*/
inline void PotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double& Virial )
// ???better calc Virial, when molecule forces are calculated:
//    summing up molecule virials instead of site virials???
{ // Force Calculation
  double f[3];
  double u;
  double drs[3],dr2;  // site distance vector & length^2
  // LJ centers
#ifdef COMPLEX_POTENTIAL_SET
  // no LJ interaction between solid atoms of the same component
  const unsigned int nt1 = mi.numTersoff();
  if((mi.componentid() != mj.componentid()) || !nt1)
  {
#endif
    const unsigned int nc1 = mi.numLJcenters();
    const unsigned int nc2 = mj.numLJcenters();
    for(unsigned int si=0;si<nc1;++si)
    {
      const double* dii=mi.ljcenter_d(si);
      for(unsigned int sj=0;sj<nc2;++sj)
      {
        const double* djj=mj.ljcenter_d(sj);
        SiteSiteDistance(drm,dii,djj,drs,dr2);
        double eps24;
        params >> eps24;
        double sig2;
        params >> sig2;
	cout.flush();
        PotForceLJ(drs,dr2,eps24,sig2,f,u);
#ifdef TRUNCATED_SHIFTED
        double shift6;
        params >> shift6;
        u += shift6;
#endif

// even for interactions within the cell a neighbor might try to add/subtract
// better use atomic...
// and even better use a order where critical sections occure only at some boundary cells...
#ifdef _OPENMP
#pragma omp critical
#endif
        {
          mi.Fljcenteradd(si,f);
          mj.Fljcentersub(sj,f);
        }
        Upot6LJ+=u;
        /*
        u/=6.;
        mi.Upotadd(u);
        mj.Upotadd(u);
        */
        for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
      }
    }
#ifdef COMPLEX_POTENTIAL_SET
  }
#endif

  double m1[3], m2[3];  // angular momenta

  const unsigned ne1 = mi.numCharges();
  const unsigned ne2 = mj.numCharges();
  const unsigned int nq1=mi.numQuadrupoles();
  const unsigned int nq2=mj.numQuadrupoles();
#ifdef COMPLEX_POTENTIAL_SET
  const unsigned int nd1=mi.numDipoles();
  const unsigned int nd2=mj.numDipoles();
#endif
  for(unsigned si = 0; si < ne1; si++)
  {
    const double* dii = mi.charge_d(si);
    // Charge-Charge
    for(unsigned sj = 0; sj < ne2; sj++)
    {
      const double* djj = mj.charge_d(sj);
      double q1q2per4pie0;  // 4pie0 = 1 in reduced units
      params >> q1q2per4pie0;
      SiteSiteDistance(drm, dii, djj, drs, dr2);
      PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);

      mi.Fchargeadd(si, f);
      mj.Fchargesub(sj, f);

      UpotXpoles += u;
      for(unsigned short d = 0; d < 3; d++) Virial += drm[d] * f[d];
    }
    // Charge-Quadrupole
    for(unsigned sj = 0; sj < nq2; sj++)
    {
      const double* djj = mj.quadrupole_d(sj);
      double qQ025per4pie0;  // 4pie0 = 1 in reduced units
      params >> qQ025per4pie0;
      SiteSiteDistance(drm, dii, djj, drs, dr2);
      const double* ejj = mj.quadrupole_e(sj);
      PotForceChargeQuadrupole(drs, dr2, ejj, qQ025per4pie0, f, m2, u);
      
      mi.Fchargeadd(si, f);
      mj.Fquadrupolesub(sj, f);
      mj.Madd(m2);

      UpotXpoles += u;
      for(unsigned short d = 0; d < 3; d++) Virial += drm[d] * f[d];
    }
#ifdef COMPLEX_POTENTIAL_SET
    // Charge-Dipole
    for(unsigned sj = 0; sj < nd2; sj++)
    {
      const double* djj = mj.dipole_d(sj);
      double minusqmyper4pie0;
      params >> minusqmyper4pie0;
      SiteSiteDistance(drm, dii, djj, drs, dr2);
      const double* ejj = mj.dipole_e(sj);
      PotForceChargeDipole(drs, dr2, ejj, minusqmyper4pie0, f, m2, u);
      
      mi.Fchargeadd(si, f);
      mj.Fdipolesub(sj, f);
      mj.Madd(m2);

      UpotXpoles += u;
      for(unsigned short d = 0; d < 3; d++) Virial += drm[d] * f[d];
    }
#endif
  }
  for(unsigned int si=0;si<nq1;++si)
  {
    const double* dii = mi.quadrupole_d(si);
    const double* eii = mi.quadrupole_e(si);

    // Quadrupole-Charge
    for(unsigned sj = 0; sj < ne2; sj++)
    {
      const double* djj = mj.charge_d(sj);
      double qQ025per4pie0;  // 4pie0 = 1 in reduced units
      params >> qQ025per4pie0;
      minusSiteSiteDistance(drm, dii, djj, drs, dr2);
      PotForceChargeQuadrupole(drs, dr2, eii, qQ025per4pie0, f, m1, u);
      
      mi.Fquadrupolesub(si, f);
      mj.Fchargeadd(sj, f);
      mi.Madd(m1);

      UpotXpoles += u;
      for(unsigned short d = 0; d < 3; d++) Virial -= drm[d] * f[d];
    }
    // Quadrupole-Quadrupole -------------------
    for(unsigned int sj=0;sj<nq2;++sj)
    {
      //double drs[3];
      const double* djj=mj.quadrupole_d(sj);
      double q2075;
      params >> q2075;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* ejj=mj.quadrupole_e(sj);
      PotForce2Quadrupole(drs,dr2,eii,ejj,q2075,f,m1,m2,u);

      mi.Fquadrupoleadd(si,f);
      mj.Fquadrupolesub(sj,f);
      mi.Madd(m1);
      mj.Madd(m2);

      UpotXpoles+=u;
      /*
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
#ifdef COMPLEX_POTENTIAL_SET
    // Quadrupole-Dipole -----------------------
    for(unsigned int sj=0;sj<nd2;++sj)
    {
      //double drs[3];
      const double* djj=mj.dipole_d(sj);
      double qmy15;
      params >> qmy15;
      minusSiteSiteDistance(drm, dii, djj, drs, dr2);
      const double* ejj=mj.dipole_e(sj);
      //for(unsigned short d=0;d<3;++d) drs[d]=-drs[d]; // avoid that and toggle add/sub below
      PotForceDiQuadrupole(drs,dr2,ejj,eii,qmy15,f,m2,m1,u);

      mi.Fquadrupolesub(si, f); 
      mj.Fdipoleadd(sj, f);
      mi.Madd(m1);
      mj.Madd(m2);
      UpotXpoles += u;
      for(unsigned short d=0; d<3; d++) Virial -= drm[d]*f[d];
    }
#endif
  }
#ifdef COMPLEX_POTENTIAL_SET
  for(unsigned int si=0;si<nd1;++si)
  {
    const double* dii=mi.dipole_d(si);
    const double* eii=mi.dipole_e(si);
    // Dipole-Charge
    for(unsigned sj = 0; sj < ne2; sj++)
    {
      const double* djj = mj.charge_d(sj);
      double minusqmyper4pie0;
      params >> minusqmyper4pie0;
      minusSiteSiteDistance(drm, dii, djj, drs, dr2);
      PotForceChargeDipole(drs, dr2, eii, minusqmyper4pie0, f, m1, u);
      
      mi.Fdipolesub(si, f);
      mj.Fchargeadd(sj, f);
      mi.Madd(m1);

      UpotXpoles += u;
      for(unsigned short d = 0; d < 3; d++) Virial -= drm[d] * f[d];
    }
    // Dipole-Quadrupole -----------------------
    for(unsigned int sj=0;sj<nq2;++sj)
    {
      //double drs[3];
      const double* djj=mj.quadrupole_d(sj);
      double myq15;
      params >> myq15;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* ejj=mj.quadrupole_e(sj);
      PotForceDiQuadrupole(drs,dr2,eii,ejj,myq15,f,m1,m2,u);

      mi.Fdipoleadd(si,f);
      mj.Fquadrupolesub(sj,f);
      mi.Madd(m1);
      mj.Madd(m2);
      UpotXpoles+=u;
      /*
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
    // Dipole-Dipole ---------------------------
    for(unsigned int sj=0;sj<nd2;++sj)
    {
      //double drs[3];
      const double* djj=mj.dipole_d(sj);
      double my2;
      params >> my2;
      double rffac;
      params >> rffac;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* ejj=mj.dipole_e(sj);
      PotForce2Dipole(drs,dr2,eii,ejj,my2,rffac,f,m1,m2,u,MyRF);

      mi.Fdipoleadd(si,f);
      mj.Fdipolesub(sj,f);
      mi.Madd(m1);
      mj.Madd(m2);
      UpotXpoles+=u;
      /*
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
  }
#endif

  // check whether all parameters were used
  assert(params.eos());
}

#ifdef GRANDCANONICAL
/*
 * calculates the potential energy of the mi-mj interaction,
 * if and only if both of them are fluid (no interaction between Tersoff sites is considered)
 */
inline void FluidPot(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF)
{
  double f[3];
  double u;
  double drs[3], dr2;  // site distance vector & length^2
  // LJ centers
#ifdef COMPLEX_POTENTIAL_SET
  // no LJ interaction between equal solid atoms
  const unsigned int nt1 = mi.numTersoff();
  if((mi.componentid() != mj.componentid()) || !nt1)
  {
#endif
    const unsigned int nc1 = mi.numLJcenters();
    const unsigned int nc2 = mj.numLJcenters();
    for(unsigned int si=0;si<nc1;++si)
    {
      const double* dii=mi.ljcenter_d(si);
      for(unsigned int sj=0;sj<nc2;++sj)
      {
        const double* djj=mj.ljcenter_d(sj);
        SiteSiteDistance(drm,dii,djj,drs,dr2);
        double eps24;
        params >> eps24;
        double sig2;
        params >> sig2;
        PotForceLJ(drs,dr2,eps24,sig2,f,u);
#ifdef TRUNCATED_SHIFTED
        double shift6;
        params >> shift6;
        u += shift6;
#endif
        Upot6LJ+=u;
      }
    }
#ifdef COMPLEX_POTENTIAL_SET
  }
#endif

  double m1[3], m2[3];  // angular momenta

  const unsigned ne1 = mi.numCharges();
  const unsigned ne2 = mj.numCharges();
  const unsigned int nq1=mi.numQuadrupoles();
  const unsigned int nq2=mj.numQuadrupoles();
#ifdef COMPLEX_POTENTIAL_SET
  const unsigned int nd1=mi.numDipoles();
  const unsigned int nd2=mj.numDipoles();
#endif
  for(unsigned si = 0; si < ne1; si++)
  {
    const double* dii = mi.charge_d(si);
    // Charge-Charge
    for(unsigned sj = 0; sj < ne2; sj++)
    {
      const double* djj = mj.charge_d(sj);
      double q1q2per4pie0;  // 4pie0 = 1 in reduced units
      params >> q1q2per4pie0;
      SiteSiteDistance(drm, dii, djj, drs, dr2);
      PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);
      UpotXpoles += u;
    }
    // Charge-Quadrupole
    for(unsigned sj = 0; sj < nq2; sj++)
    {
      const double* djj = mj.quadrupole_d(sj);
      double qQ025per4pie0;  // 4pie0 = 1 in reduced units
      params >> qQ025per4pie0;
      SiteSiteDistance(drm, dii, djj, drs, dr2);
      const double* ejj = mj.quadrupole_e(sj);
      PotForceChargeQuadrupole(drs, dr2, ejj, qQ025per4pie0, f, m2, u);
      UpotXpoles += u;
    }
#ifdef COMPLEX_POTENTIAL_SET
    // Charge-Dipole
    for(unsigned sj = 0; sj < nd2; sj++)
    {
      const double* djj = mj.dipole_d(sj);
      double minusqmyper4pie0;
      params >> minusqmyper4pie0;
      SiteSiteDistance(drm, dii, djj, drs, dr2);
      const double* ejj = mj.dipole_e(sj);
      PotForceChargeDipole(drs, dr2, ejj, minusqmyper4pie0, f, m2, u);
      UpotXpoles += u;
    }
#endif
  }
  for(unsigned int si=0;si<nq1;++si)
  {
    const double* dii = mi.quadrupole_d(si);
    const double* eii = mi.quadrupole_e(si);

    // Quadrupole-Charge
    for(unsigned sj = 0; sj < ne2; sj++)
    {
      const double* djj = mj.charge_d(sj);
      double qQ025per4pie0;  // 4pie0 = 1 in reduced units
      params >> qQ025per4pie0;
      minusSiteSiteDistance(drm, dii, djj, drs, dr2);
      PotForceChargeQuadrupole(drs, dr2, eii, qQ025per4pie0, f, m1, u);
      UpotXpoles += u;
    }
    // Quadrupole-Quadrupole -------------------
    for(unsigned int sj=0;sj<nq2;++sj)
    {
      //double drs[3];
      const double* djj=mj.quadrupole_d(sj);
      double q2075;
      params >> q2075;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* ejj=mj.quadrupole_e(sj);
      PotForce2Quadrupole(drs,dr2,eii,ejj,q2075,f,m1,m2,u);
      UpotXpoles+=u;
    }
#ifdef COMPLEX_POTENTIAL_SET
    // Quadrupole-Dipole -----------------------
    for(unsigned int sj=0;sj<nd2;++sj)
    {
      //double drs[3];
      const double* djj=mj.dipole_d(sj);
      double qmy15;
      params >> qmy15;
      minusSiteSiteDistance(drm, dii, djj, drs, dr2);
      const double* ejj=mj.dipole_e(sj);
      //for(unsigned short d=0;d<3;++d) drs[d]=-drs[d]; // avoid that and toggle add/sub below
      PotForceDiQuadrupole(drs,dr2,ejj,eii,qmy15,f,m2,m1,u);
      UpotXpoles += u;
    }
#endif
  }
#ifdef COMPLEX_POTENTIAL_SET
  for(unsigned int si=0;si<nd1;++si)
  {
    const double* dii=mi.dipole_d(si);
    const double* eii=mi.dipole_e(si);
    // Dipole-Charge
    for(unsigned sj = 0; sj < ne2; sj++)
    {
      const double* djj = mj.charge_d(sj);
      double minusqmyper4pie0;
      params >> minusqmyper4pie0;
      minusSiteSiteDistance(drm, dii, djj, drs, dr2);
      PotForceChargeDipole(drs, dr2, eii, minusqmyper4pie0, f, m1, u);
      UpotXpoles += u;
    }
    // Dipole-Quadrupole -----------------------
    for(unsigned int sj=0;sj<nq2;++sj)
    {
      //double drs[3];
      const double* djj=mj.quadrupole_d(sj);
      double myq15;
      params >> myq15;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* ejj=mj.quadrupole_e(sj);
      PotForceDiQuadrupole(drs,dr2,eii,ejj,myq15,f,m1,m2,u);
      UpotXpoles+=u;
    }
    // Dipole-Dipole ---------------------------
    for(unsigned int sj=0;sj<nd2;++sj)
    {
      //double drs[3];
      const double* djj=mj.dipole_d(sj);
      double my2;
      params >> my2;
      double rffac;
      params >> rffac;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* ejj=mj.dipole_e(sj);
      PotForce2Dipole(drs,dr2,eii,ejj,my2,rffac,f,m1,m2,u,MyRF);
      UpotXpoles+=u;
    }
  }
#endif

  // check whether all parameters were used
  assert(params.eos());
}
#endif

#ifdef COMPLEX_POTENTIAL_SET

//!
//! drij should contain the distance from j to i.
//!
inline double Tersoffbij( Molecule* mi, Molecule* mj,
                          double params[15], double drij[3], double drij1 )
{
   double drik[3];
   double zeta = 0.0;
   unsigned i_curTN = mi->getCurTN();
   Molecule* mk;
   for(unsigned nmk = 0; nmk < i_curTN; nmk++)
   {
      mk = mi->getTersoffNeighbour(nmk);
      if((mk->id() == mj->id()) || (mk->numTersoff() == 0)) continue;
      double drik2 = mk->dist2(*mi, drik);  // distance k->i
      if(params[12] > drik2)
      {
         double drik1 = sqrt(drik2);
         double drdr = 0.0;
         for(int d=0; d<3; d++) drdr += drij[d] * drik[d];
         double h_min_costheta = params[2] - drdr/(drij1*drik1);
         double g = params[11] - params[3]/(params[4] + h_min_costheta*h_min_costheta);
         if(drik1 > params[0])
         {
            double kC = 0.5 + 0.5*cos((drik1 - params[0]) * params[10]);
            zeta += kC * g;
         }
         else zeta += g;
#ifndef NDEBUG
         /*
         cout << "\t\t\t\t\t\tI(" << mi->r(0) << " / " << mi->r(1) << " / " << mi->r(2) << ")\n";
         cout << "\t\t\t\t\t\tJ(" << mj->r(0) << " / " << mj->r(1) << " / " << mj->r(2) << ")\n";
         cout << "\t\t\t\t\t\tK(" << mk->r(0) << " / " << mk->r(1) << " / " << mk->r(2) << ")\n";
         cout << "\t\t\t\t\t\ti-j = (" << drij[0] << "; " << drij[1] << "; " << drij[2] << ")\n";
         cout << "\t\t\t\t\t\ti-k = (" << drik[0] << "; " << drik[1] << "; " << drik[2] << ")\n";
         cout << "\t\t\t\t\tangle " << mi->id() << "-" << mj->id() << "-" << mk->id() << ": arccos " << drdr/(drij1*drik1) << "\n";
         cout << "\t\t\t\t\tg = " << g << " = " << params[11] << " - " << params[3] << "/(" << params[4] << " + " << h_min_costheta << "^2)\n";
         */
#endif
      }
   }
   double BZ = params[8]*zeta;
   double BZtoN = pow(BZ, params[9]);
   double b = pow(1.0 + BZtoN, params[14]);
#ifndef NDEBUG
   /*
   cout << "\t\t\t\tpow(" << params[8] << "*zeta, " << params[9] << ") = " << BZtoN << "\n";
   cout << "\t\t\t\tb(" << mi->id() << ", " << mj->id() << ") = " << b << " (zeta = " << zeta << ", -0.5/n = " << params[14] << ")\n";
   */
#endif
   return b;
}

//! @brief returns the sum, whereas Uij is added to the double reference
//!
//! drij should contain the distance from j to i.
//!
inline double TersoffUIJplusUJI( Molecule* mi, Molecule* mj, double params[15],
                                 double drij[3], double drij2, double& UpotTersoff )
{
   double Uij = 0.0;
   double Uji = 0.0;
   double drji[3];
   double drij1 = sqrt(drij2);
   for(int d=0; d<3; d++) drji[d] = -drij[d];  // drji: distance i->j

   double UA = params[13] * exp(params[7] * drij1);
   double UR = params[ 5] * exp(params[6] * drij1);
   double bij = Tersoffbij(mi, mj, params, drij, drij1);
   double bji = Tersoffbij(mj, mi, params, drji, drij1);
   if(drij1 > params[0])
   {
      double kChalf = 0.25 + 0.25*cos((drij1 - params[0]) * params[10]);
      Uij = kChalf * (UR + bij*UA);
      Uji = kChalf * (UR + bji*UA);
   }
   else
   {
      Uij = 0.5*(UR + bij*UA);
      Uji = 0.5*(UR + bji*UA);
   }
   UpotTersoff += Uij;
   return Uij+Uji;
}
//!
//! drij should contain the distance from j to i.
//!
inline double TersoffUIJattr( Molecule* mi, Molecule* mj,
                              double params[15], double drij[3], double drij2 )
{
   double drij1 = sqrt(drij2);
   double UA = params[13] * exp(params[7] * drij1);
   double bij = Tersoffbij(mi, mj, params, drij, drij1);
   if(drij1 > params[0])
   {
      double kChalf = 0.25 + 0.25*cos((drij1 - params[0]) * params[10]);
      return kChalf*bij*UA;
   }
   else return 0.5*bij*UA;
}

//! @brief calculate Tersoff potential for a single atom
//!
//! parameters: R, S, h, c^2, d^2, A, -lambda, -mu, beta, n_i, pi/(S-R), 1+(c/d)^2, S^2, -B, -0.5/n_i
//! A "Molecule" may have at most a single Tersoff centre.
//!
//! the function itself returns the sum of all Uij, Uji, and the attractive
//! part of Ujk, i.e. all potential energy terms that are influenced by atom i. 
//! However, only the sum over the Uij is added to UpotTersoff.
//!
//! used for computing the Tersoff potential based force on the atom
//!
inline double TersoffPotential(Molecule* mi, double params[15], double& UpotTersoff)
{
   double Ui = 0.0;

   double distanceVector[3];
   unsigned i_curTN = mi->getCurTN();
   Molecule *mj, *mk;
   for(unsigned nmj = 0; nmj < i_curTN; nmj++)
   {
      mj = mi->getTersoffNeighbour(nmj);
      if(mj->numTersoff() == 0) continue;
#ifndef NDEBUG
      // cout << "\t\tconsidering " << mi->id() << "-" << mj->id() << " interaction.\n";
#endif
      double drij2 = mj->dist2(*mi, distanceVector);  // distance j->i
      if(params[12] > drij2)
         Ui += TersoffUIJplusUJI(mi, mj, params, distanceVector, drij2, UpotTersoff);

      unsigned j_curTN = mj->getCurTN();
      for(unsigned nmk = 0; nmk < j_curTN; nmk++)
      {
         mk = mj->getTersoffNeighbour(nmk);
         if(mk->id() == mi->id() || (mk->numTersoff() == 0)) continue; 
#ifndef NDEBUG
      // cout << "\t\tconsidering " << mj->id() << "-" << mk->id() << " attraction.\n";
#endif
         double drjk2 = mk->dist2(*mj, distanceVector);  // distance k->j
         if(params[12] > drjk2)
            Ui += TersoffUIJattr(mj, mk, params, distanceVector, drjk2);
      }
   }

   return Ui;
}
inline void TersoffPotForce(Molecule* mi, double params[15], double& UpotTersoff, double delta_r)
{
   double f[3];
   if(mi->numTersoff() == 0) return;
   double offsetU = TersoffPotential(mi, params, UpotTersoff);

#ifndef NDEBUG
   // cout << "atom " << mi->id() << " offset U = " << offsetU << "\n";
#endif

   for(int d=0; d<3; d++)
   {
      mi->move(d, delta_r);

      double irrelevantU = 0.0;
      double currentU = TersoffPotential(mi, params, irrelevantU);
#ifndef NDEBUG
      // cout << "\tdim. " << d << " current U = " << currentU << "\n";
#endif
      
      f[d] = (offsetU - currentU) / delta_r;
      mi->move(d, -delta_r);
   }

   mi->Ftersoffadd(0, f);
}

/*
//! @brief calculate Tersoff Potential and Force
//!
//! this potential is not symmetric, it compares i-j and i-k bond angles.
//! a separate call is necessary for the second half of the potential, where
//! j-i and j-k bond angles are compared.
//!
//! @brief calculate Tersoff Potential and Force
//!
//! this potential is not symmetric, it compares i-j and i-k bond angles.
//! a separate call is necessary for the second half of the potential, where
//! j-i and j-k bond angles are compared.
//!
inline void TersoffPotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& UpotTersoff, double& Virial)
{
   double f[3];
   double drs[3],dr2;  // site distance vector & length^2
   // LJ centers
   const unsigned int nt1=mi.numTersoff();
   const unsigned int nt2=mj.numTersoff();

   double R, S, h, cc, dd, A, minus_B, minus_lambda, minus_mu, beta, ni;
   double fcij, dfcij, fcik, dfcik, one_plus_ccinvdd;
   double S_minus_R_inv, xi, xiik;

   double drmik[3], drsik[3];
   double drik2, rik, h_minus_costheta, zeta, g_theta, contribution, vpr, costheta;

   double rij, uattr, urep, frep, norm, minus_two_phi_a;

   double tmp_j2, tmp_k2, tmp_jk, tmp_2,
          tmp_3, tmp_4, tmp_5, tmp_6, tmp_grad, b_ij;  // Bezeichnungen aus IMD uebernommen
   double dzeta_i[3], dzeta_j[3], dcos_j[3], dcos_k[3], force_j[3], fk[3];
   map<int, double*> dzeta_k[3];

   double tmp_virial = 0.0;  // wie in IMD, hier also nur Attraktion

   unsigned i_curTN = mi.getCurTN();
   Molecule* mk_first;
   // map<Molecule*, bool>* iTn = mi.getTersoffNeighbours();
   // map<Molecule*, bool>::iterator mk;

#ifndef NDEBUG
   cout.precision(2);
   cout << "Tersoff interaction between\n\t" << mi.id() << " (" << mi.r(0)
        << " / " << mi.r(1) << " / " << mi.r(2) << ")\nand \t" << mj.id() << " (" << mj.r(0)
        << " / " << mj.r(1) << " / " << mj.r(2) << ").\n";
#endif

   for(unsigned int si=0; si < nt1; si++)
   {
      const double* dii = mi.tersoff_d(si);
      for(unsigned int sj=0; sj < nt2; sj++)
      {
         params >> R >> S >> h >> cc >> dd
                >> A >> minus_B >> minus_lambda >> minus_mu >> beta >> ni;

         S_minus_R_inv = 1.0 / (S - R);
         one_plus_ccinvdd = 1.0 + cc/dd;

         const double* djj=mj.tersoff_d(sj);
         SiteSiteDistance(drm,dii,djj,drs,dr2);
#ifndef NDEBUG
         // cout << "Distance ij: (" << -drs[0] << " / " << -drs[1] << " / " << -drs[2] << ")\n";
#endif
         rij = sqrt(dr2);

         // cutoff function: take out of loop for optimization
         if (rij < R)
         {
            fcij  = 1.0; 
            dfcij = 0.0;
         }
         else if (rij > S)
         {
            continue;
         }
         else
         {
            xi = (S - rij) * S_minus_R_inv;
            fcij = xi*xi * (3.0 - xi - xi);
            dfcij = 6.0*xi*(xi - 1.0) * S_minus_R_inv;

#ifndef NDEBUG
            cout << "fC  =  " << fcij << "\n";  //XIII//
#endif
         }

         // In IMD, the corresponding section in imd_forces_covalent.c
         // is only responsible for the attractive part of the interaction,
         // hence we first calculate the repuslive part separately
         //
         urep = A*exp(minus_lambda*rij);  // 2U_rep (wird nur zur Haelfte angerechnet)

         // je hoeher (positiv) dieser Wert ist, desto staerker die Abstossung:
         frep = -0.5 * (dfcij*urep + fcij*minus_lambda*urep);  // -dU/dr

         UpotTersoff += 0.5 * fcij * urep;
         norm = 1.0 / rij;

         // Achtung: drs ist in j-i-Richtung
         for(unsigned short d = 0; d < 3; d++) f[d] = norm*frep*drs[d];
         mi.Ftersoffadd(si,f);  // i in j-i-Richtung (Abstossung)
         mj.Ftersoffsub(sj,f);  // deshalb bei j in j-i Richtung abziehen
         for(unsigned short d = 0; d < 3; d++) Virial += drm[d]*f[d];

         zeta = 0.0;
         for(int d=0; d < 3; d++)
         {
            dzeta_i[d] = 0.0;
            dzeta_j[d] = 0.0;
            dzeta_k[d] = map<int, double*>();
         }

         // for(mk = iTn->begin(); mk != iTn->end(); mk++)
         for(unsigned mk = 0; mk < i_curTN; mk++)
         {
            mk_first = mi.getTersoffNeighbour(mk);
            if(mk_first->id() == mj.id()) continue;  // never evaluate i-j angle to i-j
            // if(mk->first->id() == mj.id()) continue;  // never evaluate i-j angle to i-j etc.
#ifndef NDEBUG
            cout << "Evaluating " << mi.id() << "-" << mj.id()
                 << "-" << mk_first->id() << " contribution. ";  //XIII//
#endif
            mk_first->dist2(mi, drmik);
            const unsigned int nt3 = mk_first->numTersoff();
            for(int d=0; d < 3; d++)
               dzeta_k[d][mk_first->id()] = new double[nt3];
            for(unsigned int sk=0; sk < nt3; sk++)
            {
               const double* dkk = mk_first->tersoff_d(sk);
               SiteSiteDistance(drmik,dii,dkk,drsik,drik2);  // drsik: k->i
               rik = sqrt(drik2);

               // cutoff function: take out of loop for optimization
               if (rik < R)
               {
	         fcik  = 1.0; 
	         dfcik = 0.0;
               }
               else if (rik > S)
               {
                  continue;
               }
               else
               {
                  xiik = (S - rik) * S_minus_R_inv;
                  fcik = xiik*xiik * (3.0 - xiik - xiik);
                  dfcik = 6.0*xiik*(xiik - 1.0) * S_minus_R_inv;
               }

               // angular term
               tmp_jk = 1.0 / (rij*rik);  
               vpr = drs[0]*drsik[0] + drs[1]*drsik[1] + drs[2]*drsik[2];
               costheta = vpr * tmp_jk;
               h_minus_costheta = h - costheta;  // in IND: tmp_1
               tmp_2 = 1.0 / (dd + h_minus_costheta*h_minus_costheta);
               g_theta = one_plus_ccinvdd - cc*tmp_2;
               contribution = fcik * g_theta;
               zeta  += contribution; 
#ifndef NDEBUG
               cout << "fC = " << fcik  //XIII//
                    << ", contrib = " << contribution << "\n";  //XIII//
#endif

               tmp_j2 = costheta / dr2;
               tmp_k2 = costheta / drik2;

               // derivatives of cos(theta), 
               // note that dcos_i + dcos_j + dcos_k = 0 
               for(int d=0; d < 3; d++)
               {
                  dcos_j[d] = tmp_j2*drs[d] - tmp_jk*drsik[d];
                  dcos_k[d] = tmp_k2*drsik[d] - tmp_jk*drs[d];
               }

               tmp_3 = 2.0 * cc * h_minus_costheta * tmp_2 * tmp_2 * fcik;
               tmp_grad = dfcik / rik * g_theta;

               // derivatives of zeta; dzeta_i is not the full derivative 
               for(int d=0; d < 3; d++)
               {
                  dzeta_k[d][mk_first->id()][sk] = - tmp_grad*drsik[d] - tmp_3*dcos_k[d];
                  dzeta_i[d] += tmp_3*dcos_k[d] + tmp_grad*drsik[d];
                  dzeta_j[d] -= tmp_3 * dcos_j[d];
               }
            }
         }

         minus_two_phi_a = minus_B * exp(minus_mu * rij);
         tmp_4 = pow(beta*zeta, ni);
         b_ij = pow(1.0+tmp_4, -0.5/ni);
         uattr = b_ij * minus_two_phi_a;
         UpotTersoff += 0.5 * fcij * uattr;
#ifndef NDEBUG
         cout << "zeta = " << zeta << ", bij = " << b_ij << ", fC = " << fcij  //XIII//
              << ", Eij_attr = " << 0.5 * fcij * uattr << "\n";  //XIII//
#endif

         tmp_6 = 0.5 * (minus_two_phi_a*fcij*minus_mu*b_ij + dfcij*b_ij*minus_two_phi_a) / rij;
         if ( zeta == 0.0 )   // only one neighbor of i
            tmp_5 = 0.0;
         else
            tmp_5 = 0.25 * fcij * uattr * tmp_4 / (zeta * (1.0+tmp_4));

         // tmp force on particle j
         for(int d=0; d < 3; d++)
         {
            force_j[d] = tmp_5*dzeta_j[d] + tmp_6*drs[d];
         }

         // update force on particle k
         // for(mk = iTn->begin(); mk != iTn->end(); mk++)
         for(unsigned mk = 0; mk < i_curTN; mk++)
         {
            mk_first = mi.getTersoffNeighbour(mk);
            if(mk_first->id() == mj.id()) continue;
            // if(mk->first->id() == mj.id()) continue; etc.

            const unsigned int nt3 = mk_first->numTersoff();
            for(unsigned int sk=0; sk < nt3; sk++)
            {
               const double* dkk = mk_first->tersoff_d(sk);
               SiteSiteDistance(drmik,dii,dkk,drsik,drik2);
               rik = sqrt(drik2);
               if(rik > S) continue;

               for(int d=0; d<3; d++)
               {
                  fk[d] = tmp_5 * dzeta_k[d][mk_first->id()][sk];
                  tmp_virial -= fk[d] * drsik[d];
               }
               mk_first->Ftersoffadd(sk, fk);
            }

            for(int d=0; d<3; d++) 
            {
               delete dzeta_k[d][mk_first->id()];
               dzeta_k[d][mk_first->id()] = (double*)0;
            }
         }

         // update force on particle j
#ifndef NDEBUG
         cout << "attractive force on " << mj.id() << ": (" << force_j[0] << " / "  //XIII
              << force_j[1] << " / " << force_j[2] << ")\n";  //XIII
#endif
         mj.Ftersoffadd(sj, force_j);
         for(int d=0; d<3; d++) tmp_virial -= drs[d]*force_j[d];

         // update force on particle i 
         for(int d=0; d<3; d++) f[d] = tmp_5*dzeta_i[d] - force_j[d];
         mi.Ftersoffadd(si, f);

#ifndef NDEBUG
         cout << "\n";  //XIII//
#endif
      }
   }

   Virial += tmp_virial;

   // check whether all parameters were accessed
   assert(params.eos());
}
*/

/*
inline void TersoffPotForce_original_cutoff(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& UpotTersoff, double& Virial)
{
   double f[3];
   double drs[3],dr2;  // site distance vector & length^2
   // LJ centers
   const unsigned int nt1=mi.numTersoff();
   const unsigned int nt2=mj.numTersoff();

   double R, S, h, cc, dd, A, minus_B, minus_lambda, minus_mu, beta, ni;
   double fcij, dfcij, fcik, dfcik, one_plus_ccinvdd;
   //APPROX double S_minus_R_inv, xi, xiik;
   double cut_tmp, cut_tmp_j, cut_tmp_k;  

   double drmik[3], drsik[3];
   double drik2, rik, h_minus_costheta, zeta, g_theta, contribution, vpr, costheta;

   double rij, uattr, urep, frep, norm, minus_two_phi_a;

   double tmp_j2, tmp_k2, tmp_jk, tmp_2,
          tmp_3, tmp_4, tmp_5, tmp_6, tmp_grad, b_ij;  // Bezeichnungen aus IMD uebernommen
   double dzeta_i[3], dzeta_j[3], dcos_j[3], dcos_k[3], force_j[3], fk[3];
   map<int, double*> dzeta_k[3];

   double tmp_virial = 0.0;  // wie in IMD, hier also nur Attraktion

   map<Molecule*, bool>* iTn = mi.getTersoffNeighbours();
   map<Molecule*, bool>::iterator mk;

#ifndef NDEBUG
   cout.precision(2);
   cout << "Tersoff interaction between\n\t" << mi.id() << " (" << mi.r(0)
        << " / " << mi.r(1) << " / " << mi.r(2) << ")\nand \t" << mj.id() << " (" << mj.r(0)
        << " / " << mj.r(1) << " / " << mj.r(2) << ").\n";
#endif

   for(unsigned int si=0; si < nt1; si++)
   {
      const double* dii = mi.tersoff_d(si);
      for(unsigned int sj=0; sj < nt2; sj++)
      {
         params >> R >> S >> h >> cc >> dd
                >> A >> minus_B >> minus_lambda >> minus_mu >> beta >> ni;

         cut_tmp = 3.1415927 / (S - R);
         one_plus_ccinvdd = 1.0 + cc/dd;

         const double* djj=mj.tersoff_d(sj);
         SiteSiteDistance(drm,dii,djj,drs,dr2);
#ifndef NDEBUG
         // cout << "Distance ij: (" << -drs[0] << " / " << -drs[1] << " / " << -drs[2] << ")\n";
#endif
         rij = sqrt(dr2);

         // cutoff function: take out of loop for optimization
         if (rij < R)
         {
            fcij  = 1.0; 
            dfcij = 0.0;
         }
         else if (rij > S)
         {
            continue;
         }
         else
         {
            cut_tmp_j = cut_tmp * (rij - R);
            fcij   = 0.5 * ( 1.0 + cos( cut_tmp_j ) );
            dfcij  = -0.5 * cut_tmp * sin( cut_tmp_j );

#ifndef NDEBUG
            cout << "fC  =  " << fcij << "\n";  //XIII//
#endif
         }

         // In IMD, the corresponding section in imd_forces_covalent.c
         // is only responsible for the attractive part of the interaction,
         // hence we first calculate the repuslive part separately
         //
         urep = A*exp(minus_lambda*rij);  // 2U_rep (wird nur zur Haelfte angerechnet)

         // je hoeher (positiv) dieser Wert ist, desto staerker die Abstossung:
         frep = -0.5 * (dfcij*urep + fcij*minus_lambda*urep);  // -dU/dr

         UpotTersoff += 0.5 * fcij * urep;
         norm = 1.0 / rij;

         // Achtung: drs ist in j-i-Richtung
         for(unsigned short d = 0; d < 3; d++) f[d] = norm*frep*drs[d];
         mi.Ftersoffadd(si,f);  // i in j-i-Richtung (Abstossung)
         mj.Ftersoffsub(sj,f);  // deshalb bei j in j-i Richtung abziehen
         for(unsigned short d = 0; d < 3; d++) Virial += drm[d]*f[d];

         zeta = 0.0;
         for(int d=0; d < 3; d++)
         {
            dzeta_i[d] = 0.0;
            dzeta_j[d] = 0.0;
            dzeta_k[d] = map<int, double*>();
         }
         for(mk = iTn->begin(); mk != iTn->end(); mk++)
         {
            if(mk->first->id() == mj.id()) continue;  // never evaluate i-j angle to i-j
#ifndef NDEBUG
            cout << "Evaluating " << mi.id() << "-" << mj.id()
                 << "-" << mk->first->id() << " contribution. ";  //XIII//
#endif
            mk->first->dist2(mi, drmik);
            const unsigned int nt3 = mk->first->numTersoff();
            for(int d=0; d < 3; d++)
               dzeta_k[d][mk->first->id()] = new double[nt3];
            for(unsigned int sk=0; sk < nt3; sk++)
            {
               const double* dkk = mk->first->tersoff_d(sk);
               SiteSiteDistance(drmik,dii,dkk,drsik,drik2);  // drsik: k->i
               rik = sqrt(drik2);

               // cutoff function: take out of loop for optimization
               if (rik < R)
               {
	         fcik  = 1.0; 
	         dfcik = 0.0;
               }
               else if (rik > S)
               {
                  continue;
               }
               else
               {
                  cut_tmp_k = cut_tmp * (rik - R);
                  fcik   = 0.5 * ( 1.0 + cos( cut_tmp_k ) );
                  dfcik  = -0.5 * cut_tmp * sin( cut_tmp_k );
               }

               // angular term
               tmp_jk = 1.0 / (rij*rik);  
               vpr = drs[0]*drsik[0] + drs[1]*drsik[1] + drs[2]*drsik[2];
               costheta = vpr * tmp_jk;
               h_minus_costheta = h - costheta;  // in IND: tmp_1
               tmp_2 = 1.0 / (dd + h_minus_costheta*h_minus_costheta);
               g_theta = one_plus_ccinvdd - cc*tmp_2;
               contribution = fcik * g_theta;
               zeta  += contribution; 
#ifndef NDEBUG
               cout << "fC = " << fcik  //XIII//
                    << ", contrib = " << contribution << "\n";  //XIII//
#endif

               tmp_j2 = costheta / dr2;
               tmp_k2 = costheta / drik2;

               // derivatives of cos(theta), 
               // note that dcos_i + dcos_j + dcos_k = 0 
               for(int d=0; d < 3; d++)
               {
                  dcos_j[d] = tmp_j2*drs[d] - tmp_jk*drsik[d];
                  dcos_k[d] = tmp_k2*drsik[d] - tmp_jk*drs[d];
               }

               tmp_3 = 2.0 * cc * h_minus_costheta * tmp_2 * tmp_2 * fcik;
               tmp_grad = dfcik / rik * g_theta;

               // derivatives of zeta; dzeta_i is not the full derivative 
               for(int d=0; d < 3; d++)
               {
                  dzeta_k[d][mk->first->id()][sk] = - tmp_grad*drsik[d] - tmp_3*dcos_k[d];
                  dzeta_i[d] += tmp_3*dcos_k[d] + tmp_grad*drsik[d];
                  dzeta_j[d] -= tmp_3 * dcos_j[d];
               }
            }
         }

         minus_two_phi_a = minus_B * exp(minus_mu * rij);
         tmp_4 = pow(beta*zeta, ni);
         b_ij = pow(1.0+tmp_4, -0.5/ni);
         uattr = b_ij * minus_two_phi_a;
         UpotTersoff += 0.5 * fcij * uattr;
#ifndef NDEBUG
         cout << "zeta = " << zeta << ", bij = " << b_ij << ", fC = " << fcij  //XIII//
              << ", Eij_attr = " << 0.5 * fcij * uattr << "\n";  //XIII//
#endif

         tmp_6 = 0.5 * (minus_two_phi_a*fcij*minus_mu*b_ij + dfcij*b_ij*minus_two_phi_a) / rij;
         if ( zeta == 0.0 )   // only one neighbor of i
            tmp_5 = 0.0;
         else
            tmp_5 = 0.25 * fcij * uattr * tmp_4 / (zeta * (1.0+tmp_4));

         // tmp force on particle j
         for(int d=0; d < 3; d++)
         {
            force_j[d] = tmp_5*dzeta_j[d] + tmp_6*drs[d];
         }

         // update force on particle k
         for(mk = iTn->begin(); mk != iTn->end(); mk++)
         {
            if(mk->first->id() == mj.id()) continue;

            const unsigned int nt3 = mk->first->numTersoff();
            for(unsigned int sk=0; sk < nt3; sk++)
            {
               const double* dkk = mk->first->tersoff_d(sk);
               SiteSiteDistance(drmik,dii,dkk,drsik,drik2);
               rik = sqrt(drik2);
               if(rik > S) continue;

               for(int d=0; d<3; d++)
               {
                  fk[d] = tmp_5 * dzeta_k[d][mk->first->id()][sk];
                  tmp_virial -= fk[d] * drsik[d];
               }
               mk->first->Ftersoffadd(sk, fk);
            }

            for(int d=0; d<3; d++) 
            {
               delete dzeta_k[d][mk->first->id()];
               dzeta_k[d][mk->first->id()] = (double*)0;
            }
         }

         // update force on particle j
#ifndef NDEBUG
         cout << "attractive force on " << mj.id() << ": (" << force_j[0] << " / "  //XIII
              << force_j[1] << " / " << force_j[2] << ")\n";  //XIII
#endif
         mj.Ftersoffadd(sj, force_j);
         for(int d=0; d<3; d++) tmp_virial -= drs[d]*force_j[d];

         // update force on particle i 
         for(int d=0; d<3; d++) f[d] = tmp_5*dzeta_i[d] - force_j[d];
         mi.Ftersoffadd(si, f);

#ifndef NDEBUG
         cout << "\n";  //XIII//
#endif
      }
   }

   Virial += tmp_virial;

   // check whether all parameters were accessed
   assert(params.eos());
}
*/

/*
inline void ancient_mh_implemented_TersoffPotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& UpotTersoff, double& Virial)
{
   double f[3];
   double u;
   double drs[3],dr2;  // site distance vector & length^2
   // LJ centers
   const unsigned int nt1=mi.numTersoff();
   const unsigned int nt2=mj.numTersoff();

   double A, minus_B, RR, minus_lambda, minus_mu, S, S_minus_R_inv,
             cc, dd, one_plus_ccinvdd, h, ni, minus_two_ni_inv, beta;
   double fChalf, minusdfChalfdr, zeta, xi, b, fCik, xiik;

   double drmik[3], drsik[3];
   double drik2, rik, h_minus_costheta, contribution, vpr, costheta;

   double rij, urep, uattr, full_potential, norm, fac;

   map<Molecule*, bool>* iTn;
   map<Molecule*, bool>::iterator mk;

#ifndef NDEBUG
   cout << "\ni: " << mi.id() << "\tj: " << mj.id() << "\tTersoff: " << nt1 << "x" << nt2 << "\n";  //XIII//
#endif

   for(unsigned int si=0; si < nt1; si++)
   {
      const double* dii = mi.tersoff_d(si);
      for(unsigned int sj=0; sj < nt2; sj++)
      {
         const double* djj=mj.tersoff_d(sj);
         SiteSiteDistance(drm,dii,djj,drs,dr2);

         params >> A;
         params >> minus_B;
         params >> minus_lambda;
         params >> minus_mu;
         params >> RR;
         params >> S;
         params >> S_minus_R_inv;
         params >> cc;
         params >> dd;
         params >> one_plus_ccinvdd;
         params >> h;
         params >> ni;
         params >> minus_two_ni_inv;
         params >> beta;

         if(dr2 > S*S)
         {
            continue;
         }
         if(RR > dr2)
         {
            fChalf = 0.5;
            minusdfChalfdr = 0.0;
         }
         else
         {
            xi = (S - sqrt(dr2)) * S_minus_R_inv;
            fChalf = xi*xi * (1.5 - xi);
#ifndef NDEBUG
            cout << "fC  =  2 * " << fChalf << "\n";  //XIII//
#endif
            minusdfChalfdr = 3.0 * xi * (1.0 - xi) * S_minus_R_inv;
         }
         rij = sqrt(dr2);

         iTn = mi.getTersoffNeighbours();
         zeta = 0;
         cout << "\n" << mi.id() << " (chk " << mj.id() << ") has " << iTn->size() << " neighbours with Tersoff sites.\n";  //XIII//
         for(mk = iTn->begin(); mk != iTn->end(); mk++)
         {
            if(mk->first->id() == mj.id()) continue;  // never evaluate i-j angle to i-j
#ifndef NDEBUG
            cout << "processing " << mi.id() << "/" << mj.id()
                 << " Tersoff interaction with " << mk->first->id() << ".\n";  //XIII//
#endif
            mk->first->dist2(mi, drmik);
            const unsigned int nt3 = mk->first->numTersoff();
            for(unsigned int sk=0; sk < nt3; sk++)
            {
               const double* dkk = mk->first->tersoff_d(sk);
               SiteSiteDistance(drmik,dii,dkk,drsik,drik2);
               rik = sqrt(drik2);
#ifndef NDEBUG
               cout << "distance = " << rik << "\n";  //XIII//
#endif
               if(rik > S) continue;  // approximation
               if(RR > drik2) fCik = 1.0;
               else
               {
                  xiik = (S - rik) * S_minus_R_inv;
                  fCik = xiik*xiik * (3.0 - 2.0*xiik);
               }

               vpr = drs[0]*drsik[0] + drs[1]*drsik[1] + drs[2]*drsik[2];
               costheta = vpr / (rij*rik);
               h_minus_costheta = h - costheta;
               contribution = fCik * (one_plus_ccinvdd - cc/(dd + h_minus_costheta*h_minus_costheta));
#ifndef NDEBUG
               cout << "h = " << h << ",  rij = " << rij << ", rik = " << sqrt(drik2)  //XIII//
                    << ",  vector product = " << vpr << ",  costheta = " << costheta << "\n";  //XIII//
               cout << "h-cos = " << h_minus_costheta << ", fCik = " << fCik << ", contribution = g = "  //XIII//
                    << contribution << "\n";  //XIII//
#endif
               zeta += contribution;
            }
         }
         b = pow(1.0 + pow(beta*zeta, ni), minus_two_ni_inv);

         urep = A*exp(minus_lambda*rij);
         uattr = minus_B*b*exp(minus_mu*rij);
         full_potential = urep + uattr;
         u = fChalf * full_potential;

         norm = 1.0 / rij;
         fac = norm * (minusdfChalfdr*full_potential - fChalf*(minus_mu*uattr + minus_lambda*urep));
         for(unsigned short d = 0; d < 3; d++) f[d] = fac*drs[d];
	
#ifdef _OPENMP
#pragma omp critical
#endif
         {
            mi.Ftersoffadd(si,f);
            mj.Ftersoffsub(sj,f);
         }
         UpotTersoff+=u;
         for(unsigned short d = 0; d < 3; d++) Virial += drm[d]*f[d];

#ifndef NDEBUG
         cout << "\nzeta: " << zeta << "\tb_ij: " << b << "\n";  //XIII//
         cout << "u_rep: " << urep << "\tu_attr: " << uattr << "\tfC: " << 2.0*fChalf << "\n";  //XIII//
         cout << "f_rep: " << (minusdfChalfdr*urep - fChalf*minus_lambda*urep)  //XIII//
              << "\tf_attr: " << (minusdfChalfdr*uattr - fChalf*minus_lambda*uattr) << "\n";  //XIII//
         cout << "u = Vij/2 = " << u << "\tf -> i: " << f[0] << " " << f[1] << " " << f[2] << "\n";  //XIII//
#endif
      }
   }

   // check whether all parameters were accessed
   assert(params.eos());
}
*/
#endif

#endif /*POTFORCE_H_*/
