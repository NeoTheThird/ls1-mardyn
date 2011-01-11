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

#include <cmath>

#include "molecules/Comp2Param.h"

using namespace std;

void Comp2Param::initialize(const vector<Component>& components, const vector<double>& mixcoeff, double epsRF, double rc, double rcLJ, bool wallLJ)
{
  m_numcomp=components.size();
  m_ssparatbl.redim(m_numcomp,m_numcomp);

  for(unsigned int compi=0;compi<m_numcomp;++compi)
  {
    for(unsigned int compj=0;compj<m_numcomp;++compj)
    {
       ParaStrm& pstrmij=m_ssparatbl(compi,compj);
       pstrmij = ParaStrm();
    }
  }  

  // interaction between LJ centers
  vector<double>::const_iterator mixpos=mixcoeff.begin();
  for(unsigned int compi=0;compi<m_numcomp;++compi)
  {
    ParaStrm& pstrmii=m_ssparatbl(compi,compi);
#ifndef NDEBUG
    cout.flush();
    cout << "\nLJ parameters for " << compi << " and " << compi << ": ";
    cout.flush();
#endif
    unsigned int nci=components[compi].numLJcenters();
    double epsi,sigi,epsj,sigj,epsilon24,sigma2;
#ifdef TRUNCATED_SHIFTED
    double shift6i;
#endif
#ifdef COMPLEX_POTENTIAL_SET
    unsigned nti = components[compi].numTersoff();
    // no LJ interaction between solid atoms belonging to the same component
    if(!nti)
    {
#endif
      // interaction between equal particles
      for(unsigned int centeri=0;centeri<nci;++centeri)
      {
        const LJcenter& ljcenteri=static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
        epsi=ljcenteri.eps();
        sigi=ljcenteri.sigma();
#ifdef TRUNCATED_SHIFTED
        shift6i = ljcenteri.shift6();
#endif
        for(unsigned int centerj=0;centerj<nci;++centerj)
        {
          const LJcenter& ljcenterj=static_cast<const LJcenter&>(components[compi].ljcenter(centerj));
          epsj=ljcenterj.eps();
          sigj=ljcenterj.sigma();
          epsilon24=24.*sqrt(epsi*epsj);
          sigma2=.5*(sigi+sigj); sigma2*=sigma2;
          pstrmii << epsilon24;
#ifndef NDEBUG
          cout << "eps24=" << epsilon24 << " ";
#endif
          pstrmii << sigma2;
#ifndef NDEBUG
          cout << "sig2=" << sigma2 << " ";
          cout.flush();
#endif
#ifdef TRUNCATED_SHIFTED
          pstrmii << shift6i;
#ifndef NDEBUG
          cout << "shift6=" << shift6i << " ";
#endif
#endif
        }
      }
#ifdef COMPLEX_POTENTIAL_SET
    }
#endif  
    // unlike LJ interaction
    for(unsigned int compj=compi+1;compj<m_numcomp;++compj)
    {
#ifdef COMPLEX_POTENTIAL_SET
       unsigned ntj = components[compj].numTersoff();
       if(wallLJ || !nti || !ntj)
       {
#endif
#ifndef NDEBUG
         cout.flush();
         cout << "\nLJ parameters for " << compi << " and " << compj << ": ";
         cout.flush();
#endif
         ParaStrm& pstrmij=m_ssparatbl(compi,compj);
         unsigned int ncj=components[compj].numLJcenters();
         double xi=*mixpos; ++mixpos;
         double eta=*mixpos; ++mixpos;
#ifdef TRUNCATED_SHIFTED
         double shift6combined, sigperrc2, sigperrc6;
#endif
         for(unsigned int centeri=0;centeri<nci;++centeri)
         {
           const LJcenter& ljcenteri=static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
           epsi=ljcenteri.eps();
           sigi=ljcenteri.sigma();
           for(unsigned int centerj=0;centerj<ncj;++centerj)
           {
             const LJcenter& ljcenterj=static_cast<const LJcenter&>(components[compj].ljcenter(centerj));
             epsj=ljcenterj.eps();
             sigj=ljcenterj.sigma();
             epsilon24=24.*xi*sqrt(epsi*epsj);
             sigma2=eta*.5*(sigi+sigj); sigma2*=sigma2;
#ifdef TRUNCATED_SHIFTED
             sigperrc2 = sigma2/(rcLJ*rcLJ);
             sigperrc6 = sigperrc2*sigperrc2*sigperrc2;
             shift6combined = epsilon24 * (sigperrc6 - sigperrc6*sigperrc6);
#endif
             pstrmij << epsilon24;
#ifndef NDEBUG
             cout << "eps24=" << epsilon24 << " ";
#endif
             pstrmij << sigma2;
#ifndef NDEBUG
             cout << "sig2=" << sigma2 << " ";
             cout.flush();
#endif
#ifdef TRUNCATED_SHIFTED
             pstrmij << shift6combined;  
#ifndef NDEBUG
             cout << "shift6=" << shift6combined << " ";
#endif
#endif
           }
         }
         ParaStrm& pstrmji=m_ssparatbl(compj,compi);
#ifndef NDEBUG
         cout.flush();
         cout << "\nLJ parameters for " << compj << " and " << compi << ": ";
         cout.flush();
#endif
         for(unsigned int centerj=0;centerj<ncj;++centerj)
         {
           const LJcenter& ljcenterj=static_cast<const LJcenter&>(components[compj].ljcenter(centerj));
           epsj=ljcenterj.eps();
           sigj=ljcenterj.sigma();
           for(unsigned int centeri=0;centeri<nci;++centeri)
           {
             const LJcenter& ljcenteri=static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
             epsi=ljcenteri.eps();
             sigi=ljcenteri.sigma();
             epsilon24=24.*xi*sqrt(epsi*epsj);
             sigma2=eta*.5*(sigi+sigj); sigma2*=sigma2;
#ifdef TRUNCATED_SHIFTED
             sigperrc2 = sigma2/(rcLJ*rcLJ);
             sigperrc6 = sigperrc2*sigperrc2*sigperrc2;
             shift6combined = epsilon24 * (sigperrc6 - sigperrc6*sigperrc6);
#endif
             pstrmji << epsilon24;
#ifndef NDEBUG
             cout << "eps24=" << epsilon24 << " ";
#endif
             pstrmji << sigma2;
#ifndef NDEBUG
             cout << "sig2=" << sigma2 << " ";
             cout.flush();
#endif
#ifdef TRUNCATED_SHIFTED
             pstrmji << shift6combined;  
#ifndef NDEBUG
             cout << "shift6=" << shift6combined << " ";
#endif
#endif
           }
         }
#ifdef COMPLEX_POTENTIAL_SET
       }
#endif
    }
  }

  for(unsigned int compi=0;compi<m_numcomp;++compi)
  {
    unsigned nei = components[compi].numCharges();
    unsigned int nqi=components[compi].numQuadrupoles();
#ifdef COMPLEX_POTENTIAL_SET
    unsigned int ndi=components[compi].numDipoles();
    unsigned nti = components[compi].numTersoff();
#endif
    for(unsigned int compj=0;compj<m_numcomp;++compj)
    {
      ParaStrm& pstrmij=m_ssparatbl(compi,compj);
#ifndef NDEBUG
         cout.flush();
         cout << "\nElectrostatics parameters for " << compi << " and " << compj << ": ";
         cout.flush();
#endif
      unsigned nej = components[compj].numCharges();
      unsigned int nqj=components[compj].numQuadrupoles();
#ifdef COMPLEX_POTENTIAL_SET
      unsigned int ndj=components[compj].numDipoles();
      unsigned ntj = components[compj].numTersoff();

      /*
       * no electrostatics between solid atoms
       */
      if(nti && ntj) continue;
#endif
      for(unsigned ei = 0; ei < nei; ei++)
      {
        const Charge& chargei = static_cast<const Charge&>(components[compi].charge(ei));
        double cvali = chargei.q();
        // Charge-Charge
        for(unsigned ej = 0; ej < nej; ej++)
        {
          const Charge& chargej = static_cast<const Charge&>(components[compj].charge(ej));
          double cvalj = chargej.q();
          double q1q2per4pie0 = cvali*cvalj;  // 4pi*epsilon0 = Einheit!!!
          pstrmij << q1q2per4pie0;
#ifndef NDEBUG
          cout << "11=" << q1q2per4pie0 << " ";
          cout.flush();
#endif
        }
        // Charge-Quadrupole
        for(unsigned qj=0; qj < nqj; qj++)
        {
          const Quadrupole& quadrupolej = static_cast<const Quadrupole&>(components[compj].quadrupole(qj));
          double absqj=quadrupolej.absQ();
          double qQ025per4pie0 = 0.25*cvali*absqj;
          pstrmij << qQ025per4pie0;
#ifndef NDEBUG
          cout << "14=" << qQ025per4pie0 << " ";
#endif
        }
#ifdef COMPLEX_POTENTIAL_SET
        // Charge-Dipole
        for(unsigned int dj = 0; dj < ndj; dj++)
        {
          const Dipole& dipolej = static_cast<const Dipole&>(components[compj].dipole(dj));
          double absmyj = dipolej.absMy();
          double minusqmyper4pie0 = -cvali*absmyj;
          pstrmij << minusqmyper4pie0;
#ifndef NDEBUG
          cout << "12=" << minusqmyper4pie0 << " ";
#endif
        }
#endif
      }
      for(unsigned int qi=0;qi<nqi;++qi)
      {
        const Quadrupole& quadrupolei=static_cast<const Quadrupole&>(components[compi].quadrupole(qi));
        double absqi=quadrupolei.absQ();
        // Quadrupole-Charge
        for(unsigned ej = 0; ej < nej; ej++)
        {
          const Charge& chargej = static_cast<const Charge&>(components[compj].charge(ej));
          double cvalj = chargej.q();
          double qQ025per4pie0 = 0.25*cvalj*absqi;
          pstrmij << qQ025per4pie0;
#ifndef NDEBUG
          cout << "41=" << qQ025per4pie0 << " ";
#endif
        }
        // Quadrupole-Quadrupole
        for(unsigned int qj=0;qj<nqj;++qj)
        {
          const Quadrupole& quadrupolej=static_cast<const Quadrupole&>(components[compj].quadrupole(qj));
          double absqj=quadrupolej.absQ();
          double q2075=.75*absqi*absqj;
          pstrmij << q2075;
#ifndef NDEBUG
          cout << "44=" << q2075 << " ";
#endif
        }
#ifdef COMPLEX_POTENTIAL_SET
        // Quadrupole-Dipole
        for(unsigned int dj=0;dj<ndj;++dj)
        {
          const Dipole& dipolej=static_cast<const Dipole&>(components[compj].dipole(dj));
          double absmyj=dipolej.absMy();
          double qmy15=1.5*absqi*absmyj;
          pstrmij << qmy15;
#ifndef NDEBUG
          cout << "42=" << qmy15 << " ";
#endif
        }
#endif
      }
#ifdef COMPLEX_POTENTIAL_SET
      for(unsigned int di=0;di<ndi;++di)
      {
        const Dipole& dipolei=static_cast<const Dipole&>(components[compi].dipole(di));
        double absmyi=dipolei.absMy();
        double epsRFInvrc3=2.*(epsRF-1.)/((rc*rc*rc)*(2.*epsRF+1.));
        // Dipole-Charge
        for(unsigned ej = 0; ej < nej; ej++)
        {
          const Charge& chargej = static_cast<const Charge&>(components[compj].charge(ej));
          double cvalj = chargej.q();
          double minusqmyper4pie0 = -cvalj*absmyi;
          pstrmij << minusqmyper4pie0;
#ifndef NDEBUG
          cout << "21=" << minusqmyper4pie0 << " ";
#endif
        }
        // Dipole-Quadrupole
        for(unsigned int qj=0;qj<nqj;++qj)
        {
          const Quadrupole& quadrupolej=static_cast<const Quadrupole&>(components[compj].quadrupole(qj));
          double absqj=quadrupolej.absQ();
          double myq15=1.5*absmyi*absqj;
          pstrmij << myq15;
#ifndef NDEBUG
          cout << "24=" << myq15 << " ";
#endif
        }
        // Dipole-Dipole
        for(unsigned int dj=0;dj<ndj;++dj)
        {
          const Dipole& dipolej=static_cast<const Dipole&>(components[compj].dipole(dj));
          double absmyj=dipolej.absMy();
          double my2=absmyi*absmyj;
          pstrmij << my2;
#ifndef NDEBUG
          cout << "22a=" << my2 << " ";
#endif
          double rffac=my2*epsRFInvrc3;
          pstrmij << rffac;
#ifndef NDEBUG
          cout << "22b=" << rffac << " ";
#endif
        }
      }
#endif
    }
  }

#ifdef COMPLEX_POTENTIAL_SET
  /*
  // Tersoff interaction parameters
  //
  m_tersofftbl.redim(m_numcomp, m_numcomp);

  for(unsigned int compi = 0; compi < this->m_numcomp; compi++)
  {
    unsigned int nti = components[compi].numTersoff();
    for(unsigned int compj = 0; compj < this->m_numcomp; compj++)
    {
      unsigned int ntj = components[compj].numTersoff();
      ParaStrm& pstrmij = m_tersofftbl(compi, compj);

      for(unsigned int centeri = 0; centeri < nti; centeri++)
      {
        const Tersoff& tersoffi = static_cast<const Tersoff&>(components[compi].tersoff(centeri));
        double Ai = tersoffi.A();
        double Bi = tersoffi.B();
        double minus_lambdai = tersoffi.minusLambda();
        double minus_mui = tersoffi.minusMu();
        double Ri = tersoffi.R();
        double Si = tersoffi.S();
        double cci = tersoffi.cSquare();
        double ddi = tersoffi.dSquare();
        double hi = tersoffi.h();
        double ni = tersoffi.n();
        double betai = tersoffi.beta();

        for(unsigned int centerj = 0; centerj < ntj; centerj++)
        {
          const Tersoff& tersoffj = static_cast<const Tersoff&>(components[compj].tersoff(centerj));
          double Aj = tersoffj.A();
          double Bj = tersoffj.B();
          double minus_lambdaj = tersoffj.minusLambda();
          double minus_muj = tersoffj.minusMu();
          double Rj = tersoffj.R();
          double Sj = tersoffj.S();

          double S = sqrt(Si*Sj);

          pstrmij << sqrt(Ri*Rj)       // R
                  << S                 // S
                  << hi                // h
                  << cci               // c^2
                  << ddi               // d^2
                  << sqrt(Ai*Aj)       // A
                  << -1.0*sqrt(Bi*Bj)  // -B
                  << 0.5 * (minus_lambdai + minus_lambdaj)  // -lambda
                  << 0.5 * (minus_mui + minus_muj)          // -mu
                  << betai             // beta
                  << ni;               // n_i
        }
      }
    }
  }
  */
#endif
}

