/*************************************************************************
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

#include "Wall.h"


using namespace std;
using Log::global_log;


Wall::Wall(){}

Wall::~Wall(){}

void Wall::initialize(const std::vector<Component>& components, 
		 double in_rhoWall, double in_sigWall, double in_epsWall, double* in_xi, double* in_eta, 
		 double in_yOffWall, double in_yWallCut, double in_yMirr)
{
  global_log->info() << "Initializing the wall function.\n";
  _rhoW = in_rhoWall;
  _yc = in_yWallCut; 
  _yOff = in_yOffWall;
  _yMirr = in_yMirr;
  
  
  /*!*** So far: only 1CLJ components allowed ****/
  _nc = components.size();
  _eps_wi = new double[_nc];
  _sig3_wi = new double[_nc];
  _uShift_9_3 = new double[_nc];
  _uPot_9_3 = new double[_nc];
  
  for (unsigned i = 0; i < _nc; i++)
  {
    _eps_wi[i] = in_xi[i]*sqrt(in_epsWall * components[i].getEps(0));
    double sig_wi;
    sig_wi = 0.5*in_eta[i]*(in_sigWall + components[i].getSigma(0));
    _sig3_wi[i] = sig_wi *sig_wi *sig_wi;
    double sig9_sf = _sig3_wi[i] *_sig3_wi[i] *_sig3_wi[i];
    double y = _yc - _yOff;
    double y3 = y*y*y;
    double y9 = y3*y3*y3;
    _uShift_9_3[i] = 4.0/3.0*M_PI*_rhoW*_eps_wi[i]*_sig3_wi[i]*(sig9_sf/15.0/y9 - _sig3_wi[i]/2.0/y3);
    _uPot_9_3[i] = 0.0;
  }
  
    
}

void Wall::calcTSLJ_9_3( ParticleContainer* partContainer, Domain* domain)
{
  
  double regionLowCorner[3], regionHighCorner[3];
  list<Molecule*> particlePtrsForRegion;
  
  /*! LJ-9-3 potential applied in y-direction */
  if(partContainer->getBoundingBoxMin(1) < _yc){ // if linked cell within the potential range (inside the potential's cutoff)
    for(Molecule* currMolec = partContainer->begin(); currMolec != partContainer->end(); currMolec = partContainer->next()){
      //! so far for 1CLJ only, several 1CLJ-components possible
      double y, y3, y9;
      unsigned cid = currMolec -> componentid();
      y = currMolec->r(1) - _yOff;
      if(y < _yc){
	y3 = y*y*y;
	y9 = y3*y3*y3;
	double f[3];
	for(unsigned d = 0; d < 3; d++)
	  f[d] = 0.0;
	
	double sig9_wi;
	sig9_wi = _sig3_wi[cid]*_sig3_wi[cid]*_sig3_wi[cid]; 
	f[1] = 4.0*M_PI* _rhoW * _eps_wi[cid] * _sig3_wi[cid] * ( sig9_wi/5.0/y9 - _sig3_wi[cid]/2.0/y3 ) / y;
	_uPot_9_3[cid] += 4.0*M_PI* _rhoW * _eps_wi[cid] * _sig3_wi[cid] * ( sig9_wi/15.0/y9 - _sig3_wi[cid]/6.0/y3 ) + _uShift_9_3[cid];
	currMolec->Fljcenteradd(0, f);
      } // end if()
    }
  }
  
  double u_pot;
  u_pot = _uPot_9_3[0] + domain -> getLocalUpotCompSpecific();
  domain->setLocalUpotCompSpecific(u_pot);
  for(unsigned cid = 0; cid < _nc; cid++){
    _uPot_9_3[cid] = 0.0;
  }
  
  /*!*** Mirror boundary in y-direction on top of the simulation box****/
  if(partContainer->getBoundingBoxMax(1) > _yMirr){ // if linked cell in the region of the mirror boundary
    for(unsigned d = 0; d < 3; d++){
     regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
     regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
    }
    regionLowCorner[1] = partContainer->getBoundingBoxMin(1) > _yMirr ? partContainer->getBoundingBoxMin(1) : _yMirr;
    partContainer->getRegion(regionLowCorner, regionHighCorner, particlePtrsForRegion);
    
    std::list<Molecule*>::iterator particlePtrIter;
    for(particlePtrIter = particlePtrsForRegion.begin(); particlePtrIter != particlePtrsForRegion.end(); particlePtrIter++){
      if((*particlePtrIter)->v(1) > 0){
	(*particlePtrIter)->setv(1, -(*particlePtrIter)->v(1)); 
      }
    }
    
      
    /*  
    for(Molecule* currMolec = particlePtrsForRegion -> begin(); currMolec != particlePtrsForRegion->end(); currMolec = particlePtrsForRegion->next()){
      if(currMolec->r(1) > _yMirr && currMolec->v(1) > 0){
	currMolec->setv(1, -1.0*currMolec->v(1)); 
      }
    }
    */
    
  } // end Mirror boundary

} // end mthod calcTSLJ_9_3(...)

