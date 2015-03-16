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

/* by Stefan Becker, <stefan.becker@mv.uni-kl.de>*/
 
 
#include "Bubble.h"


using namespace std;
using Log::global_log;


Bubble::Bubble(){}
Bubble::~Bubble(){}

void Bubble::initialize(double in_A, double in_B, double in_C, double in_D, double in_rc, double* in_BubbleCentre){
  _A = in_A;
  _B = in_B;
  _C = in_C;
  _D = in_D;
  _rc2 = in_rc*in_rc;
  for(unsigned i = 0; i<3; i++){
    _bubbleCentre[i] = in_BubbleCentre[i];
  }
  global_log->info() << "Initializing the Bubble function.\n";
//  global_log->info() << "Bubble centre: " << _bubbleCentre[0] << _bubbleCentre[1] << _bubbleCentre[2]<< "\n";
}

void Bubble::calcBubbleForce( ParticleContainer* partContainer, Domain* domain)
{
  
  double regionLowCorner[3], regionHighCorner[3];
  list<Molecule*> particlePtrsForRegion;
  
  double rFromCOB2[3];
  
  for(unsigned i = 0; i< 3; i++){
    double helpDist1 = (partContainer->getBoundingBoxMin(i) - _bubbleCentre[i])*(partContainer->getBoundingBoxMin(i) - _bubbleCentre[i]);
    double helpDist2 = (partContainer->getBoundingBoxMax(i) - _bubbleCentre[i])*(partContainer->getBoundingBoxMax(i) - _bubbleCentre[i]); 
    rFromCOB2[i] = helpDist1 < helpDist2 ? helpDist1 : helpDist2;
  }
    
//  cout<< "Calc Bubble Fun\ndist2 = " << rFromCOB2[0] + rFromCOB2[1] + rFromCOB2[2]<< "\n";
    
  if( ( rFromCOB2[0] + rFromCOB2[1] + rFromCOB2[2] < _rc2) ){ // if linked cell within the potential range (inside the potential's cutoff)
      for(unsigned d = 0; d < 3; d++){
	regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
	regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
	}
      partContainer->getRegion(regionLowCorner, regionHighCorner, particlePtrsForRegion);
      
      std::list<Molecule*>::iterator particlePtrIter;
      
      for(particlePtrIter = particlePtrsForRegion.begin(); particlePtrIter != particlePtrsForRegion.end(); particlePtrIter++){
	double rx,ry,rz;
	double r2;
	unsigned cid = (*particlePtrIter) -> componentid();
	if(cid == 0){
	  rx = (*particlePtrIter)->r(0) - _bubbleCentre[0];
	  ry = (*particlePtrIter)->r(1) - _bubbleCentre[1];
	  rz = (*particlePtrIter)->r(2) - _bubbleCentre[2];
	  r2 = rx*rx + ry*ry + rz*rz;
	  if(r2 < _rc2){
	    double normR = sqrt(r2);
	    double r3 = r2*normR;
	    double r4 = r2*r2;
	    double r11 = r4*r4*r3;
	    double r12 = r11*normR;
	    double f[3];
	    for(unsigned d = 0; d < 3; d++)
	      f[d] = 0.0;
	    double normF = ( 3.*_B*r2*(_C + _D*r12) + 12.*_D*r11*(_A - _B*r3) ) / (_C + _D*r12) / (_C + _D*r12);
	    f[0] = rx/normR *normF;
	    f[1] = ry/normR *normF;
	    f[2] = rz/normR *normF;
	    (*particlePtrIter)->Fljcenteradd(0, f);
	  } // end if (r2<rc*rc)
      } // end if(cid==0)
    }
  }
    particlePtrsForRegion.clear();
} // end mthod calcTSLJ_9_3(...)

double Bubble::getCentre(unsigned direction){
  return _bubbleCentre[direction];
}
