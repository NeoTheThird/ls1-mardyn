/***************************************************************************
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


/* by Stefan Becker, <stefan.becker@mv.uni-kl.de>*/
//Calculation of the Fluid-Bubbleinteraction 

#ifndef BUBBLE_H_
#define BUBBLE_H_

#include "molecules/Molecule.h"
#include "molecules/Component.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Domain.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>


using namespace std; 

class Bubble{
public:
  
  // constructor and destructor
  Bubble();
 ~Bubble();
  void initialize(double in_A, double in_B, double in_C, double in_D, double in_rc, double* in_BubbleCentre);
  void calcBubbleForce( ParticleContainer* partContainer, Domain* domain );

  
 	
private:
  double _A;
  double _B;
  double _C;
  double _D;
  double _rc2;
  double _bubbleCentre[3];
  double _xCenter;
  double _yCenter;
  double _zCenter;  
};


#endif /*WALL_H_*/
