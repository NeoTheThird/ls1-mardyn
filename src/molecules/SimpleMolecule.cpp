/***************************************************************************
 *   Copyright (C) 2008 by Martin Bernreuther and colleagues               *
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

#include "molecules/SimpleMolecule.h"
#include <iostream>
#include <cmath>

SimpleMolecule::SimpleMolecule(int id, int type,
                   double xPos, double yPos, double zPos,
                   double xVel, double yVel, double zVel){
  this->id = id;
  this->type = type;
  this->position[0] = xPos;
  this->position[1] = yPos;
  this->position[2] = zPos;
  this->velocity[0] = xVel;
  this->velocity[1] = yVel;
  this->velocity[2] = zVel;
  //this->nextPotentialNeighbour = NULL;
}

int SimpleMolecule::getId(){
  return this->id;
}

int SimpleMolecule::getType(){
  return this->type;
}

double SimpleMolecule::getPosition(int dimension){
  return this->position[dimension];
}

double SimpleMolecule::getVelocity(int dimension){
  return this->velocity[dimension];
}

double SimpleMolecule::getForce(int dimension){
  return this->force[dimension];
}

void SimpleMolecule::setId(int id){
  this->id = id;
}

void SimpleMolecule::setType(int type){
  this->type = type;
}

void SimpleMolecule::setPosition(double position, int dimension){
  this->position[dimension] = position;	
}

void SimpleMolecule::setVelocity(double velocity, int dimension){
  this->velocity[dimension] = velocity;
}

void SimpleMolecule::setForce(double force, int dimension){
  this->force[dimension] = force;
}

double SimpleMolecule::calcDistanceSquare(SimpleMolecule& molecule1, SimpleMolecule& molecule2){
  double distance_square = 0;
  for(int i=0; i<3; i++){
    distance_square += pow(molecule1.position[i]-molecule2.position[i],2);
  }
  return distance_square;
}

void SimpleMolecule::calcForce(SimpleMolecule& molecule1, SimpleMolecule& molecule2){
  double sigma = 1.0;
  double epsilon = 5.0;
  double r_2 = calcDistanceSquare(molecule1, molecule2);
  double s = sigma*sigma / r_2;
  s = pow(s,3);
  double f = 24 * epsilon * s / r_2 * (1-2*s);
  for(int i=0; i<3; i++){
  	double temp = f*(molecule2.position[i]-molecule1.position[i]);
    molecule1.force[i] += temp;
    molecule2.force[i] -= temp;
  }
}
