/*
 * vis_multisphere.h
 *
 *  Created on: 30.05.2012
 *      Author: eck
 */

#ifndef VIS_MULTISPHERE_H_
#define VIS_MULTISPHERE_H_

#include <vector>
#include "molecule.h"
#include "component.h"
#include <iostream>
#include <fstream>

int numberFrames(ifstream input);
int readFrame(ifstream &filestr, vector<Molecule> &vMolecules);
int writeFrameVis(ofstream &filestr, vector<Molecule> vMolecules);
int writeFrameVisMulti(ofstream &filestr, vector<Molecule> vMolecules, vector<Component> vComp);
int writeInit(ofstream &filestr);
int writeIntermis(ofstream &filestr);


#endif /* VIS_MULTISPHERE_H_ */
