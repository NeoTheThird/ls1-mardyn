/*
 * PhaseSpaceWriter.cpp
 *
 *  Created on: 21.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 */

#include"PhaseSpaceWriter.h"

extern double LATTICE_CONST_WALL_LJTS;
extern const string WALL_CU_LJ;

// constructor

PhaseSpaceWriter::PhaseSpaceWriter(string in_prefix, double in_Temperature, double in_densFac, unsigned in_nFluid, string in_fluidComponent,
		string in_wallComponent, unsigned in_wallLayers, double in_xi12, double in_xi13, double in_eta12, double in_alpha, double in_beta, double in_gamma, double in_edgeProp, bool in_stripes,
		unsigned in_numberOfStripes, bool in_LJShifted, bool in_LJunits)
{
	cout << "\n\n**********************************\nPhaseSpaceWriter opened\n";
	cout << "writing the header ...\n**********************************\n";
	_fileName = in_prefix + ".template.inp";
	_fileNameXyz = in_prefix+".xyz";
	_fluidComponentName = in_fluidComponent;
	_wallComponentName = in_wallComponent;
	_wallLayers = in_wallLayers;
	_nFluid = in_nFluid;
	_temperature = in_Temperature;
	_densFac = in_densFac;
	_xi12 = in_xi12;
	_xi13 = in_xi13;
	_eta12 = in_eta12;
	_alpha = in_alpha;
	_beta = in_beta;
	_gamma = in_gamma;
	_edgeProp = in_edgeProp;
	_stripes = in_stripes;
	_numberOfStripes = in_numberOfStripes;
	_LJShifted = in_LJShifted;
	_LJunits = in_LJunits;
//	cout << "\n**********************************\nPhaseSpaceWriter constructor finished\n**********************************\n";

/*	// generating the ofstream, i.e. the phase space file
	psstrm(_fileName.c_str(), ios::binary|ios::out);
	if(!psstrm){ // simple error handling: ofstream MUST exist
		cerr << "Error in constructor of PhaseSpaceWriter:\nPhase space file \""<< _fileName <<"\" cannot be opened.";
		exit(-50);
	}
*/
}

//destructor
PhaseSpaceWriter::~PhaseSpaceWriter(){
};

void PhaseSpaceWriter::write()
{
//	cout <<"\n**********************************\n PhaseSpaceWriter write() method started\n**********************************\n";
	//@brief: wirting the header
	// generating the ofstream, i.e. the phase space file.
	 ofstream psstrm(_fileName.c_str(), ios::binary|ios::out);
		if(!psstrm){ // simple error handling: ofstream MUST exist
			cerr << "Error in 'write()'-method of PhaseSpaceWriter:\nPhase space file \""<< _fileName <<"\" cannot be opened.";
			exit(101);
		}
	// xyz file format, to be used for VMD-visualisation
	ofstream xyzstrm(_fileNameXyz.c_str(), ios::binary|ios::out);
	if(!xyzstrm){ // simple error handling: ofstream MUST exist
				cerr << "Error in constructor of PhaseSpaceWriter:\nPhase space file \""<< _fileNameXyz <<"\" cannot be opened.";
				exit(101);
			}

	Component fluidComp(_fluidComponentName, _LJunits);
	GlobalStartGeometry geometry(_nFluid, fluidComp.calculateLiquidDensity(_temperature),
									fluidComp.calculateVaporDensity(_temperature, _densFac), _alpha, _beta, _gamma);
//	cout << "\n**********************************\n objects of classes 'Component' and 'GlobalStartGeometry' generated\n**********************************\n";
	unsigned numberOfComponents;
	if(_stripes){
		numberOfComponents = 3;
	}
	else {
		numberOfComponents = 2;
	}

	double hWall, latticeConst;
	if(_wallComponentName == WALL_CU_LJ){
		latticeConst= LATTICE_CONST_WALL_LJTS;
	}
	else{
		cerr << "error in PhaseSpaceWriter.write(): no lattice const available for this wall model.";
		exit(102);
	}
	hWall = latticeConst * (_wallLayers-0.5);
	double shielding = fluidComp.gSigma(0);
	geometry.calculateBoxFluidOffset(hWall, shielding, _edgeProp);

	psstrm << "mardyn trunk 20120208\n#mardyn input file, ls1 project\n#generated by the mkSD tool\n#lattice constant solid: "<<LATTICE_CONST_WALL_LJTS<< "\n";
	psstrm << "#-a " << _alpha << " \t -b " << _beta << " \t-g " << _gamma << " \t-d " << _densFac << "\n";
	psstrm << "# gross fluid density: " << geometry.gGrossFluidDens() << "\n";
	psstrm << "t\t0.0\n";
	// assignments of thermostats
	if(_stripes){
		psstrm << "CT\t2 1\n";
		psstrm << "CT\t3 2\n";
		psstrm << "ThT\t1 "<< _temperature << "\n";
		psstrm << "ThT\t2 "<< _temperature << "\n";
		psstrm << "CT\t1 3\nThT\t3 "<< _temperature << "\n";
	}
	else{
		psstrm << "CT\t2 1\n";
		psstrm << "ThT\t1 "<< _temperature << "\n";
		psstrm << "CT\t1 2\nThT\t2 "<< _temperature << "\n";
	}
	_boxLengthY = geometry.gBoxLength(1);
	psstrm << "L\t" << geometry.gBoxLength(0)<<"\t"<<geometry.gBoxLength(1) <<"\t" << geometry.gBoxLength(2) <<"\n";
	psstrm << "C\t"<< numberOfComponents <<"\n";
	psstrm << fluidComp.gNumberLJCenters() << " "<<fluidComp.gNumberCharges() <<" "<<fluidComp.gNumberQuadrupoles() <<" "
			<<fluidComp.gNumberDipoles() <<" "<<fluidComp.gNumberTersoff()<<"\n";
//***************************************************************************************************************************************
	// writing the parameters of the LJ-fluid
	for (unsigned i = 0; i < fluidComp.gNumberLJCenters(); i++){
		psstrm << "0.0 0.0 0.0\t";  // center of mass in the local coordinate system
		psstrm << fluidComp.gMass(i) <<" "<<fluidComp.gEps(i) << " " << fluidComp.gSigma(i) <<" " <<fluidComp.gRCutLJ() <<" "<< _LJShifted<<"\t";
	}
	psstrm << "0.0 0.0 0.0\n";  // Hauptträgheitsachsen des Gesamtmoleküls
//***************************************************************************************************************************************
	// writing the parameters of the LJ type wall
	Component wallComp(_wallComponentName, _LJunits);
	for(unsigned i = 0; i < numberOfComponents-1; i++){ // for-loop counter depends on wheter or not there's a stripes shaped wall
		psstrm << wallComp.gNumberLJCenters() << " " << wallComp.gNumberCharges() << " " << wallComp.gNumberQuadrupoles() << " "
				<< wallComp.gNumberDipoles() << " " << wallComp.gNumberTersoff() << "\n";
		// @brief: The LJ type wall has the same cut off radius as the fluid component, since ls1 is capable of only handling
		// a sigle LJ cut off radius
		for(unsigned i = 0; i < wallComp.gNumberLJCenters(); i++){
			psstrm << "0.0 0.0 0.0\t"; // center of mass in the local coordinate system
			psstrm << wallComp.gMass(i) <<" "<< wallComp.gEps(i) << " " << wallComp.gSigma(i) << " "<< wallComp.gRCutLJ() << " " << _LJShifted << "\t";
		}
		psstrm << "0.0 0.0 0.0\n"; // Hauptträgheitsachsen des Gesamtmoleküls
	}

//	cout << "\n**********************************\nobject 'wallComp' of the 'Component' class generated.\n**********************************\n";

	// so far no other models than LJ implemented => no scheme for writing the model parameters corresponding to charges, dipoles, etc.
//***************************************************************************************************************************************
	// writing Lorentz Berthelot xi and eta: each line contains (xi,eta) of one particular component with all the other components, e.g xi12 eta12 xi13 eta13 ...
	//@todo: Zweites xi muss vorgegeben werden können!!!
	if(_stripes){
		psstrm <<  _xi12 << " "<<_eta12 <<"\t "<<  _xi13 << " "<<_eta12 <<"\n1 1\n";
	}
	else{
		psstrm << _xi12 << " "<<_eta12 <<"\n";
	}
	// writing the time parameter for reaction  field method
	psstrm << "1.0e+10\n";
//***************************************************************************************************************************************

	// The call of _wallMolecule with the second argument wallComp.gMass(0) only works for the copper LJ-model!!! Multicenter LJ-models have mor than one single mass (gMass(0))!!!
	Molecule wallMolecule(_temperature, wallComp.gMass(0));


//	cout << "\n**********************************\nobject walMolecule of the class 'Molecule' generated\n**********************************\n";

	/*wallMolecule.calculateCoordinatesOfWallMolecule(_wallLayers, geometry.gBoxLength(0), geometry.gBoxLength(1), geometry.gBoxLength(2),
			geometry.gOffset(0), geometry.gOffset(1), geometry.gOffset(2), latticeConst);*/
	if(_stripes){
		wallMolecule.calculateCoordinatesOfWallMolecule(geometry.gBoxLength(0),geometry.gBoxLength(2),
										   geometry.gOffsetLiq(1), latticeConst, shielding, _numberOfStripes);
	}
	else{
		wallMolecule.calculateCoordinatesOfWallMolecule(geometry.gBoxLength(0),geometry.gBoxLength(2),
															geometry.gOffsetLiq(1), latticeConst, shielding);
	}
//	cout<< "\n**********************************\n coordinates of wall molecules calculated\n**********************************\n";

	// total number of particles
	geometry.calculateLiqFillProbabilityArray(); // within this method the number of actually filled liquid particles is calculated, too => already called here
	geometry.calculateVapFillProbabilityArray();
	psstrm << "N\t"<< wallMolecule.gNumberOfMolecules() + geometry.gNFilledLiqSlots() + geometry.gNFilledVapSlots() <<"\n";
	xyzstrm << wallMolecule.gNumberOfMolecules() + geometry.gNFilledLiqSlots() + geometry.gNFilledVapSlots()<<"\nComment\n";
	psstrm << "M\tICRVQD\n\n";

/***************************************************************************************************************************************************
 ****************************************************************************************************************************************************
 ****************************************************************************************************************************************************
 */

	// @brief: wirting the body, i.e. the particle ID, component-ID, position in the simulation box, velocities, quaternions and angular velocities
	unsigned cid = 1;	// component ID
	unsigned id = 1;	// molecule ID

	// liquid phase slots are filled
	unsigned ii[4];		// four different counters used in the for loops following
	double positionVec[3];


	//cout << "fluidUnit[0]: "<<geometry.gFluidUnit(0) << " fluidUnits[0]: " << geometry.gFluidUnits(0) <<"\n";
	//cout << "fluidUnit[1]: "<<geometry.gFluidUnit(1) << " fluidUnits[1]: " << geometry.gFluidUnits(1) <<"\n";
	//cout << "fluidUnit[2]: "<<geometry.gFluidUnit(2) << " fluidUnits[2]: " << geometry.gFluidUnits(2) <<"\n";
	cout << "\n**********************************\nWriting the body ...\n**********************************\n";
	RandomNumber rdm;
	//liquid particles
	for(ii[0] = 0; ii[0] < geometry.gLiqUnits(0); ii[0]++){
		for(ii[1] = 0; ii[1] < geometry.gLiqUnits(1); ii[1]++){
			for(ii[2] = 0; ii[2] < geometry.gLiqUnits(2); ii[2]++){
				for(ii[3] = 0; ii[3] < 3; ii[3]++){			// three slots per elementary fluid box
					if(geometry.gFillLiqArray(ii[0], ii[1], ii[2], ii[3])){
						for(unsigned short j = 0; j < 3; j++){
							positionVec[j] =  geometry.gOffsetLiq(j)+geometry.gLiqUnit(j)*(ii[j] + 0.02*rdm.randNum() + ((j==ii[3])? 0.24 : 0.74));
						}	// end for(j...)
						for(unsigned j = 0; j < 3; j++){ // in case the molecule is placed outside the simulation box it is shifted (according to periodic boundary conditions)
							if(j == 0 || j == 2){
								if(geometry.gBoxLength(j) < positionVec[j]) positionVec[j] -= geometry.gBoxLength(j);
								else if(positionVec[j] < 0) positionVec[j] += geometry.gBoxLength(j);
							}
							else{ // PBC in y-direction => collision with the wall!
								if(geometry.gBoxLength(j) < positionVec[j]){
									cerr << "Severe error in PhaseSpaceWriter::write() => writing the fluid positions:\n"
											<< "Fluid particle placed within wall due to PBC in y-direction!!!\n";
											exit(104);
								}
								else if(positionVec[j] < 0) positionVec[j] += geometry.gBoxLength(j);
							}
						}// end for(j...)
						double absVelocity = sqrt(3.0*_temperature/fluidComp.gMass(0));
						double phi = 2.0*PI*rdm.randNum();
						double omega = 2.0*PI*rdm.randNum();
						psstrm << id << " " << cid <<" \t" << positionVec[0] << " " << positionVec[1] <<" " << positionVec[2]
						        << " \t" << absVelocity * cos(phi)*cos(omega) << " " << absVelocity *cos(phi)*sin(omega)<< " " << absVelocity *sin(phi)
								<<"\t 1.0 0.0 0.0 0.0\t0 0 0\n";
						xyzstrm << "H \t" << positionVec[0] << " \t" << positionVec[1] <<" \t"<< positionVec[2] << "\n";
						id++;
					}	// end if(_fill)
					else psstrm << "\n";
				}	// end for ii[4]
			}	// end for ii[2]
		} 	// end for ii[1]
	}	// end for ii[0]
	//cid ++;
	
	
	for(unsigned i = 0; i<4;i++) // resetting the counter to zero
	  ii[i] = 0;
	// vapour particles
	for(ii[0] = 0; ii[0] < geometry.gVapUnits(0); ii[0]++){
		for(ii[1] = 0; ii[1] < geometry.gVapUnits(1); ii[1]++){
			for(ii[2] = 0; ii[2] < geometry.gVapUnits(2); ii[2]++){
				for(ii[3] = 0; ii[3] < 3; ii[3]++){			// three slots per elementary fluid box
					if(geometry.gFillVapArray(ii[0], ii[1], ii[2], ii[3])){
						for(unsigned short j = 0; j < 3; j++){
							positionVec[j] =  geometry.gOffsetVap(j)+geometry.gVapUnit(j)*(ii[j] + 0.02*rdm.randNum() + ((j==ii[3])? 0.24 : 0.74));
						}	// end for(j...)
						for(short j = 0; j < 3; j++){ // in case the molecule is placed outside the simulation box it is shifted (according to periodic boundary conditions)
							if(j == 0 || j == 2){
								if(geometry.gBoxLength(j) < positionVec[j]) positionVec[j] -= geometry.gBoxLength(j);
								else if(positionVec[j] < 0) positionVec[j] += geometry.gBoxLength(j);
							}
							else{ // PBC in y-direction => collision with the wall!
								if(geometry.gBoxLength(j) < positionVec[j]){
									cerr << "Severe error in PhaseSpaceWriter::write() => writing the fluid positions:\n"
											<< "Fluid particle placed within wall due to PBC in y-direction!!!\n";
											exit(104);
								}
								else if(positionVec[j] < 0) positionVec[j] += geometry.gBoxLength(j);
							}
						}// end for(j...)
						double absVelocity = sqrt(3.0*_temperature/fluidComp.gMass(0));
						double phi = 2.0*PI*rdm.randNum();
						double omega = 2.0*PI*rdm.randNum();
						psstrm << id << " " << cid <<" \t" << positionVec[0] << " " << positionVec[1] <<" " << positionVec[2]
						        << " \t" << absVelocity * cos(phi)*cos(omega) << " " << absVelocity *cos(phi)*sin(omega)<< " " << absVelocity *sin(phi)
								<<"\t 1.0 0.0 0.0 0.0\t0 0 0\n";
						xyzstrm << "C \t" << positionVec[0] << " \t" << positionVec[1] <<" \t"<< positionVec[2] << "\n";
						id++;
					}	// end if(_fill)
					else psstrm << "\n";
				}	// end for ii[4]
			}	// end for ii[2]
		} 	// end for ii[1]
	}	// end for ii[0]
	

	// wall molecules are being filled
	wallMolecule.calculateVelocities();
	unsigned numberOfWallMolecules = wallMolecule.gNumberOfMolecules();
	cout << "Number of fluid molecules: " << geometry.gNFilledLiqSlots() + geometry.gNFilledVapSlots()<< "\n";
	cout << "Number of wall molecules: "<< numberOfWallMolecules <<"\n**********************************\n\n";
	for (unsigned i = 0; i < numberOfWallMolecules; i++){
		psstrm << id <<" "<< wallMolecule.gMoleculeCID(i) <<"\t"<< wallMolecule.gXPos(i) <<" "<< wallMolecule.gYPos(i) <<" "<< wallMolecule.gZPos(i) <<
				"\t"<< wallMolecule.gXVelocity(i) <<" "<< wallMolecule.gYVelocity(i) <<" "<< wallMolecule.gZVelocity(i)
				<<"\t1.0 0.0 0.0 0.0 \t0 0 0\n";
		if(wallMolecule.gMoleculeCID(i) == 2){
			xyzstrm << "O \t"<< wallMolecule.gXPos(i) <<" "<< wallMolecule.gYPos(i) <<" "<< wallMolecule.gZPos(i) <<"\n";
		}
		else{
			xyzstrm << "N \t"<< wallMolecule.gXPos(i) <<" "<< wallMolecule.gYPos(i) <<" "<< wallMolecule.gZPos(i) <<"\n";
		}
	if (wallMolecule.gMoleculeCID(i) != 3 && wallMolecule.gMoleculeCID(i) != 2 ){
			cerr << "!!!Error: wall cid differs from 2 or 3, respectively! cid = " << wallMolecule.gMoleculeCID(i) << "\n";
			cerr << "id = " << id <<"\n";
		}
		id++;
	} // end for(i...)
	psstrm.close();
	xyzstrm.close();
}

double PhaseSpaceWriter::gBoxLengthY(){
	return _boxLengthY;
}
