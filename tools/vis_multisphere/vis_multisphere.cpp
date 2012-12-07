/*
 * vis_multisphere.cpp
 *
 *  Created on: 30.05.2012
 *      Author: eck
 */

using namespace std;

#include "vis_multisphere.h"
#include <sstream>
#include <string.h>
#include <iomanip>

int main(int argc, char** argv) {
	const char* usage = "usage: vis_multisphere <input vis file>\n\n";
	if(argc != 2){
		cout << usage;
		return 1;
	}
	//double refLength = 3.3211; //n2 (+ ethan)
	double refLength = 1; //Azeton
	vector<Component> vComp;
	Component comp1("azeton", 1, refLength, 1);
	vComp.push_back(comp1);
	//Component comp2("n2", 1, refLength, 1);
	//vComp.push_back(comp2);

	char* infilename = argv[1];
	string outfilename = (string)infilename + "_multi";
	char* outchar = new char[outfilename.size()+1];
	strcpy(outchar, outfilename.c_str());

	//int numFrames = numberFrames(infilename);
	vector<Molecule> vMolecules;

	ifstream input;
	input.open(infilename);
	ofstream output;
	output.open(outchar);

	writeInit(output);

	while(!input.eof()){
		//beschleunigen: hier infilestream übergeben
		readFrame(input,vMolecules);
		//hier outfilestream übergeben; jeden Frame für sich lesen und dann schreiben (oder 10/x Frames)
		writeFrameVisMulti(output,vMolecules,vComp);
		writeIntermis(output);
	}

	input.close();
	output.close();

	/*
	double point[3] = {1, 0, 0};
	double image[3];
	Quaternion quat(-0.552710, -0.441034, -0.604951, 0.366106);
	quat.rotatePoint(point, image);
	cout << image[0] << " " << image[1] << " " << image[2] << "\n";
	*/

	cout << "Fertig!\n";
}

int numberFrames(ifstream input){
	string token;
	input >> token;
	int numberFrames = 1;
	while (!input.eof()){
		if (token == "#"){
			numberFrames++;
		}
		else{
			input.ignore(1024,'\n');
		}
		input >> token;
	}

	return numberFrames;
}

int readFrame(ifstream &filestr, vector<Molecule> &vMolecules){
	vMolecules.clear();

	Molecule mol;
	string line = "";

	string token;
	filestr >> token;
	while(!filestr.eof() && token != "#"){
		if(token == "id"){
			filestr.ignore(1024,'\n');
		}
		else{
			stringstream strstr(token);
			int inputi;
			strstr >> inputi;
			mol.sIdNumber(inputi);
			filestr >> inputi;
			mol.sComponentType(inputi);
			double inputd3[3];
			for (int i = 0; i < 3; i++)
				filestr >> inputd3[i];
			mol.sPosition(inputd3);
			double inputd4[4];
			for (int i = 0; i < 4; i++)
				filestr >> inputd4[i];
			mol.sOrientationQuaternion(inputd4);
			filestr.ignore(1024,'\n');
			vMolecules.push_back(mol);
		}
		filestr >> token;
	}
	return 0;
}

int writeFrameVis(ofstream &filestr, vector<Molecule> vMolecules){
	Molecule mol;
	int numMol = vMolecules.size();
	for (int i = 0; i < numMol; i++){
		mol = vMolecules[i];
		filestr << setw(8) << mol.gIdNumber() << " " << mol.gComponentType() << " "
				<< setw(10) << mol.gPosition()[0] << " " << setw(10) << mol.gPosition()[1] << " " << setw(10) << mol.gPosition()[2] << " "
				<< setw(6) << mol.gOrientationQuaternion().gCoordinates()[0] << " " << setw(6) << mol.gOrientationQuaternion().gCoordinates()[1] << " " << setw(6) << mol.gOrientationQuaternion().gCoordinates()[2] << " " << setw(6) << mol.gOrientationQuaternion().gCoordinates()[3] << " "
				<< setw(8)<< "0" << "\n";
	}
	return 0;
}

int writeFrameVisMulti(ofstream &filestr, vector<Molecule> vMolecules, vector<Component> vComp){
	Molecule mol;
	int numMol = vMolecules.size();
	for (int i = 0; i < numMol; i++){
		mol = vMolecules[i];
		int id = mol.gIdNumber(); //wird für die einzelnen LJ-Sites mit der LJID multipliziert
		int compId = mol.gComponentType();
		double pos[3] = {mol.gPosition()[0], mol.gPosition()[1], mol.gPosition()[2]};
		Quaternion quat = mol.gOrientationQuaternion();
		Component comp = vComp[compId];
		int numLJs = comp.gnumberLJCenters();
		vector<float> vLJx = comp.gvLJx();
		vector<float> vLJy = comp.gvLJy();
		vector<float> vLJz = comp.gvLJz();
		for (int j = 0; j < numLJs; j++){
			int idFaktor = j+1;
			double posLJ[3] = {(double)vLJx[j], (double)vLJy[j], (double)vLJz[j]};
			double newPos[3];
			quat.rotatePoint(posLJ,newPos);
			for (int k = 0; k < 3; k++)
				newPos[k] += pos[k];
			int ljId = 0;
			for (int k = 0; k < compId; k++)
				ljId += vComp[k].gnumberLJCenters();
			ljId += j;
			filestr.precision(3);
			filestr << setw(8) << idFaktor*id << " " << ljId << " "
					<< fixed << setw(10) << newPos[0] << " " << setw(10) << newPos[1] << " " << setw(10) << newPos[2] << " "
					<< setw(6) << quat.gCoordinates()[0] << " " << setw(6) << quat.gCoordinates()[1] << " " << setw(6) << quat.gCoordinates()[2] << " " << setw(6) << quat.gCoordinates()[3] << " "
					<< setw(8)<< "0" << "\n";
		}
	}
	return 0;
}

int writeInit(ofstream &filestr){
	filestr << "      id t          x          y          z     q0     q1     q2     q3        c\n";
	return 0;
}

int writeIntermis(ofstream &filestr){
	filestr << "#\n";
	return 0;
}
