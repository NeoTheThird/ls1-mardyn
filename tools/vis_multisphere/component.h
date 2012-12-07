/*
 * component.h
 *
 *  Created on: 27.01.11
 *      Author: eck
 */

#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <string>
#include <vector>
#include <cmath> //für pow() und sqrt()


using namespace std;

class Component{
private:
	int _numberLJCenters;
	int _numberCharges;
	int _numberQuadrupoles;
	//wenn COMPLEX_POTENTIAL_SET braucht man auch
	int _numberDipoles;
	int _numberTersoff; //Modell um Feststoffe zu simulieren, fokussiert die Verbindungen zwischen Atomen
	//LJ center position
	vector<float> _vLJx;
	vector<float> _vLJy;
	vector<float> _vLJz;
	//LJ center Masse
	vector<float> _vLJm;
	//LJ center Energie
	vector<float> _vLJeps;
	//LJ center Länge
	vector<float> _vLJsigma;
	//Dipol Position
	vector<float> _vdipoleX;
	vector<float> _vdipoleY;
	vector<float> _vdipoleZ;
	//Dipolmoment
	vector<float> _vdipoleMyx;
	vector<float> _vdipoleMyy;
	vector<float> _vdipoleMyz;
	vector<float> _vdipoleMyabs;
	//Quadrupol Position
	vector<float> _vquadrupoleX;
	vector<float> _vquadrupoleY;
	vector<float> _vquadrupoleZ;
	//Quadrupolmoment
	vector<float> _vquadrupoleQx;
	vector<float> _vquadrupoleQy;
	vector<float> _vquadrupoleQz;
	vector<float> _vquadrupoleQabs;
	//Trägheitsmoment
	double ixx; //in Richtung der x-Achse
	double iyy;
	double izz;
	//allgemeine Parameter
	double eps;
	double sigma;
	double mass;
	double length;
	double qdr;
	double cutoff;
	int resize();
public:
	Component(string substance, double refEnergy, double refLength, double refMass);
	int gnumberLJCenters() const;
	double gMass() const;
	double gIxx() const;
	double gIyy() const;
	double gIzz() const;
	vector<float> gvLJx() const;
	vector<float> gvLJy() const;
	vector<float> gvLJz() const;
/*
	vector<float> gvquadrupoleX() const;
	vector<float> svquadrupoleX(int ersetz);
	vector<float> gvquadrupoleY() const;
	vector<float> svquadrupoleY(int ersetz);
	vector<float> gvquadrupoleZ() const;
	vector<float> svquadrupoleZ(int ersetz);
	vector<float> gvquadrupoleQx() const;
	vector<float> svquadrupoleQx(int ersetz);
	vector<float> gvquadrupoleQy() const;
	vector<float> svquadrupoleQy(int ersetz);
	vector<float> gvquadrupoleQz() const;
	vector<float> svquadrupoleQz(int ersetz);
	vector<float> gvquadrupoleQabs() const;
	vector<float> svquadrupoleQabs(int ersetz);
*/
	vector<float> gComponent() const;
	unsigned size();
	float operator[](int index);
	friend ostream& operator<< (ostream& os, Component& comp);
};

#endif /* COMPONENT_H_ */
