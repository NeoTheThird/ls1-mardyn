/*
 * component.cpp
 *
 *  Created on: 09.05.2011
 *      Author: eck
 */

#include <iostream>

using namespace std;

#include "component.h"

/*
 * umrechnen von atomic units in A°, K,
 * Unit length             1 = 1 a0 (Bohr radius) = 0.0529177 nm
 * Elementary charge       1 = 1 e = 96485.3 C/mol //bezogen auf 1 mol Stoff
 * Elementary charge       1 = 1 e = 1.602176487(40)e−19
 * Unit mass               1 = 1 kg/mol = 1000 u //bezogen auf 1 mol Stoff, eigentlich s.u.
 * Unit mass               1 = 1 kg/mol = 1000 u/mol
 * Unit mass               1 = 1.660538782(83)e-24 kg
 *
 */

Component::Component(string substance, double refEnergy, double refLength, double refMass){
	if (substance == "n2"){
		_numberLJCenters = 2;
		_numberCharges = 0;
		_numberQuadrupoles = 1;
		_numberDipoles = 0;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen-0.61537
		_vLJx[0] = 0;
		_vLJy[0] = 0;
		_vLJz[0] = -0.5232/refLength;
		_vLJm[0] = 14.007/refMass;
		_vLJeps[0] = 34.897/refEnergy;
		_vLJsigma[0] = 3.3211/refLength;
		_vLJx[1] = 0;
		_vLJy[1] = 0;
		_vLJz[1] = 0.5232/refLength;
		_vLJm[1] = 14.007/refMass;
		_vLJeps[1] = 34.897/refEnergy;
		_vLJsigma[1] = 3.3211/refLength;
		_vquadrupoleX[0] = 0;
		_vquadrupoleY[0] = 0;
		_vquadrupoleZ[0] = 0;
		_vquadrupoleQx[0] = 0;
		_vquadrupoleQy[0] = 0;
		_vquadrupoleQz[0] = 1;
		//vquadrupoleQabs[0] = -1.07038/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//_vquadrupoleQabs[0] = -122.5263/(pow(refLength,2)*sqrt(refLength*refEnergy)); //*0.5291772^2*sqrt(0.5291772)*sqrt(315774.5)=*114.46991
		_vquadrupoleQabs[0] = 1.4397*85.105589165/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//Trägheitsmoment
		ixx = 7.6685/(pow(refLength,2)*refMass);//\Sum_Einzelmassen mass_i * r_i^2//TODO etwas falsch
		iyy = 7.6685/(pow(refLength,2)*refMass);//Referenzwert: 1/(pow(refLength,2)*refMass) = 1/refEnergy/refTime/refTime
		izz = 0;
		//allgemeine Parameter
		//eps = 34.897/refEnergy; //*315774.5
		//sigma = 3.3211/refLength;//*0.5291772
		mass = 28.014/refMass;//*1000
		//length = 1.0464/refLength;//*0.5291772
		//qdr = -1.07038/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//cutoff = 38.150/refLength;
	}
	else if (substance == "c2h6"){//aus Vrabec01
		//doppelt überprüft; (fast) gleiche Werte wie bei Martin
		_numberLJCenters = 2;
		_numberCharges = 0;
		_numberQuadrupoles = 1;
		_numberDipoles = 0;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = 0;
		_vLJy[0] = 0;
		_vLJz[0] = -1.1881/refLength;
		_vLJm[0] = 15.0347/refMass;
		_vLJeps[0] = 136.99/refEnergy;
		_vLJsigma[0] = 3.4896/refLength;
		_vLJx[1] = 0;
		_vLJy[1] = 0;
		_vLJz[1] = 1.1881/refLength;
		_vLJm[1] = 15.0347/refMass;
		_vLJeps[1] = 136.99/refEnergy;
		_vLJsigma[1] = 3.4896/refLength;
		_vquadrupoleX[0] = 0;
		_vquadrupoleY[0] = 0;
		_vquadrupoleZ[0] = 0;
		_vquadrupoleQx[0] = 0;
		_vquadrupoleQy[0] = 0;
		_vquadrupoleQz[0] = 1;
		_vquadrupoleQabs[0] = -0.8277*85.105589165/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//Trägheitsmoment
		ixx = 42.445412/(pow(refLength,2)*refMass);
		iyy = 42.445412/(pow(refLength,2)*refMass);
		izz = 0;
		//allgemeine Parameter
		//eps = 136.99/refEnergy;
		//sigma = 6.59439/refLength;
		mass = 30.0694/refMass;
		//length = 2.3762014/refLength;
		//qdr = -0.61537/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//cutoff = 40.689/refLength;
	}
	else if (substance == "1clj"){
		_numberLJCenters = 1;
		_numberCharges = 0;
		_numberQuadrupoles = 0;
		_numberDipoles = 0;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = 0;
		_vLJy[0] = 0;
		_vLJz[0] = 0;
		_vLJm[0] = 1/refMass;
		_vLJeps[0] = 1/refEnergy;
		_vLJsigma[0] = 1/refLength;
		//Trägheitsmoment
		ixx = 0;
		iyy = 0;
		izz = 0;
		//allgemeine Parameter
		mass = 1/refMass;
	}
	else if (substance == "2clj_l*1"){
		_numberLJCenters = 2;
		_numberCharges = 0;
		_numberQuadrupoles = 0;
		_numberDipoles = 0;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = -0.5/refLength;
		_vLJy[0] = 0/refLength;
		_vLJz[0] = 0/refLength;
		_vLJx[1] = 0.5/refLength;
		_vLJy[1] = 0/refLength;
		_vLJz[1] = 0/refLength;
		_vLJm[0] = 1/refMass;
		_vLJm[1] = 1/refMass;
		_vLJeps[0] = 1/refEnergy;
		_vLJeps[1] = 1/refEnergy;
		_vLJsigma[0] = 1/refLength;
		_vLJsigma[1] = 1/refLength;
		//Trägheitsmoment
		ixx = 0;
		iyy = 0;
		izz = 0;
		//allgemeine Parameter
		mass = 1/refMass;
	}
	else if (substance == "ar"){
		_numberLJCenters = 1;
		_numberCharges = 0;
		_numberQuadrupoles = 0;
		_numberDipoles = 0;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = 0;
		_vLJy[0] = 0;
		_vLJz[0] = 0;
		_vLJm[0] = 39.948/refMass;
		_vLJeps[0] = 116.79/refEnergy;
		_vLJsigma[0] = 3.3952/refLength;
		//Trägheitsmoment
		ixx = 0;
		iyy = 0;
		izz = 0;
		//allgemeine Parameter
		mass = 39.948/refMass;
	}
	else if (substance == "eox"){//aus eckl08
		//doppelt überprüft; wie bei Martin, nur bei ihm ist 0 nicht der Ladungsschwerpunkt wie hier, sondern der Massenschwerpunkt
		_numberLJCenters = 3;
		_numberCharges = 0;
		_numberQuadrupoles = 0;
		_numberDipoles = 1;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = 0.78/refLength;
		_vLJy[0] = 0;
		_vLJz[0] = -0.48431/refLength;
		_vLJm[0] = 14.0268/refMass;
		_vLJeps[0] = 84.739/refEnergy;
		_vLJsigma[0] = 3.5266/refLength;
		_vLJx[1] = -0.78/refLength;
		_vLJy[1] = 0;
		_vLJz[1] = -0.48431/refLength;
		_vLJm[1] = 14.0268/refMass;
		_vLJeps[1] = 84.739/refEnergy;
		_vLJsigma[1] = 3.5266/refLength;
		_vLJx[2] = 0;
		_vLJy[2] = 0;
		_vLJz[2] = 0.73569/refLength;
		_vLJm[2] = 15.999/refMass;
		_vLJeps[2] = 62.126/refEnergy;
		_vLJsigma[2] = 3.0929/refLength;
		_vdipoleX[0] = 0;
		_vdipoleY[0] = 0;
		_vdipoleZ[0] = 0;
		_vdipoleMyx[0] = 0;
		_vdipoleMyy[0] = 0;
		_vdipoleMyz[0] = 1;
		_vdipoleMyabs[0] = -2.459 * 85.105589165/refLength/sqrt(refLength*refEnergy);//TODO warum negativ? ls1 anders definiert als ms2
		//Trägheitsmoment
		ixx = 15.23944/(pow(refLength,2)*refMass);
		//iyy = 18.95885/(pow(refLength,2)*refMass);// TODO: falsch
		iyy = 32.30725/(pow(refLength,2)*refMass);
		izz = 17.06781/(pow(refLength,2)*refMass);
		//allgemeine Parameter
		mass = 44.0526/refMass;
	}
	else if (substance == "azeton"){//von Jadran geschickt 22.02.2012
		_numberLJCenters = 4;
		_numberCharges = 0;
		_numberQuadrupoles = 1;
		_numberDipoles = 1;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = 0;
		_vLJy[0] = 0;
		_vLJz[0] = 0;
		_vLJm[0] = 12/refMass;
		_vLJeps[0] = 9.82159997/refEnergy;
		_vLJsigma[0] = 2.93069681319184/refLength;
		_vLJx[1] = 0;
		_vLJy[1] = 1.20953190922449/refLength;
		_vLJz[1] = 0;
		_vLJm[1] = 16/refMass;
		_vLJeps[1] = 106.9872519625/refEnergy;
		_vLJsigma[1] = 3.37041500923298/refLength;
		_vLJx[2] = 0;
		_vLJy[2] = -0.803097279216327/refLength;
		_vLJz[2] = 1.28531461746939/refLength;
		_vLJm[2] = 15/refMass;
		_vLJeps[2] = 111.97954865097/refEnergy;
		_vLJsigma[2] = 3.62249459146099/refLength;
		_vLJx[3] = 0;
		_vLJy[3] = -0.803097279216327/refLength;
		_vLJz[3] = -1.28531461746939/refLength;
		_vLJm[3] = 15/refMass;
		_vLJeps[3] = 111.97954865097/refEnergy;
		_vLJsigma[3] = 3.62249459146099/refLength;
		_vdipoleX[0] = 0;
		_vdipoleY[0] = 0;
		_vdipoleZ[0] = 0;
		_vdipoleMyx[0] = 0;
		_vdipoleMyy[0] = 1;
		_vdipoleMyz[0] = 0;
		_vdipoleMyabs[0] = 3.44483983941682 * 85.105589165/refLength/sqrt(refLength*refEnergy);
		//Quadrupole
		_vquadrupoleX[0] = 0;
		_vquadrupoleY[0] = -0.803097279216327/refLength;
		_vquadrupoleZ[0] = 0;
		_vquadrupoleQx[0] = 0;
		_vquadrupoleQy[0] = 1;
		_vquadrupoleQz[0] = 0;
		_vquadrupoleQabs[0] = 62.198291205748674/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//Trägheitsmoment
		ixx = 92.317446203870375/(pow(refLength,2)*refMass);//\Sum_Einzelmassen mass_i * r_i^2
		iyy = 49.561009976414525/(pow(refLength,2)*refMass);//2. geht um x- und z-Achse
		izz = 42.756436227455858/(pow(refLength,2)*refMass);//3. und 4. gehen um alle 3 Achsen
		//allgemeine Parameter
		mass = 58/refMass;
	}
	else if (substance == "2cljdq"){//von Jadran geschickt 22.02.2012
		_numberLJCenters = 2;
		_numberCharges = 0;
		_numberQuadrupoles = 1;
		_numberDipoles = 1;
		_numberTersoff = 0;
		//Vektoren auf richtige Größen bringen
		resize();
		//LJ, charges etc. setzen
		_vLJx[0] = 0;
		_vLJy[0] = 0;
		_vLJz[0] = 0.5/refLength;
		_vLJm[0] = 1/refMass;
		_vLJeps[0] = 1/refEnergy;
		_vLJsigma[0] = 1/refLength;
		_vLJx[1] = 0;
		_vLJy[1] = 0;
		_vLJz[1] = -0.5/refLength;
		_vLJm[1] = 1/refMass;
		_vLJeps[1] = 1/refEnergy;
		_vLJsigma[1] = 1/refLength;
		//Dipol
		_vdipoleX[0] = 0;
		_vdipoleY[0] = 0;
		_vdipoleZ[0] = 0;
		_vdipoleMyx[0] = 0;
		_vdipoleMyy[0] = 0;
		_vdipoleMyz[0] = 1;
		_vdipoleMyabs[0] = 2.4494897427831779/refLength/sqrt(refLength*refEnergy);
		//Quadrupol
		_vquadrupoleX[0] = 0;
		_vquadrupoleY[0] = 0;
		_vquadrupoleZ[0] = 0;
		_vquadrupoleQx[0] = 0;
		_vquadrupoleQy[0] = 0;
		_vquadrupoleQz[0] = 1;
		_vquadrupoleQabs[0] = 2/(pow(refLength,2)*sqrt(refLength*refEnergy));
		//Trägheitsmoment
		ixx = 1/(pow(refLength,2)*refMass);
		iyy = 1/(pow(refLength,2)*refMass);
		izz = 0;
		//allgemeine Parameter
		mass = 2/refMass;
	}
	else
		cerr << "Not a valid component. Aborting.\n";
}

unsigned Component::size(){
	unsigned n = 5;
	for(int i = 0; i <= _numberLJCenters-1; i++){
		n += 6;/////////
	}
	for(int i = 0; i <= _numberDipoles-1; i++){
		n += 7;
	}
	for(int i = 0; i <= _numberQuadrupoles-1; i++){
		n += 7;
	}
	return n;
}

float Component::operator[](int index){
	vector<float> comp;
	comp.push_back(_numberLJCenters);
	comp.push_back(_numberCharges);
	comp.push_back(_numberDipoles);
	comp.push_back(_numberQuadrupoles);
	comp.push_back(_numberTersoff);
	for (int i = 0; i <= _numberLJCenters-1; i++){
		comp.push_back(_vLJx[i]);
		comp.push_back(_vLJy[i]);
		comp.push_back(_vLJz[i]);
		comp.push_back(_vLJm[i]);
		comp.push_back(_vLJeps[i]);
		comp.push_back(_vLJsigma[i]);
	}
	for (int i = 0; i <= _numberDipoles-1; i++){
		comp.push_back(_vdipoleX[i]);
		comp.push_back(_vdipoleY[i]);
		comp.push_back(_vdipoleZ[i]);
		comp.push_back(_vdipoleMyx[i]);
		comp.push_back(_vdipoleMyy[i]);
		comp.push_back(_vdipoleMyz[i]);
		comp.push_back(_vdipoleMyabs[i]);
	}
	for (int i = 0; i <= _numberQuadrupoles-1; i++){
		comp.push_back(_vquadrupoleX[i]);
		comp.push_back(_vquadrupoleY[i]);
		comp.push_back(_vquadrupoleZ[i]);
		comp.push_back(_vquadrupoleQx[i]);
		comp.push_back(_vquadrupoleQy[i]);
		comp.push_back(_vquadrupoleQz[i]);
		comp.push_back(_vquadrupoleQabs[i]);
	}
	return comp[index];
}

ostream& operator<< (ostream& os, Component& comp){
	os.precision(9);
	os << "component output:" << endl;
	os << comp._numberLJCenters << " " << comp._numberCharges << " " << comp._numberQuadrupoles << " "
		<< comp._numberDipoles << " " << comp._numberTersoff << endl;
	for(int i = 0; i <= comp._numberLJCenters-1; i++){
		os << comp._vLJx[i] << " " << comp._vLJy[i] << " " << comp._vLJz[i] << "\t"
			<< comp._vLJm[i] << " " << comp._vLJeps[i] << " " << comp._vLJsigma[i] << endl;
	}
	for(int i = 0; i <= comp._numberQuadrupoles-1; i++){
		os << comp._vquadrupoleX[i] << " " << comp._vquadrupoleY[i] << " " << comp._vquadrupoleZ[i] << "\t"
			<< comp._vquadrupoleQx[i] << " " << comp._vquadrupoleQy[i] << " " << comp._vquadrupoleQz[i] << " " << comp._vquadrupoleQabs[i] << "\n";
	}
	float dummy1 = 0, dummy2 = 0, dummy3 = 0;
	os << dummy1 << " " << dummy2 << " " << dummy3 << endl;

	return os;
}

vector<float> Component::gComponent() const{
	vector<float> comp;
	int size = 5+_numberLJCenters*(6+2)+_numberDipoles*7+_numberQuadrupoles*7+3;
	comp.resize(size);

	comp[0] = _numberLJCenters;
	comp[1] = _numberCharges;
	comp[2] = _numberDipoles;
	comp[3] = _numberQuadrupoles;
	comp[4] = _numberTersoff;
	int i = 5;
	int n = 0;
	while (n <= _numberLJCenters-1){
		comp[i] = _vLJx[n];
		comp[i+1] = _vLJy[n];
		comp[i+2] = _vLJz[n];
		comp[i+3] = _vLJm[n];
		comp[i+4] = _vLJeps[n];
		comp[i+5] = _vLJsigma[n];
		comp[i+6] = 0;
		comp[i+7] = 0;
		i = i+8;
		n++;
	}
	n = 0;
	while (n <= _numberDipoles-1){
		comp[i] = _vdipoleX[n];
		comp[i+1] = _vdipoleY[n];
		comp[i+2] = _vdipoleZ[n];
		comp[i+3] = _vdipoleMyx[n];
		comp[i+4] = _vdipoleMyy[n];
		comp[i+5] = _vdipoleMyz[n];
		comp[i+6] = _vdipoleMyabs[n];
		i = i+7;
		n++;
	}
	n = 0;
	while (n <= _numberQuadrupoles-1){
		comp[i] = _vquadrupoleX[n];
		comp[i+1] = _vquadrupoleY[n];
		comp[i+2] = _vquadrupoleZ[n];
		comp[i+3] = _vquadrupoleQx[n];
		comp[i+4] = _vquadrupoleQy[n];
		comp[i+5] = _vquadrupoleQz[n];
		comp[i+6] = _vquadrupoleQabs[n];
		i = i+7;
		n++;
	}
	comp[i] = 0;
	comp[i+1] = 0;
	comp[i+2] = 0;
	i = i+3;
	return comp;
}

double Component::gMass() const{
	return mass;
}

double Component::gIxx() const{
	return ixx;
}

double Component::gIyy() const{
	return iyy;
}

double Component::gIzz() const{
	return izz;
}

vector<float> Component::gvLJx() const{
	return _vLJx;
}

vector<float> Component::gvLJy() const{
	return _vLJy;
}
vector<float> Component::gvLJz() const{
	return _vLJz;
}

int Component::gnumberLJCenters() const{
	return _numberLJCenters;
}

int Component::resize(){
	//LJ
	_vLJx.resize(_numberLJCenters);
	_vLJy.resize(_numberLJCenters);
	_vLJz.resize(_numberLJCenters);
	_vLJm.resize(_numberLJCenters);
	_vLJeps.resize(_numberLJCenters);
	_vLJsigma.resize(_numberLJCenters);
	//Dipol
	_vdipoleX.resize(_numberDipoles);
	_vdipoleY.resize(_numberDipoles);
	_vdipoleZ.resize(_numberDipoles);
	_vdipoleMyx.resize(_numberDipoles);
	_vdipoleMyy.resize(_numberDipoles);
	_vdipoleMyz.resize(_numberDipoles);
	_vdipoleMyabs.resize(_numberDipoles);
	//Quadrupol
	_vquadrupoleX.resize(_numberQuadrupoles);
	_vquadrupoleY.resize(_numberQuadrupoles);
	_vquadrupoleZ.resize(_numberQuadrupoles);
	_vquadrupoleQx.resize(_numberQuadrupoles);
	_vquadrupoleQy.resize(_numberQuadrupoles);
	_vquadrupoleQz.resize(_numberQuadrupoles);
	_vquadrupoleQabs.resize(_numberQuadrupoles);
	return 0;
}
