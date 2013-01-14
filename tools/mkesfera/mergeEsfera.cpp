
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <math.h>
#include <string.h>
#include <vector>
#include "limits.h"
#include "stdlib.h"
#include "Random.h"

const unsigned PRECISION = 10;
const double N_IDEAL = 150000.0;

using namespace std;

int main(int argc, char** argv){
  
  if(argc != 4 ){
   cout << "\n\nProgramm zum Zusammenführen einer esfera-restart Datei\nund der Wandpartikel aus mkSD *.template.inp-Datei\n\nGebrauchsanweisung: <Prefix Quelle mkSD-file> < Prefix Quelle esfera-Datei> <y_cut>\n\n";
   return 0;
  }
  
  string quelleSDPrefix, quelleSDName, quelleEsfPrefix, quelleEsfName, zielName, xyzName;
  ifstream quelleSD, quelleEsf;
  ofstream ziel, xyz;
  double y_cut;
  
  
  quelleSDPrefix = argv[1];
  quelleSDName = quelleSDPrefix + ".template.inp";
  zielName = quelleSDPrefix + ".inp";
  xyzName = quelleSDPrefix + ".xyz";
  quelleEsfPrefix = argv[2];
  quelleEsfName = quelleEsfPrefix + ".restart.xdr";
  y_cut = atof(argv[3]);
   
   
  
  unsigned int NWall, NFluid, NTotal;
  double TemperaturSD;//, TemperaturEsf;
  double boxSD[3];
  double berthelotXi;
  vector<double> wallX, wallY, wallZ, wallU, wallV, wallW;
  double hWall = 0.0;
  double hFluid = 0.0;
  vector<double> fluidX, fluidY, fluidZ, fluidU, fluidV, fluidW;
  

  
  
   
  
//********************************************************************** SD-file ***********************************************************

  quelleSD.open(quelleSDName.c_str(),ios::in);
  if(quelleSD == NULL){
      cout << "Could not open: " << quelleSDName << endl;
      return 0;
    }   
  
  //! Auslesen header SD-file
  bool header = true;
  while(header){
   
   char c;
   quelleSD >> c;
   if(c == '#') {
   quelleSD.ignore( INT_MAX,'\n' );
   continue;
   }
   quelleSD.putback(c);

   string token;
   //token.clear();
   quelleSD >> token;
   if(token == "ThT"){
     int kaese;
     quelleSD >> kaese >> TemperaturSD;
   }
   else if(token == "L"){
     quelleSD >> boxSD[0] >> boxSD[1] >> boxSD[2];
   }
   else if(token == "C"){
     int numComp;
     double eta;
     quelleSD >> numComp;
     for(unsigned short i = 0; i<=4; i++){
       string egal;
       getline(quelleSD, egal, '\n');// zur richtigen Leseposition: 4 Zeilen weiter
       egal.clear();
     }
     quelleSD >> berthelotXi >> eta;
   }
   else if(token == "M"){
     string molFormat;
    quelleSD >> molFormat; // zur richtigen Leseposition
    header = false; 
   }
   else {
     token.clear();
     /*
      * cout << "Unzulaessiger token: " << token <<"\nBisher ausgelesen:\nTemperatur SD = " << TemperaturSD << "\nBoxlänge SD: " << boxSD[0] << " "
	  << boxSD[1] << " " << boxSD[2] << "\nBerthelot xi = " << berthelotXi << "\n\n";
      */
   }
  }
  cout << "\nTemperatur SD = " << TemperaturSD << "\nBoxlänge SD: " << boxSD[0] << " "
	  << boxSD[1] << " " << boxSD[2] << "\nBerthelot xi = " << berthelotXi << "\n\n";
  
  //! Auslesen body SD
  while(!quelleSD.eof()){
    int id, cid;
    double x,y,z,u,v,w,q0,q1,q2,q3,bla1,bla2,bla3;
    quelleSD >> id >> cid >> x >> y >> z >> u >> v >> w >> q0 >> q1 >> q2 >> q3 >> bla1 >> bla2 >> bla3;
    if(cid > 1 ){
      wallX.push_back(x);
      wallY.push_back(y);
      wallZ.push_back(z);
      wallU.push_back(u);
      wallV.push_back(v);
      wallW.push_back(w);
      if(y>hWall) hWall = y; 
    }
  }
  
  double ws = wallX.size();
  if(wallX[ws-1] == wallX[ws-2] && wallY[ws-1] == wallY[ws-2] && wallZ[ws-1] == wallZ[ws-2]) ws = ws -1;    
  NWall = ws;
  quelleSD.close();
  cout << "Body SD-Datei ausgelesen\nN_Wall = " << NWall << "\n";
  
//************************************************************** esfera-file ***************************************************************
  
  quelleEsf.open(quelleEsfName.c_str(),ios::in);
  if(quelleEsf == NULL){
      cout << "Could not open: " << quelleEsfName << endl;
      return 0;
    }
    
    //! header esfera-datei
    header = true;
    while(header){
      string token, token0;
      getline(quelleEsf, token, '\n');
      token0 = token[1];
      if(token0 == "M"){
	header = false;
      }
      token.clear();
    }
    
    cout << "Header esfera-Datei abgearbeitet\n";
    while(!quelleEsf.eof()){
      int id, cid;
      double x,y,z,u,v,w,q0,q1,q2,q3,bla1,bla2,bla3;
      quelleEsf >> id >> cid >> x >> y >> z >> u >> v >> w >> q0 >> q1 >> q2 >> q3 >> bla1 >> bla2 >> bla3;
      if(y > y_cut){
	fluidX.push_back(x+1.0);
	fluidY.push_back(y-y_cut + hWall + 1.2);
	fluidZ.push_back(z+1.0);
	fluidU.push_back(u);
	fluidV.push_back(v);
	fluidW.push_back(w);
	if(y > hFluid) hFluid = y;
      }      
    }
    hFluid = hFluid - y_cut +hWall + 2.2;
    double fs = fluidX.size();
    if(fluidX[fs-1] == fluidX[fs-2] && fluidY[fs-1] == fluidY[fs-2] && fluidZ[fs-1] == fluidZ[fs-2]) fs = fs -1;
    NFluid = fs;
    quelleEsf.close();
    cout << "Body esfera ausgelesen\nN_Fluid_total = "<< NFluid << "\n";
        
    //! fill probability in order to obtain the precise Number of 150 000 fluid particles to be filled
    vector<bool> filled(NFluid);
    double nFilled  = double(NFluid);
    double totalNSlots = double(NFluid);
    double pswap;
    bool tswap;
    Random* rdm = new Random(); 
    rdm->init(12);
    for(unsigned i = 0; i <filled.size(); i++){
     filled[i] = true; 
    }
     for(unsigned short m = 0; m < PRECISION; m++){
      tswap = (nFilled < N_IDEAL);
      pswap = (N_IDEAL - nFilled) / ((tswap ? totalNSlots : 0.0) -nFilled);
      for(unsigned i = 0; i< NFluid; i++){
	if(pswap >= rdm->rnd()){
	 if(filled[i]) nFilled --;
	 filled[i] = tswap;
	 if(tswap) nFilled++;
	}
      }
       
     }
    cout << "Nach Entfernen überzaehliger Fluidteiclchen\nN_Fluid = "<< nFilled << "\n";
    
//***************************************************************** - neue input Datei *****************************************************
  
    NTotal = NWall + nFilled;

  ziel.open(zielName.c_str(), ios::trunc);
  if(ziel == NULL){
      cout << "Could not generate a merged output file.\n";
      return 0;
    }
    
  xyz.open(xyzName.c_str(), ios::trunc);
  if(xyz == NULL){
      cout << "Could not generate a merged xyz-output file.\n";
      return 0;
    }
  xyz << NTotal << "\n comment\n";
  
  //! Schreibe header
  ziel << "mardyn trunk  20130111\n#mardyn input file, ls1 project\n#generated by merging a equilibrated fluid drop (esfera) with a plane wall\n";
  ziel << "t\t0.0\n";
  ziel << "CT 2 1\nThT 1 " << TemperaturSD << "\nCT 1 2\nThT 2 "<< TemperaturSD << "\n";
  ziel << "L\t" << boxSD[0] << " " <<  hFluid << " " << boxSD[2] << "\n";
  ziel << "C\t2 \n1 0 0 0 0\n0.0 0.0 0.0 \t1 1 1 2.5 1 \t0.0 0.0 0.0\n1 0 0 0 0\n0.0 0.0 0.0 \t1.59104 100 1 2.5 1 \t0.0 0.0 0.0\n ";
  ziel << berthelotXi << " 1 \n1e+10\n";
  ziel << "N\t" << NTotal << "\nM\tICRVQD\n\n\n"; 
  
  //! Schreibe body
  unsigned id, cid;
  cid = 1;
  id = 1;
  for(unsigned int i = 0; i < NFluid; i++){
    if(filled[i]){
      ziel << id++ << " " << cid << " " << fluidX[i] << " " << fluidY[i] << " " << fluidZ[i] << " " << fluidU[i] << " " << fluidV[i] << " " << fluidW[i] << "\t1.0 0.0 0.0 0.0        0 0 0\n";
      xyz << "O " << fluidX[i] << " " << fluidY[i] << " " << fluidZ[i] << "\n";
    }
  }
  cid+=1;
  for(unsigned int i = 0; i < NWall; i++){
   ziel << id++ << " " << cid << " " << wallX[i] << " " << wallY[i] << " " << wallZ[i] << " " << wallU[i] << " " << wallV[i] << " " << wallW[i] << "\t1.0 0.0 0.0 0.0        0 0 0\n";
   xyz << "H " << wallX[i] << " " << wallY[i] << " " << wallZ[i] << "\n";
  }
  
  ziel.close();
  
  return 0;
}