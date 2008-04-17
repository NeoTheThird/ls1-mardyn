#include "Domain.h"

#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include <cmath>

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
using namespace std;

utils::Log Domain::_log("Domain");

// used by Domain::init_Corr() -----------
//! @todo WHAT DOES THIS METHOD DO?
double TICCu(int n,double rc,double sigma2);
//! @todo WHAT DOES THIS METHOD DO?
double TICSu(int n,double rc,double sigma2,double tau);
//! @todo WHAT DOES THIS METHOD DO?
double TISSu(int n,double rc,double sigma2,double tau1,double tau2);
//! @todo WHAT DOES THIS METHOD DO?
double TICCv(int n,double rc,double sigma2);
//! @todo WHAT DOES THIS METHOD DO?
double TICSv(int n,double rc,double sigma2,double tau);
//! @todo WHAT DOES THIS METHOD DO?
double TISSv(int n,double rc,double sigma2,double tau1,double tau2);
//----------------------------------------


Domain::Domain(int rank){
  this->_localRank = rank;
  this->_localUpot = 0;
  this->_localVirial = 0;   
  this->_globalUpot = 0;
  this->_globalVirial = 0; 
  this->_globalRho = 0;
  this->_universalThermostatID = map<int, int>();
  this->_localThermostatN = map<int, unsigned long>();
  this->_localThermostatN[-1] = 0;
  this->_localThermostatN[0] = 0;
  this->_universalThermostatN = map<int, unsigned long>();
  this->_universalThermostatN[-1] = 0;
  this->_universalThermostatN[0] = 0;
  this->_localRotationalDOF = map<int, unsigned long>();
  this->_localRotationalDOF[-1] = 0;
  this->_localRotationalDOF[0] = 0;
  this->_universalRotationalDOF = map<int, unsigned long>();
  this->_universalRotationalDOF[-1] = 0;
  this->_universalRotationalDOF[0] = 0;
  this->_globalLength[0] = 0;
  this->_globalLength[1] = 0;
  this->_globalLength[2] = 0;
  this->_universalBTrans = map<int, double>();
  this->_universalBTrans[0] = 1.0;
  this->_universalBRot = map<int, double>();
  this->_universalBRot[0] = 1.0;
  this->_universalTargetTemperature = map<int, double>();
  this->_universalTargetTemperature[0] = 1.0;
  this->_globalTemperatureMap = map<int, double>();
  this->_globalTemperatureMap[0] = 1.0;
  this->_local2KETrans[0] = 0.0;
  this->_local2KERot[0] = 0.0; 
  this->_currentTime = 0.0;
#ifdef COMPLEX_POTENTIAL_SET
  this->_universalConstantAccelerationTimesteps = 30;
  if(!rank)
    for(unsigned d=0; d < 3; d++)
      this->_globalVelocitySum[d] = map<unsigned, long double>();
#endif
  this->_universalComponentwiseThermostat = false;
#ifdef COMPLEX_POTENTIAL_SET
  this->_universalUndirectedThermostat = map<int, bool>();
  for(int d = 0; d <3; d++)
  {
    this->_universalThermostatDirectedVelocity[d] = map<int, double>();
    this->_localThermostatDirectedVelocity[d] = map<int, double>();
  }
#endif
}

void Domain::setLocalUpot(double Upot) {_localUpot = Upot;}

double Domain::getLocalUpot() const {return _localUpot; }

void Domain::setLocalVirial(double Virial) {_localVirial = Virial;}

double Domain::getLocalVirial() const {return _localVirial; }

double Domain::getGlobalBetaTrans() 
{
   // if(!this->_universalComponentwiseThermostat) return this->_universalBTrans[0];
   // else return this->_universalBTrans[1];
   return this->_universalBTrans[0];
}
double Domain::getGlobalBetaTrans(int thermostat) 
{
   // if(thermostat < 0) return 1.0;
   // if(!this->_universalComponentwiseThermostat) return this->_universalBTrans[0];
   // else return this->_universalBTrans[thermostat];
   return this->_universalBTrans[thermostat];
}

double Domain::getGlobalBetaRot() 
{
   // if(!this->_universalComponentwiseThermostat) return this->_universalBRot[0];
   // else return this->_universalBRot[1];
   return this->_universalBRot[0];
}
double Domain::getGlobalBetaRot(int thermostat)
{
   // if(thermostat < 0) return 1.0;
   // if(!this->_universalComponentwiseThermostat) return this->_universalBRot[0];
   // else return this->_universalBRot[thermostat];
   return this->_universalBRot[thermostat];
}

double Domain::getGlobalPressure()
{
   double globalTemperature = this->_globalTemperatureMap[0];
   return globalTemperature*_globalRho+_globalRho*getAverageGlobalVirial()/3.;
}

double Domain::getAverageGlobalVirial() const { return _globalVirial/_globalNumMolecules; }

double Domain::getAverageGlobalUpot() const { return _globalUpot/_globalNumMolecules; }

void Domain::setLocalSummv2(double summv2, int thermostat)
{
#ifndef NDEBUG
   if(this->_localRank == 0)
   {
      cout << "      * local data from proc. 0 thermostat " << thermostat
           << ":  mvv = " << summv2 << "\n";
   }
#endif
   this->_local2KETrans[thermostat] = summv2;
}
void Domain::setLocalSummv2(double summv2) { setLocalSummv2(summv2, 0); }
    
void Domain::setLocalSumIw2(double sumIw2, int thermostat)
{
   _local2KERot[thermostat] = sumIw2;
} 
void Domain::setLocalSumIw2(double sumIw2) { setLocalSumIw2(sumIw2, 0); }

double Domain::getCurrentTime(){ return _currentTime;}

void Domain::advanceTime(double timestep){ _currentTime += timestep;}

vector<Component>& Domain::getComponents(){
  return _components; 
}

Comp2Param& Domain::getComp2Params(){
  return _comp2params; 
}

double Domain::getGlobalLength(int index) const {
  return _globalLength[index];
}

#ifdef COMPLEX_POTENTIAL_SET
void Domain::calculateGlobalValues( parallel::DomainDecompBase* domainDecomp,
                                    datastructures::ParticleContainer<Molecule>*
                                       particleContainer,
                                    bool collectThermostatVelocities )
#else
void Domain::calculateGlobalValues( parallel::DomainDecompBase* domainDecomp,
                                    datastructures::ParticleContainer<Molecule>*
                                       particleContainer )
#endif
{
   double Upot = _localUpot;
   double Virial = _localVirial;
   domainDecomp->reducevalues(&Upot, &Virial);

   // Process 0 has to add the charge and dipole correction:
   // m_UpotCorr and m_VirialCorr already contain constant (internal) dipole correction
   _globalUpot = Upot + _UpotCorr;
   _globalVirial = Virial + _VirialCorr;

   /*
    * thermostat ID 0 represents the entire system
    */
   map<int, unsigned long>::iterator thermit;
   if(this->_universalComponentwiseThermostat)
   {
#ifndef NDEBUG
      if(this->_localRank == 0) cout << "      * applying a componentwise thermostat\n";
#endif
      this->_localThermostatN[0] = 0;
      this->_localRotationalDOF[0] = 0;
      this->_local2KETrans[0] = 0;
      this->_local2KERot[0] = 0;
      for(thermit = _localThermostatN.begin(); thermit != _localThermostatN.end(); thermit++)
      {
         if(thermit->first == 0) continue;
         this->_localThermostatN[0] += thermit->second;
         this->_localRotationalDOF[0] += this->_localRotationalDOF[thermit->first];
         this->_local2KETrans[0] += this->_local2KETrans[thermit->first];
         this->_local2KERot[0] += this->_local2KERot[thermit->first];
      }
   }

   for(thermit = _universalThermostatN.begin(); thermit != _universalThermostatN.end(); thermit++)
   {
      // number of molecules on the local process. After the reduce operation
      // num_molecules will contain the global number of molecules
      unsigned long numMolecules = this->_localThermostatN[thermit->first];
#ifndef NDEBUG
      cout << "    [[[ " << thermit->first << " ] N=" << numMolecules << " ] localRank=" << _localRank << " ]   ";
#endif
      double summv2 = this->_local2KETrans[thermit->first];
      double sumIw2 = this->_local2KERot[thermit->first];
      unsigned long rotDOF = this->_localRotationalDOF[thermit->first];
      domainDecomp->reducevalues(&summv2, &sumIw2, &numMolecules, &rotDOF);
      this->_universalThermostatN[thermit->first] = numMolecules;
      this->_universalRotationalDOF[thermit->first] = rotDOF;
      if(numMolecules > 0)
         this->_globalTemperatureMap[thermit->first] =
            (summv2+sumIw2) / (3.0*numMolecules + rotDOF);
      else
         _globalTemperatureMap[thermit->first] = _universalTargetTemperature[thermit->first];
      double Ti = this->_universalTargetTemperature[thermit->first];
      if(Ti > 0.0)
      {
         this->_universalBTrans[thermit->first] = sqrt( sqrt(3.0*numMolecules*Ti / summv2) );
         if(sumIw2 == 0.0) this->_universalBRot[thermit->first] = 1.0;
         else this->_universalBRot[thermit->first] = sqrt( sqrt(rotDOF*Ti / sumIw2) );
      }
      else
      {
         this->_universalBTrans[thermit->first] = 1.0;
         this->_universalBRot[thermit->first] = 1.0;
      }

#ifdef COMPLEX_POTENTIAL_SET
      if(collectThermostatVelocities && this->_universalUndirectedThermostat[thermit->first])
      {
         double sigv[3];
         for(int d=0; d < 3; d++)
            sigv[d] = this->_localThermostatDirectedVelocity[d][thermit->first];
         domainDecomp->reducevalues(&sigv[0], (double*)0);
         domainDecomp->reducevalues(&sigv[1], &sigv[2]);
         for(int d=0; d < 3; d++)
         {
            this->_localThermostatDirectedVelocity[d][thermit->first] = 0.0;
            if(numMolecules > 0)
               _universalThermostatDirectedVelocity[d][thermit->first] = sigv[d] / numMolecules;
            else _universalThermostatDirectedVelocity[d][thermit->first] = 0.0;
         }

#ifndef NDEBUG
         cout << "\n         * thermostat " << thermit->first
              << " directed velocity: ("
              << this->_universalThermostatDirectedVelocity[0][thermit->first]
              << " / " << this->_universalThermostatDirectedVelocity[1][thermit->first]
              << " / " << this->_universalThermostatDirectedVelocity[2][thermit->first] << ")\n";
#endif
      }
#endif

#ifndef NDEBUG
      if(this->_localRank == 0)
      {
         cout << "      * Th" << thermit->first << " N=" << numMolecules
              << " DOF=" << rotDOF + 3.0*numMolecules << "\n"
              << "        Tcur=" << _globalTemperatureMap[thermit->first]
              << " Ttar=" << _universalTargetTemperature[thermit->first]
              << " bt=" << _universalBTrans[thermit->first]
              << " br=" << _universalBRot[thermit->first] << "\n";
      }
#endif
   }
}

#ifdef COMPLEX_POTENTIAL_SET
void Domain::calculateThermostatDirectedVelocity(datastructures::ParticleContainer<Molecule>* partCont)
{
   Molecule* tM;
   if(this->_universalComponentwiseThermostat)
   {
      for( map<int, bool>::iterator thit = _universalUndirectedThermostat.begin();
           thit != _universalUndirectedThermostat.end();
           thit ++ )
      {
         if(thit->second)
            for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][thit->first] = 0.0;
      }
      for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
      {
         int cid = tM->componentid();
         int thermostat = this->_universalThermostatID[cid];
         if(this->_universalUndirectedThermostat[thermostat])
         {
            for(int d=0; d < 3; d++)
               _localThermostatDirectedVelocity[d][thermostat] += tM->v(d);
         }
      }
   }
   else if(this->_universalUndirectedThermostat[0])
   {
      for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][0] = 0.0;
      for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
      {
         for(int d=0; d < 3; d++)
               _localThermostatDirectedVelocity[d][0] += tM->v(d);
      }
   }
}
#endif

void Domain::calculateVelocitySums(datastructures::ParticleContainer<Molecule>* partCont)
{
   Molecule* tM;
   if(this->_universalComponentwiseThermostat)
   {
      for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
      {
         int cid = tM->componentid();
         int thermostat = this->_universalThermostatID[cid];
         this->_localThermostatN[thermostat]++;
         this->_localRotationalDOF[thermostat] += _components[cid].rot_dof();
#ifdef COMPLEX_POTENTIAL_SET
         if(this->_universalUndirectedThermostat[thermostat])
         {
            tM->calculate_mv2_Iw2( this->_local2KETrans[thermostat],
                                   this->_local2KERot[thermostat],
                                   this->_universalThermostatDirectedVelocity[0][thermostat],
                                   this->_universalThermostatDirectedVelocity[1][thermostat],
                                   this->_universalThermostatDirectedVelocity[2][thermostat]  );
         }
         else
         {
#endif
            tM->calculate_mv2_Iw2(_local2KETrans[thermostat], _local2KERot[thermostat]);
#ifdef COMPLEX_POTENTIAL_SET
         }
#endif
      }
   }
   else
   {
      for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
      {
         this->_localThermostatN[0]++;
         this->_localRotationalDOF[0] += _components[ tM->componentid() ].rot_dof();
#ifdef COMPLEX_POTENTIAL_SET
         if(this->_universalUndirectedThermostat[0])
         {
            tM->calculate_mv2_Iw2( this->_local2KETrans[0],
                                   this->_local2KERot[0],
                                   this->_universalThermostatDirectedVelocity[0][0],
                                   this->_universalThermostatDirectedVelocity[1][0],
                                   this->_universalThermostatDirectedVelocity[2][0]  );
         }
         else
         {
#endif
            tM->calculate_mv2_Iw2(_local2KETrans[0], _local2KERot[0]);
#ifdef COMPLEX_POTENTIAL_SET
         }
#endif
      }
#ifndef NDEBUG
      if(this->_localRank == 0)
      {
         cout << "      * local data from proc. 0: N = " << this->_localThermostatN[0]
              << "   rotDOF = " << this->_localRotationalDOF[0] << "   mvv = "
              << _local2KETrans[0] << " Iww = " << _local2KERot[0] << "\n";
      }
#endif
   }
}

void Domain::setPhaseSpaceFile(string filename){
  _phaseSpaceFileStream.open(filename.c_str());
}  

#ifdef COMPLEX_POTENTIAL_SET
void Domain::readPhaseSpaceHeader(double timestep, double cutoffTersoff)
#else
void Domain::readPhaseSpaceHeader(double timestep)
#endif
{
  string token;
  _phaseSpaceFileStream >> token;
  _inpversion=0;
  if((token != "MDProject") && (token != "MOLDY") && (token != "ls1r1") && (token != "mrdyn") && (token != "mardyn"))
  {
    if(_localRank == 0) cerr << "Input: NOT A MOLDY INPUT! (starts with " << token << ")" << endl;
    exit(1);
  }
  else
  {
    if((token != "mardyn") && (_localRank == 0))
       cerr << "Warning: phase space file should begin with 'mardyn' instead of '" << token << "'.\n";

    unsigned REQUIRED_INPUT_VERSION = 20070702;
#ifdef COMPLEX_POTENTIAL_SET
    REQUIRED_INPUT_VERSION = 20070702;
#endif

    _phaseSpaceFileStream >> _inpversion;
    if(_inpversion < REQUIRED_INPUT_VERSION)
    {
      if(_localRank == 0) cerr << "Input: OLD VERSION (" << _inpversion << ")" << endl;
      exit(1);
    }
  }

  char c;
  double x,y,z;
  unsigned int numcomponents=0;
  unsigned int j;
  double m,sigma,eps;
  double xi,eta;
  unsigned long i;

  // When the last header element is reached, "header" is set to false
  bool header = true;
  
  while(header) {
    _phaseSpaceFileStream >> c;
    if(c=='#') {
      _phaseSpaceFileStream.ignore(INT_MAX,'\n');
      continue;
    }
    else {
      _phaseSpaceFileStream.putback(c);
    }
    token.clear();
    _phaseSpaceFileStream >> token;
    if((token == "currentTime") || (token == "t"))
    {
      _phaseSpaceFileStream >> _currentTime;
    }
    else if((token == "Temperature") || (token == "T"))
    {
      _phaseSpaceFileStream >> _universalTargetTemperature[0];
      this->_universalComponentwiseThermostat = false;
#ifdef COMPLEX_POTENTIAL_SET
      this->_universalUndirectedThermostat[0] = false;
#endif
    }
    else if((token == "ThermostatTemperature") || (token == "ThT") || (token == "h"))
    {
       int i;
       _phaseSpaceFileStream >> i;
       _phaseSpaceFileStream >> _universalTargetTemperature[i];
#ifdef COMPLEX_POTENTIAL_SET
       this->_universalUndirectedThermostat[i] = false;
#endif
       if(i == 0)
       {
          this->_universalComponentwiseThermostat = false;
       }
       else
       {
          if(!this->_universalComponentwiseThermostat)
          {
             this->_universalComponentwiseThermostat = true;
             this->_universalTargetTemperature.erase(0);
#ifdef COMPLEX_POTENTIAL_SET
             this->_universalUndirectedThermostat.erase(0);
             for(int d=0; d < 3; d++) this->_universalThermostatDirectedVelocity[d].erase(0);
#endif
             for( vector<Component>::iterator tc = this->_components.begin();
                  tc != this->_components.end();
                  tc ++ )
                if(!(this->_universalThermostatID[ tc->ID() ] > 0))
                   this->_universalThermostatID[ tc->ID() ] = -1;
          }
       }
    }
    else if((token == "ComponentThermostat") || (token == "CT") || (token == "o"))
    {
       if(!this->_universalComponentwiseThermostat)
       {
          this->_universalComponentwiseThermostat = true;
          this->_universalTargetTemperature.erase(0);
          for( vector<Component>::iterator tc = this->_components.begin();
               tc != this->_components.end();
               tc ++ )
             if(!(this->_universalThermostatID[ tc->ID() ] > 0))
                this->_universalThermostatID[ tc->ID() ] = -1;
       }
       unsigned cid;
       _phaseSpaceFileStream >> cid;
       cid--;
       unsigned th;
       _phaseSpaceFileStream >> th;
       this->_universalThermostatID[cid] = th;
       this->_universalThermostatN[th] = 0;
    }
#ifdef COMPLEX_POTENTIAL_SET
    else if((token == "Undirected") || (token == "U"))
    {
       int tst;
       _phaseSpaceFileStream >> tst;
       this->_universalUndirectedThermostat[tst] = true;
       for(int d=0; d < 3; d++)
       {
          this->_universalThermostatDirectedVelocity[d][tst] = 0.0;
          this->_localThermostatDirectedVelocity[d][tst] = 0.0;
       }
    }
#endif
    else if((token == "Length") || (token == "L"))
    {
      _phaseSpaceFileStream >> _globalLength[0] >> _globalLength[1] >> _globalLength[2];
      for(int i=0; i < 3; i++) _globalLength[i] *= 1.0000003;   // to avoid rounding problems with coordinates
    }
    else if((token == "NumberOfComponents") || (token == "C"))
    {
      _phaseSpaceFileStream >> numcomponents;
      _components.resize(numcomponents);
      for(i=0;i<numcomponents;++i)
      {
        _components[i].setID(i);
        unsigned int numljcenters=0;
        unsigned int numcharges=0;
        unsigned int numquadrupoles=0;
        _phaseSpaceFileStream >> numljcenters >> numcharges >> numquadrupoles;
#ifdef COMPLEX_POTENTIAL_SET
        unsigned int numdipoles = 0;
        unsigned int numtersoff = 0;
        _phaseSpaceFileStream >> numdipoles >> numtersoff;
#ifndef NDEBUG
        cout << "component with " << numljcenters << " " << numcharges << " " << numquadrupoles << " "
             << numdipoles << " " << numtersoff << "\n";
#endif
#endif
        for(j=0;j<numljcenters;++j)
        {
#ifndef NDEBUG
          if(!_localRank) cout << "reading LJ data\n";
#endif
          _phaseSpaceFileStream >> x >> y >> z >> m >> eps >> sigma;
          _components[i].addLJcenter(x,y,z,m,eps,sigma);
        }
        for(j = 0; j < numcharges; j++)
        {
#ifndef NDEBUG
          if(!_localRank) cout << "reading partial charge data\n";
#endif
          double q;
          _phaseSpaceFileStream >> x >> y >> z >> m >> q;
          _components[i].addCharge(x, y, z, m, q);
        }
        for(j=0;j<numquadrupoles;++j)
        {
#ifndef NDEBUG
          if(!_localRank) cout << "reading Q data\n";
#endif
          double eQx,eQy,eQz,absQ;
          _phaseSpaceFileStream >> x >> y >> z >> eQx >> eQy >> eQz >> absQ;
          _components[i].addQuadrupole(x,y,z,eQx,eQy,eQz,absQ);
        }
#ifdef COMPLEX_POTENTIAL_SET
        for(j=0;j<numdipoles;++j)
        {
#ifndef NDEBUG
          if(!_localRank) cout << "reading D data\n";
#endif
          double eMyx,eMyy,eMyz,absMy;
          _phaseSpaceFileStream >> x >> y >> z >> eMyx >> eMyy >> eMyz >> absMy;
          _components[i].addDipole(x,y,z,eMyx,eMyy,eMyz,absMy);
        }
        for(j = 0; j < numtersoff; j++)
        {
#ifndef NDEBUG
          if(!_localRank) cout << "reading Tersoff data\n";
#endif
          double x, y, z, m, A, B, lambda, mu, R, S, c, d, h, n, beta;
          _phaseSpaceFileStream >> x >> y >> z;
          _phaseSpaceFileStream >> m >> A >> B;
          _phaseSpaceFileStream >> lambda >> mu >> R >> S;
          /*
           * note that S = rcT is strongly recommended
           * and rcT > S would lead to WRONG results
           */
          if(S > cutoffTersoff)
          {
             if(!_localRank)
                cout << "severe error:   S = " << S << "  >  rcT = " << cutoffTersoff << "\n";
             exit(1);
          }
          else if(2.0*S < cutoffTersoff)
             if(!_localRank) 
                cout << "warning:   S = " << S << ",   rcT = " << cutoffTersoff << "\n";
          _phaseSpaceFileStream >> c >> d >> h >> n >> beta;
          _components[i].addTersoff(
            x, y, z,
            m, A, B,
            lambda, mu, R, S,
            c, d, h, n, beta
          );
#ifndef NDEBUG
          if(!_localRank) cout << x << " " << y << " " << z << " " << m << " " << A << " " << B << " " << lambda << " " << mu << " " << R << " " << S << " " << c << " " << d << " " << h << " " << n << " " << beta << "\n";
#endif
        }
#endif
        double IDummy1,IDummy2,IDummy3;
        _phaseSpaceFileStream >> IDummy1 >> IDummy2 >> IDummy3;
        if(IDummy1>0.) _components[i].setI11(IDummy1);
        if(IDummy2>0.) _components[i].setI22(IDummy2);
        if(IDummy3>0.) _components[i].setI33(IDummy3);
      }
      _mixcoeff.clear();
      for(i=0;i<numcomponents-1;++i)
      {
        for(j=i+1;j<numcomponents;++j)
        {
          _phaseSpaceFileStream >> xi >> eta;
          _mixcoeff.push_back(xi);
          _mixcoeff.push_back(eta);
        }
      }
      _phaseSpaceFileStream >> _epsilonRF;
      #ifndef NDEBUG
          cout << "eRF:\t" << _epsilonRF << "\n";
      #endif
    }
    else if((token == "NumberOfMolecules") || (token == "N"))
    {
      _phaseSpaceFileStream >> _globalNumMolecules;
      header = false;
    }
#ifdef COMPLEX_POTENTIAL_SET
    else if((token == "AssignCoset") || (token == "S"))
    {
      unsigned cid, cosetid;
      _phaseSpaceFileStream >> cid >> cosetid;
      cid--;
      this->assignCoset(cid, cosetid);
    }
    else if((token == "Accelerate") || (token == "A"))
    {
       unsigned cosetid;
       _phaseSpaceFileStream >> cosetid;
       double v[3];
       for(unsigned d = 0; d < 3; d++) _phaseSpaceFileStream >> v[d];
       double tau;
       _phaseSpaceFileStream >> tau;
       double ainit[3];
       for(unsigned d = 0; d < 3; d++) _phaseSpaceFileStream >> ainit[d];
       this->specifyComponentSet(cosetid, v, tau, ainit, timestep);
    }
#endif
  }
}

void Domain::writeCheckpoint(string filename, datastructures::ParticleContainer<Molecule>* particleContainer,
                             parallel::DomainDecompBase* domainDecomp, double dt)
{
  if(!this->_localRank)
  {
    ofstream checkpointfilestream(filename.c_str());
    checkpointfilestream << "mardyn\t20070717"<< endl;
    //! @todo use real current time
    checkpointfilestream << " t\t"  << this->_currentTime << "\n";
    checkpointfilestream << " dt\t" << dt << "\n";
    checkpointfilestream << " L\t" << _globalLength[0] << " " << _globalLength[1] << " " << _globalLength[2] << endl;
    checkpointfilestream << "# rho\t" << this->_globalRho << "\n";
    checkpointfilestream << "# rc\t" << particleContainer->getCutoff() << "\n";
#ifdef COMPLEX_POTENTIAL_SET
    checkpointfilestream << "# rcT\t" << particleContainer->getTersoffCutoff() << "\n";
#endif
    checkpointfilestream << "C\t" << _components.size() << endl;
    for(vector<Component>::const_iterator pos=_components.begin();pos!=_components.end();++pos){
      pos->write(checkpointfilestream);
    }
    unsigned int numperline=_components.size();
    unsigned int iout=0;
    for(vector<double>::const_iterator pos=_mixcoeff.begin();pos!=_mixcoeff.end();++pos){
      checkpointfilestream << *pos;
      iout++;
      // 2 parameters (xi and eta)
      if(iout/2>=numperline) {
        checkpointfilestream << endl;
        iout=0;
        --numperline;
      }
      else if(!(iout%2)) {
        checkpointfilestream << "\t";
      }
      else {
        checkpointfilestream << " ";
      }
    }
    checkpointfilestream << _epsilonRF << endl;
#ifdef COMPLEX_POTENTIAL_SET
    for( map<int, bool>::iterator uutit = this->_universalUndirectedThermostat.begin();
         uutit != this->_universalUndirectedThermostat.end();
         uutit++ )
    {
      if(uutit->second) checkpointfilestream << " U\t" << uutit->first << "\n";
    }
    for( map<unsigned, unsigned>::const_iterator uCSIDit = this->_universalComponentSetID.begin();
         uCSIDit != this->_universalComponentSetID.end();
         uCSIDit++ )
    { 
      if(uCSIDit->first > 100) continue;
      checkpointfilestream << " S\t" << 1+uCSIDit->first << "\t" << uCSIDit->second << "\n";
    }
    for( map<unsigned, double>::const_iterator gTit = _universalTau.begin();
         gTit != this->_universalTau.end();
         gTit++ )
    {
      unsigned cosetid = gTit->first;
      checkpointfilestream << " A\t" << cosetid << "\t"
                           << this->_globalTargetVelocity[0][cosetid] << " "
                           << this->_globalTargetVelocity[1][cosetid] << " "
                           << this->_globalTargetVelocity[2][cosetid] << "\t"
                           << gTit->second << "\t"
                           << this->_universalAdditionalAcceleration[0][cosetid] << " "
                           << this->_universalAdditionalAcceleration[1][cosetid] << " "
                           << this->_universalAdditionalAcceleration[2][cosetid] << "\n";
    }
#endif
    if(this->_universalComponentwiseThermostat)
    {
       for( map<int, int>::iterator thermit = this->_universalThermostatID.begin();
            thermit != this->_universalThermostatID.end();
            thermit++ )
       {
          checkpointfilestream << " CT\t" << 1+thermit->first
                               << "\t" << thermit->second << "\n";
       }
       for( map<int, double>::iterator Tit = this->_universalTargetTemperature.begin();
            Tit != this->_universalTargetTemperature.end();
            Tit++ )
       {
          if((0 >= Tit->first) || (0.0 >= Tit->second)) continue;
          checkpointfilestream << " ThT " << Tit->first << "\t" << Tit->second << "\n";
       }
    }
    else
    {
       checkpointfilestream << " T\t" << _universalTargetTemperature[0] << endl;
    }
    checkpointfilestream << " N\t" << _globalNumMolecules << endl;

    checkpointfilestream << " M\t" << "ICRVQD" << endl;
    checkpointfilestream.close();
  }
  
  domainDecomp->writeMoleculesToFile(filename, particleContainer); 
}
 
void Domain::readPhaseSpaceData(datastructures::ParticleContainer<Molecule>* particleContainer) {
  
  string token;

  double x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz;
  double Fx,Fy,Fz,Mx,My,Mz;
  unsigned int numcomponents=_components.size();
  unsigned long i,id;
  int componentid;
  
  x=y=z=vx=vy=vz=q1=q2=q3=Dx=Dy=Dz=0.;
  q0=1.;
  Fx=Fy=Fz=Mx=My=Mz=0.;
  _phaseSpaceFileStream >> token;
  if((token=="MoleculeFormat") || (token == "M")) {
    string ntypestring("ICRVQD");
    enum Ndatatype { ICRVQD, IRV, ICRV, ICRVFQDM } ntype=ICRVQD;

    if(_localRank==0) cout << "reading " << _globalNumMolecules << " molecules" << flush;
    if(_inpversion >= 51129) _phaseSpaceFileStream >> ntypestring;
    ntypestring.erase(ntypestring.find_last_not_of(" \t\n")+1);
    ntypestring.erase(0,ntypestring.find_first_not_of(" \t\n"));
    if (ntypestring=="ICRVFQDM")
      ntype=ICRVFQDM;
    else if (ntypestring=="ICRV")
      ntype=ICRV;
    else if (ntypestring=="IRV")
      ntype=IRV;
    if(_localRank==0) cout << " (" << ntypestring << ")" << flush;
    if(!numcomponents)
    {
      if(_localRank==0) cout << endl << "No components defined! Setting up single one-centered LJ" << endl;
      numcomponents=1;
      _components.resize(numcomponents);
      _components[0].setID(0);
      _components[0].addLJcenter(0.,0.,0.,1.,1.,1.);
    }
    //m_molecules.clear();
    for(i=0;i<_globalNumMolecules;++i)
    {
      if(ntype==ICRVFQDM)
        _phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz >> Fx >> Fy >> Fz
              >> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz >> Mx >> My >> Mz;
      else if(ntype==ICRV)
        _phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz;
      else if(ntype==IRV)
        _phaseSpaceFileStream >> id >> x >> y >> z >> vx >> vy >> vz;
      else
        _phaseSpaceFileStream >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
              >> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz;
      if((x<0.0 || x>=_globalLength[0] || y<0.0 || y>=_globalLength[1] || z<0.0 || z>=_globalLength[2]) && _localRank == 0) {
        cerr << endl << id << ": Molecule " << x << ";" << y << ";" << z << " out of box! " << flush;
      }
      // if(_localRank==0) cout << "processing molecule " << id << ".\n";

      if(((int)componentid > (int)numcomponents))
      {
        if (_localRank == 0) cerr << "Molecule id "
                                  << id << " has wrong componentid: "
                                  << componentid << ">" << numcomponents << endl;
        exit(1);
      }
      // @todo why do componetids start with 0?
      --componentid;
        
      //  store only those molecules within the domain of this process
      
      // @todo Pointer!!! new!!!  
      Molecule m1 = Molecule(id,componentid,x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz,&_components);
      particleContainer->addParticle(m1);
      //(_molecules.back()).setFM(Fx,Fy,Fz,Mx,My,Mz);
      _components[componentid].incrnumMolecules();
      // _globalRotDOF+=_components[componentid].rot_dof();
      
      if(!(i%2048)) if(_localRank==0) cout << '.' << flush;
    }
    if(_localRank==0) cout << " done" << endl;
      
    if(!_globalRho){
      _globalRho=_globalNumMolecules/(_globalLength[0]*_globalLength[1]*_globalLength[2]);
      if(_localRank==0) cout << "calculated global Rho:\t" << _globalRho << endl;
    }
    _phaseSpaceFileStream.close();
  }
  else {
    cout << token << "\n";
    _log.error("read input file", "Error in the PhaseSpace File"); 
    exit(1);
  }
}      
      
void Domain::initParameterStreams(double cutoffRadius){
  _comp2params.initialize(_components, _mixcoeff, _epsilonRF, cutoffRadius); 
}

void Domain::initFarFieldCorr(double cutoffRadius) {
  double UpotCorrLJ=0.;
  double VirialCorrLJ=0.;
  double MySelbstTerm=0.;
  unsigned int numcomp=_components.size();
  unsigned long nummolecules=0;
  for(unsigned int i=0;i<numcomp;++i) {
    Component& ci=_components[i];
    nummolecules+=ci.numMolecules();
    unsigned int numljcentersi=ci.numLJcenters();
    unsigned int numchargesi = ci.numCharges();
#ifdef COMPLEX_POTENTIAL_SET
    unsigned int numtersoffi = ci.numTersoff();
    unsigned int numdipolesi = ci.numDipoles();
#endif

    // effective dipoles computed from point charge distributions
    double chargeBalance[3];
    for(unsigned d = 0; d < 3; d++) chargeBalance[d] = 0;
    for(unsigned int si = 0; si < numchargesi; si++)
    {
       double tq = ci.charge(si).q();
       for(unsigned d = 0; d < 3; d++) chargeBalance[d] += tq * ci.charge(si).r()[d];
    }
#ifdef COMPLEX_POTENTIAL_SET
    // point dipoles
    for(unsigned int si=0;si<numdipolesi;++si)
    {
      double tmy = ci.dipole(si).absMy();
      double evect = 0;
      for(unsigned d = 0; d < 3; d++) evect += ci.dipole(si).e()[d] * ci.dipole(si).e()[d];
      double norm = 1.0 / sqrt(evect);
      for(unsigned d = 0; d < 3; d++) chargeBalance[d] += tmy * ci.dipole(si).e()[d] * norm;
    }
#endif
    double my2 = 0;
    for(unsigned d = 0; d < 3; d++) my2 += chargeBalance[d] * chargeBalance[d];
    MySelbstTerm += my2 * ci.numMolecules();

    for(unsigned int j=0;j<numcomp;++j) {
      Component& cj=_components[j];
#ifdef COMPLEX_POTENTIAL_SET
      unsigned numtersoffj = cj.numTersoff();
      // no LJ interaction between Tersoff components
      if(numtersoffi && numtersoffj) continue;
#endif
      unsigned int numljcentersj=cj.numLJcenters();
      ParaStrm& params=_comp2params(i,j);
      params.reset_read();
      // LJ centers
      for(unsigned int si=0;si<numljcentersi;++si) {
        double xi=ci.ljcenter(si).rx();
        double yi=ci.ljcenter(si).ry();
        double zi=ci.ljcenter(si).rz();
        double tau1=sqrt(xi*xi+yi*yi+zi*zi);
        for(unsigned int sj=0;sj<numljcentersj;++sj) {
          double xj=cj.ljcenter(sj).rx();
          double yj=cj.ljcenter(sj).ry();
          double zj=cj.ljcenter(sj).rz();
          double tau2=sqrt(xj*xj+yj*yj+zj*zj);
          if(tau1+tau2>=cutoffRadius){
            cerr << "Error calculating cutoff corrections, rc too small" << endl;
          }
          double eps24;
          params >> eps24;
          double sig2;
          params >> sig2;
          double fac=double(ci.numMolecules())*double(cj.numMolecules())*eps24;
          if(tau1==0. && tau2==0.){
            UpotCorrLJ+=fac*(TICCu(-6,cutoffRadius,sig2)-TICCu(-3,cutoffRadius,sig2));
            VirialCorrLJ+=fac*(TICCv(-6,cutoffRadius,sig2)-TICCv(-3,cutoffRadius,sig2));
          }
          else if(tau1!=0. && tau2!=0.) {
            UpotCorrLJ+=fac*(TISSu(-6,cutoffRadius,sig2,tau1,tau2)-TISSu(-3,cutoffRadius,sig2,tau1,tau2));
            VirialCorrLJ+=fac*(TISSv(-6,cutoffRadius,sig2,tau1,tau2)-TISSv(-3,cutoffRadius,sig2,tau1,tau2));
          }
          else {
            if(tau2==0.) {
              tau2=tau1;
            }
            UpotCorrLJ+=fac*(TICSu(-6,cutoffRadius,sig2,tau2)-TICSu(-3,cutoffRadius,sig2,tau2));
            VirialCorrLJ+=fac*(TICSv(-6,cutoffRadius,sig2,tau2)-TICSv(-3,cutoffRadius,sig2,tau2));
          }
        }
      }
    }
  }

  double fac=M_PI*_globalRho/(3.*_globalNumMolecules);
  UpotCorrLJ*=fac;
  VirialCorrLJ*=-fac;
        
  double epsRFInvrc3=2.*(_epsilonRF-1.)/((cutoffRadius*cutoffRadius*cutoffRadius)*(2.*_epsilonRF+1.));
  MySelbstTerm*=-0.5*epsRFInvrc3;

  _UpotCorr=UpotCorrLJ+MySelbstTerm;
  _VirialCorr=VirialCorrLJ+3.*MySelbstTerm;

  if(_localRank == 0)
  {
    cout << "far field terms:\nU\t" << _UpotCorr << "\nvirial\t" << _VirialCorr << "\n";
  }
}

#ifdef COMPLEX_POTENTIAL_SET
void Domain::specifyComponentSet(unsigned cosetid, double v[3], double tau, double ainit[3], double timestep)
{
   this->_localN[cosetid] = 0;
   for(unsigned d = 0; d < 3; d++)
   {
      this->_universalAdditionalAcceleration[d][cosetid] = ainit[d];
      this->_localVelocitySum[d][cosetid] = 0.0;
   }
   this->_universalTau[cosetid] = tau;
   if(this->_localRank == 0)
   {
      this->_globalN[cosetid] = 0;
      for(unsigned d = 0; d < 3; d++)
      {
         this->_globalTargetVelocity[d][cosetid] = v[d];
         this->_globalVelocitySum[d][cosetid] = 0.0;
      }
      this->_globalVelocityQueuelength[cosetid] = (unsigned)ceil(
         sqrt(this->_universalTau[cosetid] / (timestep*this->_universalConstantAccelerationTimesteps))
      );
      cout << "coset " << cosetid << " will receive "
           << this->_globalVelocityQueuelength[cosetid] << " velocity queue entries.\n";
   }
}

void Domain::determineAdditionalAcceleration
(
   parallel::DomainDecompBase* domainDecomp,
   datastructures::ParticleContainer<Molecule>* molCont, double dtConstantAcc )
{
   for( map<unsigned, double>::iterator uAAit = _universalAdditionalAcceleration[0].begin();
        uAAit != _universalAdditionalAcceleration[0].end();
        uAAit++ )
   {
      this->_localN[uAAit->first] = 0;
      for(unsigned d = 0; d < 3; d++)
         this->_localVelocitySum[d][uAAit->first] = 0.0;
   }
   for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
   {
      unsigned cid = thismol->componentid();
      map<unsigned, unsigned>::iterator uCSIDit = this->_universalComponentSetID.find(cid);
      if(uCSIDit == _universalComponentSetID.end()) continue;
      unsigned cosetid = uCSIDit->second;
      this->_localN[cosetid]++;
      for(unsigned d = 0; d < 3; d++)
         this->_localVelocitySum[d][cosetid] += thismol->v(d);
   }

   map<unsigned, long double>::iterator gVSit;
   //--- map<unsigned, long double> priorVelocitySum[3];
   if(!this->_localRank)
      for(gVSit = _globalVelocitySum[0].begin(); gVSit != _globalVelocitySum[0].end(); gVSit++)
      {
#ifndef NDEBUG
         cout << "required entries in velocity queue: " << _globalVelocityQueuelength[gVSit->first] << "\n";
         cout << "entries in velocity queue: " << _globalPriorVelocitySums[0][gVSit->first].size() << "\n";
#endif
         for(unsigned d = 0; d < 3; d++)
         {
            while(_globalPriorVelocitySums[d][gVSit->first].size() < _globalVelocityQueuelength[gVSit->first])
               _globalPriorVelocitySums[d][gVSit->first].push_back(_globalVelocitySum[d][gVSit->first]);
            //--- priorVelocitySum[d][gVSit->first] = _globalVelocitySum[d][gVSit->first];
         }
      }

   // cout << "local rank " << _localRank << " starting domainDecomp->collectCosetVelocity.\n";
   domainDecomp->collectCosetVelocity(&_localN, _localVelocitySum, &_globalN, _globalVelocitySum);
   // cout << "local rank " << _localRank << " returning from domainDecomp->collectCosetVelocity.\n";

   if(!this->_localRank)
   {
      for(gVSit = _globalVelocitySum[0].begin(); gVSit != _globalVelocitySum[0].end(); gVSit++)
      {
         double invgN = 1.0 / this->_globalN[gVSit->first];
         double invgtau = 1.0 / this->_universalTau[gVSit->first];
         double invgtau2 = invgtau * invgtau;
         double velocityDifferencePerUCAT[3];
         for(unsigned d = 0; d < 3; d++)
         {
            velocityDifferencePerUCAT[d] = invgN * ( _globalPriorVelocitySums[d][gVSit->first].front()
                                              - _globalVelocitySum[d][gVSit->first] )
                                                 / (double)(_globalVelocityQueuelength[gVSit->first]);
            _globalPriorVelocitySums[d][gVSit->first].pop_front();
            this->_universalAdditionalAcceleration[d][gVSit->first]
               += dtConstantAcc * invgtau2
                  * (_globalTargetVelocity[d][gVSit->first] - _globalVelocitySum[d][gVSit->first]*invgN)
                     + invgtau * velocityDifferencePerUCAT[d];
         }
#ifndef NDEBUG
         cout << "z-velocity difference per UCAT: " << velocityDifferencePerUCAT[2] << "\n";
#endif
      }
   }

   domainDecomp->broadcastCosetAcceleration(_universalAdditionalAcceleration);
   domainDecomp->broadcastVelocitySum(_globalVelocitySum);
}

double Domain::getDirectedVelocity(unsigned cosetid)
{
   double vv = 0.0;
   if(!this->_localRank)
   {
      for(unsigned d = 0; d < 3; d++)
      {
         double vd = this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
         vv += vd*vd;
      }
   }
   return sqrt(vv);
}
double Domain::getDirectedVelocity(unsigned cosetid, unsigned d)
{
   if(!this->_localRank) 
      return this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
   else return 0.0;
}

double Domain::getUniformAcceleration(unsigned cosetid)
{
   double aa = 0.0;
   for(unsigned d = 0; d < 3; d++)
      aa += this->_universalAdditionalAcceleration[d][cosetid]
            * this->_universalAdditionalAcceleration[d][cosetid];
   return sqrt(aa);
}
double Domain::getUniformAcceleration(unsigned cosetid, unsigned d)
{
   return this->_universalAdditionalAcceleration[d][cosetid];
}

double Domain::getMissingVelocity(unsigned cosetid, unsigned d)
{
   double vd = this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
   return this->_globalTargetVelocity[d][cosetid] - vd;
}

void Domain::setupProfile(unsigned xun, unsigned yun, unsigned zun)
{
   this->_universalNProfileUnits[0] = xun;
   this->_universalNProfileUnits[1] = yun;
   this->_universalNProfileUnits[2] = zun;
   for(unsigned d = 0; d < 3; d++)
   {
      _universalInvProfileUnit[d] = _universalNProfileUnits[d] / _globalLength[d];
   }
   this->resetProfile();
}

void Domain::considerComponentInProfile(int cid)
{
   this->_universalProfiledComponents[cid] = true;
}

void Domain::recordProfile(datastructures::ParticleContainer<Molecule>* molCont)
{
   int cid;
   unsigned xun, yun, zun, unID;
   double mv2, Iw2;
   for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
   {
      cid = thismol->componentid();
      if(this->_universalProfiledComponents[cid])
      {
         xun = (unsigned)floor(thismol->r(0) * this->_universalInvProfileUnit[0]);
         yun = (unsigned)floor(thismol->r(1) * this->_universalInvProfileUnit[1]);
         zun = (unsigned)floor(thismol->r(2) * this->_universalInvProfileUnit[2]);
         unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
                + yun * this->_universalNProfileUnits[2] + zun;
         this->_localNProfile[unID] += 1.0;
         for(int d=0; d<3; d++) this->_localvProfile[d][unID] += thismol->v(d);
         this->_localDOFProfile[unID] += 3.0 + (long double)(this->_components[cid].rot_dof());

         // record _twice_ the total (ordered + unordered) kinetic energy
         mv2 = 0.0;
         Iw2 = 0.0;
         thismol->calculate_mv2_Iw2(mv2, Iw2);
         this->_localKineticProfile[unID] += mv2+Iw2;
      }
   }
   this->_globalAccumulatedDatasets++;
}

void Domain::collectProfile(parallel::DomainDecompBase* domainDecomp)
{
   unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
                                                     * this->_universalNProfileUnits[2];
   for(unsigned unID = 0; unID < unIDs; unID++)
   {
      this->_globalNProfile[unID] = (double)this->_localNProfile[unID];
      for(int d=0; d<3; d++) this->_globalvProfile[d][unID] = (double)this->_localvProfile[d][unID];
      this->_globalDOFProfile[unID] = (double)this->_localDOFProfile[unID];
      this->_globalKineticProfile[unID] = (double)this->_localKineticProfile[unID];
   
      domainDecomp->reducevalues(&(this->_globalNProfile[unID]), &(this->_globalvProfile[0][unID]));
      domainDecomp->reducevalues(&(this->_globalvProfile[1][unID]), &(this->_globalvProfile[2][unID]));
      domainDecomp->reducevalues(&(this->_globalDOFProfile[unID]), &(this->_globalKineticProfile[unID]));
   }
}

void Domain::outputProfile(const char* prefix)
{
   if(this->_localRank) return;
   
   string vzpryname(prefix);
   string Tpryname(prefix);
   string rhpryname(prefix);
   rhpryname += ".rhpry";
   vzpryname += ".vzpry";
   Tpryname += ".Tpry";
   ofstream rhpry(rhpryname.c_str());
   ofstream vzpry(vzpryname.c_str());
   ofstream Tpry(Tpryname.c_str());
   if (!(vzpry && Tpry && rhpry))
   {
      return;
   }
   rhpry.precision(4);
   rhpry << "# y\trho\ttotal DOF\n# \n";
   vzpry.precision(4);
   vzpry << "# y\tvz\tv\n# \n";
   Tpry.precision(5);
   Tpry << "# y\t2Ekin/#DOF\n# \n";

   double layerVolume = this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]
                                               / this->_universalNProfileUnits[1];
   for(unsigned y = 0; y < this->_universalNProfileUnits[1]; y++)
   {
      double yval = (y + 0.5) / this->_universalInvProfileUnit[1];

      long double Ny = 0.0;
      long double DOFy = 0.0;
      long double twoEkiny = 0.0;
      long double velocitysumy[3];
      for(unsigned d = 0; d < 3; d++) velocitysumy[d] = 0.0;
      for(unsigned x = 0; x < this->_universalNProfileUnits[0]; x++)
      {
         for(unsigned z = 0; z < this->_universalNProfileUnits[2]; z++)
         {
            unsigned unID = x * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
                          + y * this->_universalNProfileUnits[2] + z;
            Ny += this->_globalNProfile[unID];
            DOFy += this->_globalDOFProfile[unID];
            twoEkiny += this->_globalKineticProfile[unID];
            for(unsigned d = 0; d < 3; d++) velocitysumy[d] += this->_globalvProfile[d][unID];
         }
      }

      if(Ny >= 64.0)
      {
         double vvdir = 0.0;
         for(unsigned d = 0; d < 3; d++)
         {
            double vd = velocitysumy[d] / Ny;
            vvdir += vd*vd;
         }
         rhpry << yval << "\t" << (Ny / (layerVolume * this->_globalAccumulatedDatasets))
               << "\t" << DOFy << "\n";
         vzpry << yval << "\t" << (velocitysumy[2] / Ny) << "\t" << sqrt(vvdir) << "\n";
         Tpry << yval << "\t" << (twoEkiny / DOFy) << "\n";
      }
      else
      {
         rhpry << yval << "\t0.000\t" << DOFy << "\n";
      }
   }

   rhpry.close();
   vzpry.close();
   Tpry.close();
}

void Domain::resetProfile()
{
   unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
                                                     * this->_universalNProfileUnits[2];
   for(unsigned unID = 0; unID < unIDs; unID++)
   {
      this->_localNProfile[unID] = 0.0;
      this->_globalNProfile[unID] = 0.0;
      for(int d=0; d<3; d++)
      {
         this->_localvProfile[d][unID] = 0.0;
         this->_globalvProfile[d][unID] = 0.0;
      }
      this->_localDOFProfile[unID] = 0.0;
      this->_globalDOFProfile[unID] = 0.0;
      this->_localKineticProfile[unID] = 0.0;
      this->_globalKineticProfile[unID] = 0.0;
   }
   this->_globalAccumulatedDatasets = 0;
}
#endif

// Helper functions ================================================================================
// used by Domain::init_Corr()

double TICCu(int n,double rc,double sigma2) {
  return -pow(rc,2*n+3) / (pow(sigma2,n)*(2*n+3));
}

double TICSu(int n,double rc,double sigma2,double tau){
  return -( pow(rc+tau,2*n+3) - pow(rc-tau,2*n+3) ) * rc / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3) ) +  ( pow(rc+tau,2*n+4) - pow(rc-tau,2*n+4) ) / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3)*(2*n+4) );
}
    
double TISSu(int n,double rc,double sigma2,double tau1,double tau2){
  double tauMinus,tauPlus;
  tauPlus = tau1+tau2;
  tauMinus = tau1-tau2;
  return -(   pow(rc+tauPlus,2*n+4) - pow(rc+tauMinus,2*n+4) - pow(rc-tauMinus,2*n+4) + pow(rc-tauPlus,2*n+4) ) * rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) ) +  (   pow(rc+tauPlus,2*n+5) - pow(rc+tauMinus,2*n+5) - pow(rc-tauMinus,2*n+5) + pow(rc-tauPlus,2*n+5) ) / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) );
}

double TICCv(int n,double rc,double sigma2){
  return 2*n * TICCu(n,rc,sigma2);
}

double TICSv(int n,double rc,double sigma2,double tau){
  return -( pow(rc+tau,2*n+2) - pow(rc-tau,2*n+2) ) * rc*rc / ( 4*pow(sigma2,n)*tau*(n+1) ) - 3*TICSu(n,rc,sigma2,tau);
}

double TISSv(int n,double rc,double sigma2,double tau1,double tau2){
  double tauMinus,tauPlus;
  tauPlus = tau1+tau2;
  tauMinus = tau1-tau2;
  return -(   pow(rc+tauPlus,2*n+3) - pow(rc+tauMinus,2*n+3) - pow(rc-tauMinus,2*n+3) + pow(rc-tauPlus,2*n+3) ) * rc*rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3) ) - 3*TISSu(n,rc,sigma2,tau1,tau2);
}
