#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <string>
#include <fstream>
#include <map>
#include <queue>
#include <list>

#include "molecules/Comp2Param.h"
#include "molecules/Component.h"
#include "utils/Log.h"

class Molecule;

namespace datastructures{
  template<class ParticleType>
  class ParticleContainer;
}

#ifdef GRANDCANONICAL
namespace ensemble
{
   class ChemicalPotential;
}
#endif

namespace parallel{
  class DomainDecompBase; 
}


using namespace std;

//! @brief This class is used to read in the phasespace and to handle macroscopic values
//!
//! The main purpose of this class is responsible for all macroscopic value (Potential,...).
//! It is important to differentiate between local and global version of those values
//! Macroscopic values are values that aggregate some property of a set of molecules.
//! As this program is designed to run on a parallel computer, there are typically 
//! several processes. Each process has an instance of this class, but with a different 
//! subset of molecules. Therefore, also the macroscopic values are only representative 
//! for the "local" domain of that process. Thus, local member variables that represent 
//! "local" macroscopic values, begin with _local...
//! 
//! At some points of the simulation, macroscopic values for the whole set of molecules 
//! have to be calculated. Those values are stored in member variables beginning with _global. 
//!
//!  
//! The phasespace file consists of two parts, a header (including component description)
//! and the molecule data. Header and data have to be read in with two different methods. 
//! Thats' because the molecules from the data part have to be stored in a molecule 
//! container, which can only be created with information from the header part.
//! 
//! Example for the Phasespace file:
//! <TABLE><TR><TD><CODE>
//! MDProject 20070123\n
//! currenttime "t (double)" \n
//! Temperature "T (double)" \n
//! Length "x (double)" "y (double)" "z (double)" \n
//! NumberOfComponents "C (int)" \n
//!    [#LJ]_C1 [#DP]_C1 [#QP]_C1
//!    [xlj_1]_C1 [ylj_1]_C1 [zlj_1]_C1 [m_1]_C1 [eps_1]_C1 [sigma_1]_C1 \n
//!    ... \n
//!    [xlj_#LJ]_C1 [ylj_#LJ]_C1 [zlj_#LJ]_C1 [m_#LJ]_C1 [eps_#LJ]_C1 [sigma_#LJ]_C1 \n
//!    [xdp_1]_C1 [ydp_1]_C1 [zdp_1]_C1 [eMyx_1]_C1 [eMyy_1]_C1 [eMyz_1]_C1 [absMy_1]_C1 \n
//!    ... \n
//!    [xdp1_#DP]_C1 [ydp1_#DP]_C1 [zdp1_#DP]_C1 [eMyx_#DP]_C1 [eMyy_#DP]_C1 [eMyz_#DP]_C1 [absMy_#DP]_C1 \n
//!    [xqp_1]_C1 [yqp_1]_C1 [zqp_1]_C1 [eMyx_1]_C1 [eMyy_1]_C1 [eMyz_1]_C1 [absMy_1]_C1 \n
//!    ... \n
//!    [xqp1_#QP]_C1 [yqp1_#QP]_C1 [zqp1_#QP]_C1 [eQx_#QP]_C1 [eQy_#QP]_C1 [eQz_#QP]_C1 [absQ_#QP]_C1 \n
//!    [I11]_C1 [I22]_C1 [I33]_C1 \n
//!    
//!    ... \n
//!  
//!    [#LJ]_Cn [#DP]_Cn [#QP]_Cn
//!    [xlj_1]_Cn [ylj_1]_Cn [zlj_1]_Cn [m_1]_Cn [eps_1]_Cn [sigma_1]_Cn \n
//!    ... \n
//!    [xlj_#LJ]_Cn [ylj_#LJ]_Cn [zlj_#LJ]_Cn [m_#LJ]_Cn [eps_#LJ]_Cn [sigma_#LJ]_Cn \n
//!    [xdp_1]_Cn [ydp_1]_Cn [zdp_1]_Cn [eMyx_1]_Cn [eMyy_1]_Cn [eMyz_1]_Cn [absMy_1]_Cn \n
//!    ... \n
//!    [xdp1_#DP]_Cn [ydp1_#DP]_Cn [zdp1_#DP]_Cn [eMyx_#DP]_Cn [eMyy_#DP]_Cn [eMyz_#DP]_Cn [absMy_#DP]_Cn \n
//!    [xqp_1]_Cn [yqp_1]_Cn [zqp_1]_Cn [eMyx_1]_Cn [eMyy_1]_Cn [eMyz_1]_Cn [absMy_1]_Cn \n
//!    ... \n
//!    [xqp1_#QP]_Cn [yqp1_#QP]_Cn [zqp1_#QP]_Cn [eQx_#QP]_Cn [eQy_#QP]_Cn [eQz_#QP]_Cn [absQ_#QP]_Cn \n
//!    [I11]_Cn [I22]_Cn [I33]_Cn \n
//!
//! NumberOfMolecules "N" \n
//! [id_1] [type_1] [x_1] [y_1] [z_1] [vx_1] [vy_1] [vz_1] [q0_1] [q1_1] [q2_1] [q3_1] [Dx_1] [Dy_1] [Dz_1]\n
//! [id_2] [type_2] [x_2] [y_2] [z_2] [vx_2] [vy_2] [vz_2] [q0_2] [q1_2] [q2_2] [q3_2] [Dx_2] [Dy_2] [Dz_2]\n
//! ... \n
//! [id_N] [type_N] [x_N] [y_N] [z_N] [vx_N] [vy_N] [vz_N] [q0_N] [q1_N] [q2_N] [q3_N] [Dx_N] [Dy_N] [Dz_N]\n 
//! </CODE></TD></TR></TABLE>
class Domain{
  public:
    //! The constructor sets _localRank to rank and initializes all member variables 
    Domain(int rank);

    //! @brief gets a filename and opens an ifstream associated with the given file
    //! 
    //! As the reading of the phasespace file is separated into two parts,
    //! but each line of the file should only be parsed once, not the filename
    //! itself is stored, but a stream (_phaseSpaceFileStream) which is associated with
    //! the file
    //! @param filename full path to the input file
    void setPhaseSpaceFile(string filename);  
        
    //! @brief reads in header of the input file (including component description)
    //!
    //! The Header in the input file consists of several elements. An element starts
    //! in a new line with the element name followed by some whitespace (e.g. "Temperature ").
    //! After that, there can be any number of tokens belonging to that element. 
    //! A description of the different elements follows below. But first
    //! some rules dealing with the order of the elements:
    //! \li The first header element must start with the string "MDProject" 
    //!     followed by a timestamp.
    //! \li The last header element must start with the string NumberOfMolecules.
    //!     After that, the header is over and the molecule data follows
    //! \li The order of the remaining header lines is not important.
    //! \li Header elements beginning with "#" are ignored. 
    //!
    //! Now the description of the different elements:
    //! \li MDProject: One token follows with a version number (see description of _inpversion);
    //! \li currenttime: One token follows with the start time
    //! \li Temperature: One token follows with the temperature
    //! \li Length: Three tokens follow with the length of the simulation box
    //!     in the three dimensions (x, y and z)
    //! \li NumberOfComponents: Here follow several tokens, the first one is the
    //!     actual Number of Components.
    //!     Then the values describing the components have to follow, seperated by 
    //!     whitespace. For each component, the following values have to be provided:
    //!     - Number of Lennard-Jones-Centers, Number of Dipoles, Number of Quadrupoles (all int)
    //!     - For each LJ-Center: x-coord., y-coord., z-coord., mass, epsilon, sigma (all double)
    //!     - For each Dipole: x-coord., y-coord., z-coord., eMyx, eMyy, eMyz, absMy (all double)
    //!     - For each Quadrupole: x-coord., y-coord., z-coord., eQx, eQy, eQz, absQ (all double)
    //!     - moments of inertia for principal axes: I11, I22, I33 (all double)
    //!     - For each pair of different components: xi, eta (both double)
    //!     - epsilonRF (double)
    //! \li NumberOfMolecules: One token follows with the number of molecules
    //!
    //! An example can be seen in the documentation of this class
    void readPhaseSpaceHeader(
       double timestep
#ifdef TRUNCATED_SHIFTED
          ,
       double cutoff
#endif
#ifdef COMPLEX_POTENTIAL_SET
          ,
       double cutoffTersoff
#endif
    );

    //! @brief reads in the data of all molecules
    //! 
    //! The Molecule Data starts in a new line with the string "MoleculeFormat"
    //! followed by whitespace and a string representing the format.
    //! In the standard case (format ICRVQD), the following values are provided for each molecule:
    //! \li id of the molecule (int)
    //! \li id of the component of the molecule (int)
    //! \li Coordinates: x, y, z (all double)
    //! \li velocities: vx, vy, vz (all double)
    //! \li Orientation (quaternion): q0, q1, q2, q3 (all double)
    //! \li Angular Momentum: Dx, Dy, Dz (all double)
    //!
    //! returns the maximal molecule ID
    //! 
    //! An example can be seen in the documentation of this class
    //! @param particleContainer Here the Molecules from the input file are stored 
    unsigned long readPhaseSpaceData(datastructures::ParticleContainer<Molecule>* particleContainer
#ifdef TRUNCATED_SHIFTED
          ,
       double cutoffRadius
#endif
#ifdef GRANDCANONICAL
          ,
       list< ensemble::ChemicalPotential >* lmu
#endif
);

    //! @brief writes a checkpoint file that can be used to continue the simulation
    //!
    //! The format of the checkpointfile written by this method is the same as the format
    //! of the input file. 
    //! @param filename Name of the checkpointfile (including path)
    //! @param particleContainer The molecules that have to be written to the file are stored here
    //! @param domainDecomp In the parallel version, the file has to be written by more than one process.
    //!                     Methods to achieve this are available in domainDecomp 
    void writeCheckpoint(string filename, datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, double dt);

    //! @brief initialize far field correction parameters
    //! 
    //! By limiting the calculation to pairs of particles which have
    //! less distance than the given cutoff radius, an error is made
    //! By calculating approximations for the neglected pairs, the
    //! error can be reduced
    //! @param cutoffRadius cutoff radius
    //! @todo How does the Correction work? Give reference to some paper,
    //!       documentation in the implementation
    void initFarFieldCorr(double cutoffRadius); // CHECKED

    //! @brief initialize parameter streams
    //!
    //! This method should only be called, after the the component information
    //! and all molecule data have been read in
    //! @param cutoffRadius cutoff radius
    void initParameterStreams(double cutoffRadius);

    //! @brief set the potential of the local process
    void setLocalUpot(double Upot); 

    //! @brief get the potential of the local process
    double getLocalUpot() const;

    //! @brief set the virial of the local process
    void setLocalVirial(double Virial); 

    //! @brief get the virial of the local process
    double getLocalVirial() const;   
    
    //! @brief get thermostat scaling for translations 
    double getGlobalBetaTrans();
    double getGlobalBetaTrans(int thermostat);

    //! @brief get thermostat scaling for rotations 
    double getGlobalBetaRot();
    double getGlobalBetaRot(int thermostat);

    double T(int thermostat) { return this->_globalTemperatureMap[thermostat]; }
    double targetT(int thermostat) { return this->_universalTargetTemperature[thermostat]; }

    //! @brief return the length of the domain
    //!
    //! @param index dimension for which the length should be returned
    double getGlobalLength(int index) const;

    //! @brief get the global pressure
    //!
    //! @todo provide justification for the formula
    double getGlobalPressure();

    //! @brief get the global average potential per particle
    //!
    //! Before this method is called, it has to be sure that the
    //! global potential has been calculated (method calculateGlobalValues)
    double getAverageGlobalUpot() const;

    //! @brief get the global average virial per particle
    //!
    //! Before this method is called, it has to be sure that the
    //! global virial has been calculated (method calculateGlobalValues)
    double getAverageGlobalVirial() const;

    //! @brief sets _localSummv2 to the given value
    void setLocalSummv2(double summv2);
    void setLocalSummv2(double summv2, int thermostat);

    //! @brief sets _localSumIw2 to the given value
    void setLocalSumIw2(double sumIw2);
    void setLocalSumIw2(double sumIw2, int thermostat);

    //! @brief sets _localThermostatN(i) and _localRotationalDOF(i)
    void setLocalNrotDOF(int i, unsigned long N, unsigned long rotDOF)
    {
       this->_localThermostatN[i] = N;
       this->_localRotationalDOF[i] = rotDOF;
    }

    unsigned getComponentRotDOF(int cid) { return this->_components[cid].rot_dof(); }

    //! @brief get the current time
    double getCurrentTime();
    
    //! @brief advance the current time by timestep
    void advanceTime(double timestep);

    //! @brief get a reference to the vector of Components
    vector<Component>& getComponents();
    
    //! @brief get the parameter streams
    Comp2Param& getComp2Params();
       
    //! @brief calculate the global macroscopic values
    //!
    //! @param domainDecomp domain decomposition
    //! @param particleContainer particle Container
    //! @todo more detailed description
#ifdef COMPLEX_POTENTIAL_SET
    void calculateGlobalValues( parallel::DomainDecompBase* domainDecomp,
                                datastructures::ParticleContainer<Molecule>* particleContainer,
                                bool collectThermostatVelocity, double Tfactor );
#else
    void calculateGlobalValues( parallel::DomainDecompBase* domainDecomp,
                                datastructures::ParticleContainer<Molecule>* particleContainer,
                                double Tfactor );
#endif
    
    //! @brief calculate _localSummv2 and _localSumIw2
    //! @todo more detailed description
#ifdef COMPLEX_POTENTIAL_SET
    void calculateThermostatDirectedVelocity(datastructures::ParticleContainer<Molecule>* partCont);
    bool thermostatIsUndirected(int th) { return this->_universalUndirectedThermostat[th]; }
    double thermostatv(int th, int d) { return this->_universalThermostatDirectedVelocity[d][th]; }
#endif
    void calculateVelocitySums(datastructures::ParticleContainer<Molecule>* partCont);

    int ownrank() { return this->_localRank; }

    bool severalThermostats() { return this->_universalComponentwiseThermostat; }
    int getThermostat(int cid) { return this->_universalThermostatID[cid]; }
    unsigned maxThermostat()
    {
       return (_universalComponentwiseThermostat)? (_universalThermostatN.size() - 2): 0;
    }

#ifdef COMPLEX_POTENTIAL_SET
    /// assigns a coset ID to a component (ID)
    void assignCoset(unsigned cid, unsigned cosetid) { _universalComponentSetID[cid] = cosetid; }
    /// sets the information on the acceleration model for one coset
    void specifyComponentSet(unsigned cosetid, double v[3], double tau, double ainit[3], double timestep);
    /// sets the number of timesteps between two updates of the uniform acceleration
    void setUCAT(unsigned uCAT) { this->_universalConstantAccelerationTimesteps = uCAT; }
    /// returns the number of timesteps between two updates of the uniform acceleration
    unsigned getUCAT() { return this->_universalConstantAccelerationTimesteps; }
    /// are there any cosets?
    bool isAcceleratingUniformly() { return ( (this->_universalTau.size() > 0)
                                              && (this->_universalConstantAccelerationTimesteps > 0) ); }
    /// updates the intensity and direction of the uniform acceleration
    void determineAdditionalAcceleration
    (
        parallel::DomainDecompBase* domainDecomp,
        datastructures::ParticleContainer<Molecule>* molCont, double dtConstantAcc
    );
    /// returns the acceleration map (necessary for passing data to the integrator)
    map<unsigned, double>* getUAA() { return this->_universalAdditionalAcceleration; }
    /// returns the cosetid of a component (0 for unaccelerated components)
    unsigned getComponentSet(unsigned cid)
    {
       if(_universalComponentSetID.find(cid) == _universalComponentSetID.end()) return 0;
       else return this->_universalComponentSetID[cid];
    }
    double getDirectedVelocity(unsigned cosetid);
    double getDirectedVelocity(unsigned cosetid, unsigned d);
    double getUniformAcceleration(unsigned cosetid);
    double getUniformAcceleration(unsigned cosetid, unsigned d);
    /// returns the difference between the desired velocity and the global average velocity
    double getMissingVelocity(unsigned cosetid, unsigned d);
    double getCosetN(unsigned cosetid) { return this->_globalN[cosetid]; }

    void setupProfile(unsigned xun, unsigned yun, unsigned zun);
    void considerComponentInProfile(int cid); 
    void recordProfile(datastructures::ParticleContainer<Molecule>* molCont);
    void collectProfile(parallel::DomainDecompBase* domainDecomp);
    void outputProfile(const char* prefix);
    void resetProfile();

    unsigned maxCoset() { return this->_universalTau.size(); }
#endif

    double N() {return this->_globalNumMolecules;}
    double N(unsigned cid) { return this->_components[cid].numMolecules(); }
#ifdef GRANDCANONICAL
    void Nadd(unsigned cid, int N, int localN);
#endif

   double getGlobalLength(int d) { return this->_globalLength[d]; }

   void setupRDF(double interval, unsigned bins);
   void resetRDF();
   void collectRDF(parallel::DomainDecompBase* domainDecomp);
   void outputRDF(const char* prefix, unsigned i, unsigned j);
   void accumulateRDF();
   void tickRDF() { this->_universalRDFTimesteps++; }
   inline void observeRDF(unsigned i) { this->_localCtr[i] ++; }
   inline void observeRDF(double dd, unsigned i, unsigned j)
   {
      if(dd > this->ddmax) return;
      if(i > j) { this->observeRDF(dd, j, i); return; }
      unsigned l = (unsigned)(sqrt(dd)/this->_universalInterval);
      this->_localDistribution[i][j-i][l] ++;
   }
   void thermostatOff() { this->_universalNVE = true; }
   void thermostatOn() { this->_universalNVE = false; }
   bool NVE() { return this->_universalNVE; }
   bool thermostatWarning() { return (this->_universalSelectiveThermostatWarning > 0); }

#ifdef GRANDCANONICAL
   void evaluateRho(unsigned long localN, parallel::DomainDecompBase* comm);
#endif

  private:
    //! Logging interface
    static utils::Log _log;  
  
    //! rank of the local process
    int _localRank;

    //! file stream associatted to the input file
    ifstream _phaseSpaceFileStream;

    //! @brief Version of the input file
    //!
    //! even though it is desirable, that the format of the input file
    //! doesn't change, is sometimes does change. When that happens,
    //! the code which reads in the input file (parser) has to be changed as well.
    //! old versions of the input file then can't be read any more.
    //! So whenever the parser is changed, _inpversion is set to the
    //! date of the change (YYYYMMDD) (hard-coded). Only input files
    //! with the same version are sure to be processed correctly
    unsigned long _inpversion;
    
    //! Potential of the local process
    double _localUpot;
    //! Virial of the local process
    double _localVirial;
    //! global Potential
    double _globalUpot;
    //! global virial
    double _globalVirial;
    // //! @todo WHAT IS THIS? 
    // double _globalBetaTrans;
    // //! @todo WHAT IS THIS?
    // double _globalBetaRot;
    //! global density
    double _globalRho;
    // //! global Temperature
    // double _globalTemperature;
    // //! global Number of rotational degrees of freedom
    // unsigned long _globalRotDOF;
    //! global Number of Molecules
    unsigned long _globalNumMolecules;
    //! side length of the cubic simulation box
    double _globalLength[3];

    //! does a componentwise thermostat apply?
    bool _universalComponentwiseThermostat;
    //! thermostat IDs. negative: no thermostat, 0: global, positive: componentwise
    //! in the case of a componentwise thermostat, all components are assigned
    //! a thermostat ID different from zero.
    map<int, int> _universalThermostatID;
    //! _localThermostatN[0] and _universalThermostatN[0] are always the total number
    //! of particles in the subdomain and, respectively, the entire domain
    map<int, unsigned long> _localThermostatN;
    map<int, unsigned long> _universalThermostatN;
    map<int, unsigned long> _localRotationalDOF;
    map<int, unsigned long> _universalRotationalDOF;
    //! _globalTemperatureMap[0] is always the temperature of the whole system,
    //! including components to which no thermostat is applied.
    //! The _globalTemperatureMap stores actual CURRENT temperatures, whereas
    //! the temperature objective of the thermostat is stored in _universalTargetTemperature
    map<int, double> _globalTemperatureMap;
    map<int, double> _universalTargetTemperature;
    map<int, double> _universalBTrans;
    map<int, double> _universalBRot;
#ifdef COMPLEX_POTENTIAL_SET
    //! should the directed movement be subtracted when calculating the temperature?
    map<int, bool> _universalUndirectedThermostat;
    //! stores the velocity of the directed movement
    map<int, double> _universalThermostatDirectedVelocity[3];
    map<int, double> _localThermostatDirectedVelocity[3];
#endif
    bool _universalNVE;

#ifdef COMPLEX_POTENTIAL_SET
    /// calculate new value of the uniform acceleration each # timesteps
    unsigned _universalConstantAccelerationTimesteps;
    /// assigns a component set ID to some of the components
    map<unsigned, unsigned> _universalComponentSetID;
    /// local number of molecules that belong to a given component set ID
    map<unsigned, unsigned long> _localN;
    /// global number of molecules that belong to a given component set ID
    map<unsigned, unsigned long> _globalN;
    /// local sum of the velocity vectors corresponding to a given component set ID
    map<unsigned, long double> _localVelocitySum[3];
    /// global sum of the velocity vectors corresponding to a given component set ID
    map<unsigned, long double> _globalVelocitySum[3];
    /// uniform acceleration
    map<unsigned, double> _universalAdditionalAcceleration[3];
    /// target average velocity for the molecules of a coset
    map<unsigned, double> _globalTargetVelocity[3];
    /// delay variable tau of a coset
    map<unsigned, double> _universalTau;
    /// queue of previously recorded velocity sums
    map<unsigned, deque<long double> > _globalPriorVelocitySums[3];
    /// number of items in the velocity queue
    map<unsigned, unsigned> _globalVelocityQueuelength;

    //! 1 / dimension of a profile cuboid
    double _universalInvProfileUnit[3];
    //! number of successive profile cuboids in x/y/z direction
    unsigned _universalNProfileUnits[3];
    //! local N profile map
    map<unsigned, long double> _localNProfile;
    //! global N profile map
    map<unsigned, double> _globalNProfile;
    //! local directed velocity profile map
    map<unsigned, long double> _localvProfile[3];
    //! global directed velocity  profile map
    map<unsigned, double> _globalvProfile[3];
    //! local kinetic energy profile map
    map<unsigned, long double> _localKineticProfile;
    //! global kinetic energy profile map
    map<unsigned, double> _globalKineticProfile;
    //! local counter w. r. t. degrees of freedom
    map<unsigned, long double> _localDOFProfile;
    //! global counter w. r. t. degrees of freedom
    map<unsigned, double> _globalDOFProfile;
    //! how many _evaluated_ timesteps are currently accumulated in the profile?
    unsigned _globalAccumulatedDatasets;
    //! which components should be considered?
    map<unsigned, bool> _universalProfiledComponents;
#endif

    bool _doCollectRDF;
    double _universalInterval;
    unsigned _universalBins;
    unsigned _universalRDFTimesteps, _universalAccumulatedTimesteps;
    double ddmax;
    unsigned long *_localCtr, *_globalCtr, *_globalAccumulatedCtr;
    unsigned long ***_localDistribution, ***_globalDistribution, ***_globalAccumulatedDistribution;

    int _universalSelectiveThermostatCounter;
    int _universalSelectiveThermostatWarning;
    int _universalSelectiveThermostatError;

    //! local sum (over all molecules) of the mass multiplied with the squared velocity
    map<int, double> _local2KETrans;
    //! local sum (over all molecules) of the moment of inertia 
    //! multiplied with the squared  rotational velocity
    map<int, double> _local2KERot; 
    
    //! reaction field
    //! @todo WHAT IS THIS?
    //! @todo local or global?
    double _epsilonRF;

    //! Potential correction for the error made by the cutoff
    //! @todo local or global?
    double _UpotCorr;
    //! Virial correction for the error made by the cutoff
    //! @todo local or global?
    double _VirialCorr;

    //! @todo comment
    double _currentTime;

    //! Components resp. molecule types
    vector<Component> _components;
    //! parameter streams for each possible pair of molecule-types
    Comp2Param _comp2params;
    //! modified Lorentz-Berthelot mixing rule parameters
    //! @todo more explanation
    vector<double> _mixcoeff;
};


#endif /*DOMAIN_H_*/
