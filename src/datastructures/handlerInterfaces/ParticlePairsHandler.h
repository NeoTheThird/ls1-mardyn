#ifndef PARTICLEPAIRSHANDLER_H_
#define PARTICLEPAIRSHANDLER_H_

namespace datastructures {
  template<class ParticleType>
  class ParticlePairsHandler; 
}

//! @brief interface for defining the action performed when processing a pair
//! 
//! The idea of a ParticleContainer is, that the container itself only knows
//! about how to efficiently store and access Particles (or neighbouring pairs
//! of Particles) It doesn't know anything (exception follows) about the 
//! internal structure of a Particle, or what action should be performed for 
//! a pair of particles (an action could e.g. be to calculate the force).
//! The only thing that has to be known from a particle is it's position.
//! 
//! An application e.g. wants to find all neighbouring particle pairs and
//! calculate the forces between them. The retrieval of the pairs has to be
//! done by the ParticleContainer, but the force calculation has to be
//! performed somewhere else. That's where this interface comes into play.
//! There are typically three things to be done:
//! \li Do some initial stuff before the pair processing
//! \li Do something for each pair
//! \li Do something after all pairs have been processed
//!
//! The ParticleContainer has an instance of a class implementing this interface
//! as a member variable. When the method of the ParticleContainer that traverses
//! the pairs is called, it first has to call the method "init()" of it's
//! ParticlePairsHandler, then for each pair "processPair(...)" and at the end
//! "finish()".
//! A class implementing this interface now serves as an adapter between the
//! particleContainer an some other part of the programm (e.g. force calculation)
template<class ParticleType>
class datastructures::ParticlePairsHandler{
  public:
    //! Constructor
    ParticlePairsHandler(){
    }
  
    //! Destructor
    virtual ~ParticlePairsHandler(){
    }
    
    //! @brief things to be done before particle pairs are processed
    virtual void init() = 0;
    
    //! @brief things to be done after particle pairs are processed
    virtual void finish() = 0;
    
    //! @brief returns Lennard-Jones + polarity intermolecular interaction energy
    //!
    //! @param particle1 first particle
    //! @param particle2 second particle
    //! @param distanceVector[3] distance between the two particles
    //! @param pairType describes whether the pair is am original pair(0), a duplicated pair(1),
    //!                 or a pair containing a grand canonical test particle (2)
    //!                 for details about pair types see comments on traversePairs() in ParticleContainer
    virtual double processPair(ParticleType& particle1, ParticleType& particle2, double distanceVector[3], int pairType, double dd, bool cLJ) = 0;
#ifdef COMPLEX_POTENTIAL_SET
    virtual void preprocessTersoffPair(Molecule& particle1, Molecule& particle2, bool pairType) = 0;
    // virtual void processTersoffPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], int pairType) = 0;
    virtual void processTersoffAtom(Molecule& particle1, double params[15], double delta_r) = 0;
#endif
    virtual void recordRDF() = 0;
    virtual void enableWallLJ() = 0;
    virtual void disableWallLJ() = 0;
};

#endif /*PARTICLEPAIRSHANDLER_H_*/
