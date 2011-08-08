/*
 * Martin Horsch, ls1 project moderated by Martin Bernreuther
 * (C)2011 GNU General Public License
 */
#ifndef GRANDCANONICAL_H
#define GRANDCANONICAL_H

#include <list>
#include "molecules/Molecule.h"

#define DELETION_TRUE 1
#define DELETION_FALSE 0
#define DELETION_INVALID -1

namespace datastructures{
   template<class ParticleType> class ParticleContainer;
}
namespace ensemble{
   class Random;
#ifdef GRANDCANONICAL
   class ChemicalPotential;
#endif
}
namespace parallel{
  class DomainDecompBase; 
}
using namespace std;

typedef datastructures::ParticleContainer<Molecule> TMoleculeContainer;

class ensemble::Random
{
 public:
   Random();
   void init(int seed, int rank);

   float rnd_muVT();

   int getIX() { return this->ix_muVT; }

 private:
   int ix_muVT, iy_muVT;
   float am_muVT;
   int ownrank;  // only for debugging purposes (indicate rank in console output)
};

#ifdef GRANDCANONICAL
class ensemble::ChemicalPotential
{
 public:
   ChemicalPotential();

   void setMu(int cid, double chempot) { this->mu = chempot; this->componentid = cid; }
   unsigned getInterval() { return this->interval; }
   void setInterval(unsigned delta) { this->interval = delta; }
   void setInstances(unsigned n) { this->instances = n; }
   void setSystem(double x, double y, double z, double m);
   void setGlobalN(unsigned long N) { this->globalN = N; }
   void setNextID(unsigned long id) { this->nextid = id; }
   void setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1);
   void setIncrement(unsigned idi) { this->id_increment = idi; }

   void prepareTimestep(TMoleculeContainer* cell, parallel::DomainDecompBase* comm);  // C must not contain the halo!

   // false if no deletion remains for this subdomain
   int getDeletion(TMoleculeContainer* cell, double* minco, double* maxco);
   unsigned long getInsertion(double* ins);  // 0 if no insertion remains for this subdomain
   bool decideDeletion(double deltaUTilde);
   bool decideInsertion(double deltaUTilde);

   Molecule loadMolecule();
   void storeMolecule( Molecule old )
   {
      assert(old.componentid() == componentid);
#ifndef NDEBUG
      old.check(old.id());
#endif
      this->reservoir.push_back(old);
   }
   bool hasSample() { return !this->reservoir.empty(); }

   void setPlanckConstant(double h_in) { this->h = h_in; }
   void submitTemperature(double T);
   void setControlVolume(
      double x0, double y0, double z0, double x1, double y1, double z1
   );

   unsigned long getGlobalN() { return this->globalN; }
   double getGlobalRho() { return (double)(this->globalN) / this->globalV; }

   void outputIX() { cout << "  r" << ownrank << "[IX" << rnd.getIX() << "]  "; }

   void assertSynchronization(parallel::DomainDecompBase* comm);

   double getMu() { return this->mu; }
   int getComponentID() { return this->componentid; }
   int rank() { return this->ownrank; }

   void disableWidom() { this->widom = false; }
   void enableWidom() { this->widom = true; }
   bool isWidom() { return this->widom; }

   double getLambda() { return this->lambda; }
   float getDensityCoefficient() { return this->decisive_density; }

 private:
   int ownrank;  // only for debugging purposes (indicate rank in console output)

   double h;  // Plancksches Wirkungsquantum

   double mu;
   double muTilde;
   int componentid;
   unsigned interval;  // how often?
   unsigned instances;  // how many trial insertions and deletions?
   ensemble::Random rnd;
   double system[3];  // extent of the system
   float minredco[3];  // minimal coordinates of the subdomain reduced w. r. t. the system size
   float maxredco[3];   // maximal coordinates of the subdomain reduced w. r. t. the system size

   unsigned long nextid;  // ID given to the next inserted particle
   unsigned id_increment;
   list<unsigned> remainingDeletions;  // position of the vector that should be deleted
   list<double> remainingInsertions[3];
   list<unsigned long> remainingInsertionIDs;
   list<float> remainingDecisions;  // first deletions, then insertions

   unsigned long globalN;
   double globalV;
   double molecularMass;
   double globalReducedVolume;

   bool restrictedControlVolume;
   double control_bottom[3];
   double control_top[3];

   float decisive_density;

   double lambda;

   bool widom;  // Widom method -> determine mu by test insertions which are all rejected

   list<Molecule> reservoir;
};
#endif

#endif
