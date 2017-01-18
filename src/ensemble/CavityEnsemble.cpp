
#include "CavityEnsemble.h"

#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
#include "molecules/Quaternion.h"

#define COMMUNICATION_THRESHOLD 3

using Log::global_log;

CavityEnsemble::CavityEnsemble()
{
   this->ownrank = -1;
   this->initialized = false;
   this->rotated = false;
   this->interval = 1000;
   this->componentid = (unsigned)-1;
   this->T = -1.0;

   this->globalV = 0.0;
   this->system[0] = 0.0; system[1] = 0.0; system[2] = 0.0;
   this->minredco[0] = 0.0; minredco[1] = 0.0; minredco[2] = 0.0;
   this->maxredco[0] = 1.0; maxredco[1] = 1.0; maxredco[2] = 1.0;

   this->restrictedControlVolume = false;
   this->control_bottom[0] = 0.0; control_bottom[1] = 0.0; control_bottom[2] = 0.0;
   this->control_top[0] = 1.0; control_top[1] = 1.0; control_top[2] = 1.0;

   this->active = set<unsigned long>();
   this->reservoir = map<unsigned long, Molecule*>();
   this->globalActive = 0;

   this->boundarySpecified = false;
   
   this->maxNeighbours = 1024;
   this->r2n = 0.0;
   this->idoffset = 0;
}

void CavityEnsemble::setSystem(double x, double y, double z)
{
   this->system[0] = x; this->system[1] = y; this->system[2] = z;
   if(!this->restrictedControlVolume)
   {
      this->globalV = x*y*z;
	    
      for(int d=0; d < 3; d++)
      {
         this->control_bottom[d] = 0.0;
         this->control_top[d] = this->system[d];
      }
   }
}

void CavityEnsemble::setSubdomain(int rank, double x0, double x1, double y0, double y1, double z0, double z1)
{
   this->ownrank = rank;
   this->async.init(8623);

   for(int d=0; d < 3; d++)
      if(control_bottom[d] >= control_top[d])
         this->restrictedControlVolume = false;

   if(!this->restrictedControlVolume)
   {
      this->globalV = this->system[0] * system[1] * system[2];
    
      for(int d=0; d < 3; d++)
      {
         this->control_bottom[d] = 0.0;
         this->control_top[d] = this->system[d];
      }
   }
   this->minredco[0] = (x0 - control_bottom[0]) /
                          (control_top[0] - control_bottom[0]);
   this->minredco[1] = (y0 - control_bottom[1]) /
                          (control_top[1] - control_bottom[1]);
   this->minredco[2] = (z0 - control_bottom[2]) /
                          (control_top[2] - control_bottom[2]);
   this->maxredco[0] = (x1 - control_bottom[0]) /
                          (control_top[0] - control_bottom[0]);
   this->maxredco[1] = (y1 - control_bottom[1]) /
                          (control_top[1] - control_bottom[1]);
   this->maxredco[2] = (z1 - control_bottom[2]) /
                          (control_top[2] - control_bottom[2]);
}

void CavityEnsemble::setControlVolume(double x0, double y0, double z0, double x1, double y1, double z1)
{
   if((x0 >= x1) || (y0 >= y1) || (z0 >= z1))
   {
      global_log->error() << "\nInvalid control volume (" << x0 << " / " << y0 
                          << " / " << z0 << ") to (" << x1 << " / " << y1 << " / "
                          << z1 << ")." << std::endl;
			exit(711);
   }
	 
   this->restrictedControlVolume = true;
   this->globalV = (x1-x0) * (y1-y0) * (z1-z0);
   this->control_bottom[0] = x0;
   this->control_top[0]    = x1;
   this->control_bottom[1] = y0;
   this->control_top[1]    = y1;
   this->control_bottom[2] = z0;
   this->control_top[2]    = z1;
}

void CavityEnsemble::init(Component* component, unsigned Nx, unsigned Ny, unsigned Nz)
{
   if(this->ownrank < 0)
   {
      global_log->error() << "\nInvalid rank " << ownrank << ".\n";
      exit(712);
   }
   if(this->initialized)
   {
      global_log->error() << "\nCavity ensemble initialized twice.\n";
      exit(713);
   }
   if(0.0 >= this->T)
   {
      global_log->error() << "\nInvalid temperature T = " << T << ".\n";
      exit(714);
   }
   if(0.0 >= this->globalV)
   {
      global_log->error() << "\nInvalid control volume V_ctrl = " << globalV << ".\n";
      exit(715);
   }

   this->componentid = component->ID();
   unsigned rotdof = component->getRotationalDegreesOfFreedom();
   
   double nun[3];
   nun[0] = Nx; nun[1] = Ny; nun[2] = Nz;

   // reduced lattice in each dimension given by epsilon + 0.5/N, epsilon + 1.5/N, ..., epsilon + (N-0.5)/N
   unsigned minlu[3]; unsigned maxlu[3];
   for(int d = 0; d < 3; d++)
   {
      minlu[d] = (unsigned)round(nun[d] * this->minredco[d] - 0.0009756);
      maxlu[d] = (unsigned)round(nun[d] * this->maxredco[d] - 0.0009756) - 1;
   }
   /*
   cout << "nun: " << nun[0] << " / " << nun[1] << " / " << nun[2] << "\n";
   cout << "minredco: " << minredco[0] << " / " << minredco[1] << " / " << minredco[2] << "\n";
   cout << "minlu: " << minlu[0] << " / " << minlu[1] << " / " << minlu[2] << "\n";
   cout << "maxredco: " << maxredco[0] << " / " << maxredco[1] << " / " << maxredco[2] << "\n";
   cout << "maxlu: " << maxlu[0] << " / " << maxlu[1] << " / " << maxlu[2] << "\n";
   */

   double max_spacing = 0.0;
   double grid_spacing[3];
   for(int d = 0; d < 3; d++)
   {
      grid_spacing[d] = (control_top[d] - control_bottom[d])/nun[d];
      if(grid_spacing[d] > max_spacing) max_spacing = grid_spacing[d];
   }

   unsigned tid;
   unsigned tlu[3];
   double tq[3];
   Molecule* tm;
   for(tlu[0] = minlu[0]; maxlu[0] >= tlu[0]; tlu[0]++)
   {
      for(tlu[1] = minlu[1]; maxlu[1] >= tlu[1]; tlu[1]++)
      {
         for(tlu[2] = minlu[2]; maxlu[2] >= tlu[2]; tlu[2]++)
         {
            tid = tlu[0]*Ny*Nz + tlu[1]*Nz + tlu[2] + 1; // + this->idoffset;
            for(int d=0; d < 3; d++)
               tq[d] = control_bottom[d] + (0.5009756 + tlu[d])*grid_spacing[d];

            double v[3];
            double vv = 0.0;
            for(int d=0; d < 3; d++)
            {
               v[d] = -0.5 + this->async.rnd();
               vv += v[d] * v[d];
            }
            double vnorm = sqrt(3.0*T / (vv * component->m()));

            double qtr[4];
            double qqtr = 0.0;
            for(int d=0; d < 4; d++)
            {
               qtr[d] = -0.5 + this->async.rnd();
               qqtr += qtr[d] * qtr[d];
            }
            double qtrnorm = sqrt(1.0 / qqtr);

            double D[3];
            double Dnorm = 0.0;
            if(rotdof > 0)
            {
               for(int d=0; d < 3; d++) D[d] = -0.5 + this->async.rnd();

               double w[3];
               Quaternion tqtr = Quaternion(qtr[0]*qtrnorm, qtr[1]*qtrnorm, qtr[2]*qtrnorm, qtr[3]*qtrnorm);
               tqtr.rotate(D, w);
               double Iw2 = w[0]*w[0]*component->I11()
                          + w[1]*w[1]*component->I22()
                          + w[2]*w[2]*component->I33();

               Dnorm = sqrt(T*rotdof / Iw2);
            }
            else
            {
               D[0] = 0.0; D[1] = 0.0; D[2] = 0.0;
            }

            tm = new Molecule(tid, component, tq[0], tq[1], tq[2], v[0]*vnorm, v[1]*vnorm, v[2]*vnorm, qtr[0]*qtrnorm, qtr[1]*qtrnorm, qtr[2]*qtrnorm, qtr[3]*qtrnorm, D[0]*Dnorm, D[1]*Dnorm, D[2]*Dnorm);

            this->reservoir[tid] = tm;
         }
      }
   }

   this->initialized = true;
}

unsigned long CavityEnsemble::communicateNumCavities(DomainDecompBase* comm)
{
   comm->collCommInit(1);
   comm->collCommAppendUnsLong(this->active.size());
   comm->collCommAllreduceSum();
   this->globalActive = comm->collCommGetUnsLong();
   comm->collCommFinalize();
   
   return this->globalActive;
}

void CavityEnsemble::preprocessStep()
{
   if(this->rotated) return;
   if(this->reservoir.size() == 0) return;
   
   map<unsigned long, Molecule*>::iterator resit = this->reservoir.begin();
   
   double qtr[4];
   Component* tc = resit->second->component();
   unsigned rotdof = tc->getRotationalDegreesOfFreedom();
   if(rotdof > 0)
   {
      for(resit = reservoir.begin(); resit != reservoir.end(); resit++)
      {
            double qqtr = 0.0;
            for(int d=0; d < 4; d++)
            {
               qtr[d] = -0.5 + this->async.rnd();
               qqtr += qtr[d] * qtr[d];
            }
            double qtrnorm = sqrt(1.0 / qqtr);
            resit->second->setq(Quaternion(qtrnorm*qtr[0], qtrnorm*qtr[1], qtrnorm*qtr[2], qtrnorm*qtr[3]));
      }
   }
   this->rotated = true;
}

bool CavityEnsemble::decideActivity(unsigned neighbours, unsigned long tmid)
{
   bool isActive = (this->maxNeighbours >= neighbours);
   if(isActive) this->active.insert(tmid);
   else this->active.erase(tmid);
   return isActive;
}

bool CavityEnsemble::decideActivity(double /*uPotTilde*/, unsigned long tmid)
{
   /*** ADD CAVITY ENSEMBLE ***/
   
   this->active.erase(tmid);
   
   /*** ADD CAVITY ENSEMBLE ***/
   
   return false;
}

map<unsigned long, Molecule*> CavityEnsemble::activeParticleContainer()
{
   map<unsigned long, Molecule*> retv;
   set<unsigned long>::iterator resit;
   for(resit = this->active.begin(); resit != active.end(); resit++)
   {
      retv[*resit] = this->reservoir[*resit];
   }
   return retv;
}
