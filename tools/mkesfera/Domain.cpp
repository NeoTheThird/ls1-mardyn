/*
 * GNU GPL version 2
 */

#include "Domain.h"
#include "Random.h"
#include <cmath>

#define DT 0.002
#define TIME 20111102
#define VARFRACTION 0.07
#define BOXOVERLOAD 1.3333
#define AVGBIN 0.025

void Domain::write(char* prefix, double cutoff, double T, bool do_shift, int format)
{
   ofstream xdr, txt, buchholz;
   stringstream strstrm, txtstrstrm, buchholzstrstrm;
   if(format == FORMAT_BRANCH)
   {
      strstrm << prefix << ".xdr";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      strstrm << prefix << ".inp";
   }
   if(format == FORMAT_BECKER)
   {
      strstrm << prefix << ".inp";
   }
   xdr.open(strstrm.str().c_str(), ios::trunc);
   if(format == FORMAT_BRANCH)
   {
      txtstrstrm << prefix << "_1R.txt";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      txtstrstrm << prefix << "_1R.cfg";
   }
   if(format == FORMAT_BECKER){
      txtstrstrm << prefix << "_1R.cfg";
   }
   txt.open(txtstrstrm.str().c_str(), ios::trunc);
   if(format == FORMAT_BUCHHOLZ)
   {
      buchholzstrstrm << prefix << "_1R.xml";
      buchholz.open(buchholzstrstrm.str().c_str(), ios::trunc);

      /*
       * Gesamter Inhalt der Buchholz-Datei
       */
      buchholz << "<?xml version = \'1.0\' encoding = \'UTF-8\'?>\n<mardyn version=\""
               << TIME << "\">\n   <simulation type=\"MD\">\n      <input type=\"oldstyle\">"
               << prefix << "_1R.cfg</input>\n   </simulation>\n</mardyn>";
   }

   unsigned fl_units;
   double rhomax = (rho_i > rho_o)? rho_i: rho_o;
   double N_boxes = 8.0*R_o*R_o*R_o * (BOXOVERLOAD*rhomax) / 3.0;
   fl_units = ceil(pow(N_boxes, 1.0/3.0));
   double fl_unit = 2.0*R_o / (double)fl_units;

   Random* r = new Random();
   r->init(
      (int)(10000.0*R_o) - (int)(3162.3*cutoff)
        + (int)(1000.0*T) - (int)(316.23*rho_i)
        + (int)(100.0*R_i)
   );

   bool fill[fl_units][fl_units][fl_units][3];
   unsigned slots = 3.0*fl_units*fl_units*fl_units;
   double boxdensity = (double)slots / (8.0*R_o*R_o*R_o);
   cerr << "Box density: " << boxdensity << " (unit cell: " << fl_unit << ").\n";
   double P_in = rho_i / boxdensity;
   double P_out = rho_o / boxdensity;
   cerr << "Insertion probability: " << P_in << " inside, " << P_out << " outside.\n";

   double goffset[3][3];
   goffset[0][0] = 0.0; goffset[1][0] = 0.5; goffset[2][0] = 0.5;
   goffset[0][1] = 0.5; goffset[1][1] = 0.0; goffset[2][1] = 0.5;
   goffset[0][2] = 0.5; goffset[1][2] = 0.5; goffset[2][2] = 0.0;

   unsigned N = 0;
   unsigned idx[3];
   for(idx[0]=0; idx[0] < fl_units; idx[0]++)
      for(idx[1]=0; idx[1] < fl_units; idx[1]++)
         for(idx[2]=0; idx[2] < fl_units; idx[2]++)
            for(int p=0; p < 3; p++)
            {
               double qq = 0.0;
               for(int d = 0; d < 3; d++)
               {
                  double qrel = (idx[d] + goffset[d][p])*fl_unit - R_o;
                  qq += qrel*qrel;
               }
               double tP = (qq > R_i*R_i)? P_out: P_in;
               bool tfill = (tP >= r->rnd());
               fill[idx[0]][idx[1]][idx[2]][p] = tfill;
               if(tfill) N++;
            }
   cerr << "Filling " << N << " out of " << slots << " slots (total density: " << N / (8.0*R_o*R_o*R_o) << ").\n";

   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt.precision(6);
      txt << "mardynconfig\n# \ntimestepLength\t" << DT << "\ncutoffRadius\t" << cutoff << "\nLJCutoffRadius\t" << cutoff << "\ntersoffCutoffRadius\t0.5\n";
   }
   if(format == FORMAT_BECKER){
     txt.precision(6);

	cout << "\n**********************************\nWriting the config file \n**********************************\n\n";

	txt << "MDProjectConfig\n# \ntimestepLength\t" << DT << "\ncutoffRadius\t"<< cutoff << "\nLJCutoffRadius\t"<< cutoff << "\ntersoffCutoffRadius\t0.5\n";
	txt << "initCanonical\t5000\ninitStatistics\t100000\nphaseSpaceFile\tOldStyle\t"<< prefix<<".inp\n";
	txt << "parallelization DomainDecomposition \n";
	//txt << "parallelization KDDecomposition2 200 3 \n";
	txt << "# for LinkedCells, the cellsInCutoffRadius has to be provided\ndatastructure\tLinkedCells \t1\n";
	txt << "output\tResultWriter\t"<< 40 <<"\t"<< prefix<<"\noutput\tXyzWriter\t"<< 100000<< "\t"<< prefix<<".buxyz\n";
	txt << "profile\t1 " << (unsigned)round(R_o / AVGBIN) << " 1\nprofileRecordingTimesteps\t1\nprofileOutputTimesteps\t40\n";
	txt << "profiledComponent\t1\nprofileOutputPrefix\t"<<prefix<<"\n";
	txt << "AlignCentre\t25\t1\n" ;
	txt << "nomomentum\t200\n";
//	cout << "\n**********************************\nwrite() method of Configwriter finished\n**********************************\n";

   }

   if(format == FORMAT_BRANCH) txt << "phaseSpaceFile\t" << prefix << ".xdr\n";
   if(format == FORMAT_BUCHHOLZ) txt << "phaseSpaceFile\tOldStyle\t" << prefix << ".inp\nparallelization\tDomainDecomposition\n";

   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      txt << "datastructure\tLinkedCells 1\n";

      txt << "output\tResultWriter 500\t" << prefix << "_1R\nresultOutputTimesteps\t500\noutput\tXyzWriter 60000\t" << prefix << "_1R\noutput\tXdrWriter 60000\t" << prefix << "_1R\ninitCanonical\t10\ninitStatistics\t3003003\n";
      txt << "profile\t1 " << (unsigned)round(R_o / AVGBIN) << " 1\nprofileRecordingTimesteps\t1\nprofileOutputTimesteps\t500000\nprofiledComponent\t1\nprofileOutputPrefix\t" << prefix << "_1R\nesfera\nnomomentum 16384\nAlignCentre 333 0.003\nchemicalPotential 0 component 1 conduct " << (unsigned)round(N/3.0) << " tests every 1 steps\nWidom\n";
   }

   xdr.precision(8);
   if(format == FORMAT_BECKER){
      xdr << "mardyn trunk "<< TIME << "\n#mardyn input file, ls1 project\n#generated by the mkesfera tool\n";
      xdr << "t\t0.0\n";
      xdr << "T \t"<< T << "\n";
      xdr << "L\t" << 2.0*R_o << " " << 2.0*R_o << " " << 2.0*R_o << "\n";
      xdr << "C\t1 \n1 0 0 0 0\n0.0 0.0 0.0 \t1 1 1 " << cutoff << " " << (do_shift? " 1": " 0")  << "\t0.0 0.0 0.0\n1e+10\n";
      xdr << "N\t" << N << "\nM\tICRVQD\n"; 
   }
   if(format == FORMAT_BRANCH)
   {
      xdr << "mardyn " << TIME << " tersoff\n"
          << "# mardyn input file, ls1 project\n"
          << "# generated by the mkesfera tool\n# \n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      xdr << "mardyn trunk " << TIME << "\n";
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      xdr << "t\t0\nL\t" << 2.0*R_o << " " << 2.0*R_o << " " << 2.0*R_o << "\n";
   }
   if(format == FORMAT_BUCHHOLZ)
   {
      xdr << "T\t" << T << "\n";
      xdr << "C\t1\t1 0 0 0 0\t0 0 0 1 1 1\t" << cutoff << " " << (do_shift? " 1": " 0") << "\t0 0 0 1e+10\n";
   }
   if(format == FORMAT_BRANCH)
   {
      xdr << "C\t1\t1 0 0 0 0\t0 0 0 1 1 1\t0 0 0 1e+10\n";
      xdr << "T\t" << T << "\n";
   }
   if((format == FORMAT_BRANCH) || (format == FORMAT_BUCHHOLZ))
   {
      xdr << "N\t" << N << "\nM\tICRVQD\n";
   }

   double v = sqrt(3.0 * T);

   unsigned ID = 1;
   for(idx[0]=0; idx[0] < fl_units; idx[0]++)
      for(idx[1]=0; idx[1] < fl_units; idx[1]++)
         for(idx[2]=0; idx[2] < fl_units; idx[2]++)
            for(unsigned p=0; p < 3; p++)
            {
               // cerr << idx[0] << "\t" << idx[1] << "\t" << idx[2] << "\t|\t" << fill[idx[0]][idx[1]][idx[2]][p] << "\n";
               if(fill[idx[0]][idx[1]][idx[2]][p])
               {
                  double q[3];
                  for(int d=0; d < 3; d++)
                  {
                     q[d] = (idx[d] + VARFRACTION*(r->rnd() - 0.5) + goffset[d][p])*fl_unit;
                     if(q[d] < 0.0) q[d] += 2.0*R_o;
                     else if(q[d] > 2.0*R_o) q[d] -= 2.0*R_o;
                  }
                  double phi = 6.283185 * r->rnd();
                  double omega = 6.283185 * r->rnd();

                  xdr << ID << " " << 1 << "\t" << q[0]
                      << " " << q[1] << " " << q[2]
                      << "\t" << v*cos(phi)*cos(omega) << " "
                      << v*cos(phi)*sin(omega) << " "
                      << v*sin(phi) << "\t1 0 0 0 0 0 0\n";

                  ID++;
               }
            }

   xdr.close();
   txt.close();
   if(format == FORMAT_BUCHHOLZ)
   {
      buchholz.close();
   }
}

