/*
 * GNU GPL version 2
 */

#include "Domain.h"

#include <fstream>
#include <cstdio>
#include <math.h>
#include <string.h>

const double TC_LJTS = 1.0779;
const double RHOC_LJTS = 0.3190;

using namespace std;

int main(int argc, char** argv) 
{
   const char* usage = "usage: mkesfera <prefix> [-b] [-c] [-e] -I <radius inside> -O <radius outside> [-l <liquid density>] [-v <vapour density>] [-R <cutoff>] [-r] [-S] -T <temperature> [-U] [-u]\n\n-b \tuse bubble type scenario, default: drop type scenario\n-c \tuse be-c-ker format, the default \n-e\tuse B-e-rnreuther format\n-r\tuse b-r-anch format \n-S\tshift, the dafault \n-U\tunshift\n-u\tuse B-u-chholz format\n";
   if((argc < 8) || (argc > 16))
   {
      cout << "There are " << argc
           << " arguments where 8 to 16 should be given.\n\n";
      cout << usage;
      return 1;
   }

   bool in_rhoLiq = false;
   bool in_rhoVap = false;
   
   bool do_shift = true;
   unsigned format = FORMAT_BECKER;
   unsigned type = DROP_TYPE;

   double cutoff = 2.5;
   double rhoLiq;
   double rhoVap;
   double R = 7.0;
   double R2 = 14.0;
   double T = 1.0;

   char* prefix = argv[1];

   for(int i=2; i < argc; i++)
   {
      if(*argv[i] != '-')
      {
         cout << "Flag expected where '" << argv[i]
              << "' was given.\n\n";
         cout << usage;
         return 2;
      }
      for(int j=1; argv[i][j]; j++)
      {
         if(argv[i][j] == 'e') format = FORMAT_BERNREUTHER;
	 else if(argv[i][j] == 'b') type = BUBBLE_TYPE;
	 if(argv[i][j] == 'c') format = FORMAT_BECKER;
         else if(argv[i][j] == 'I')
         {
            i++;
            R = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'l')
         {
            i++;
            rhoLiq = atof(argv[i]);
	    in_rhoLiq = true;
            break;
         }
         else if(argv[i][j] == 'O')
         {
            i++;
            R2 = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'v')
         {
            i++;
            rhoVap = atof(argv[i]);
	    in_rhoVap = true;
            break;
         }
         else if(argv[i][j] == 'R')
         {
            i++;
            cutoff = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'r') format = FORMAT_BRANCH;
         else if(argv[i][j] == 'S') do_shift = true;
         else if(argv[i][j] == 'T')
         {
            i++;
            T = atof(argv[i]);
            break;
         }
         else if(argv[i][j] == 'U') do_shift = false;
         else if(argv[i][j] == 'u') format = FORMAT_BUCHHOLZ;
         else
         {
            cout << "Invalid flag '-" << argv[i][j]
                 << "' was detected.\n\n" << usage;
            return 2; 
         }
      }
   }
   
// default density according to Kedia et al.
if(!in_rhoLiq) rhoLiq = RHOC_LJTS + 0.5649*pow((TC_LJTS - T),(1.0/3.0)) + 0.1314*(TC_LJTS - T) + 0.0413*pow((TC_LJTS - T),(3.0/2.0));
if(!in_rhoVap) rhoVap = RHOC_LJTS - 0.5649*pow((TC_LJTS - T),(1.0/3.0)) + 0.2128*(TC_LJTS - T) + 0.0702*pow((TC_LJTS - T),(3.0/2.0));
if(format == FORMAT_BERNREUTHER)
   {
      cout << "B-e-rnreuther format (flag -e) "
           << "is unavailable at present.\n\n" << usage;
      return 3;
   }

   Domain* dalet;
   if(type == DROP_TYPE){
   dalet = new Domain(R, R2, rhoLiq, rhoVap);
   }
   else{
     dalet = new Domain(R, R2, rhoVap, rhoLiq);
   }
   dalet->write(prefix, cutoff, T, do_shift, format);

   return 0;
}

