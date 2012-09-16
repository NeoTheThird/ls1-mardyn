/*
 * ConfigWriter.cpp
 *
 *  Created on: 05.01.2012
 *      Author: becker
 */


#include"ConfigWriter.h"

const double DT = 0.030620; // corresponds to 1 fs


extern const string WALL_TERSOFF;
extern const string WALL_CU_LJ;
extern double LATTICE_CONST_WALL_LJTS;
//extern const double  LATTICE_CONST_CU;


// @brief: implementing the constructor and destructor, respectively
ConfigWriter::ConfigWriter(
		char* in_prefix, string in_wall, int in_wallLays, double in_sigFluid,
		double in_refTime, double in_cutoffRadius, double in_ljCutoffRadius, double in_wallCutoffRadius,
		unsigned in_profilePhi, unsigned in_profileR, unsigned in_profile_H,
		unsigned in_profileOutputTimesteps, unsigned initCanon, bool in_movie
		)
{
	cout << "\n**********************************\nConfigwriter opened\n**********************************\n";
	sPrefix(in_prefix);
	sTimestepLength(in_refTime);
	wallLays = in_wallLays;
	sigFluid = in_sigFluid;
	//unknown: what kind of cutoff radius???
	sCutoffRadius(in_cutoffRadius);
	sLjCutoffRadius(in_ljCutoffRadius);
	// the only information needed in the config file:
	// what wall model applied => writing out the corresponding cutoff radius
	if(in_wall == WALL_TERSOFF)
	{
		wall = WALL_TERSOFF;
		sTersoffCutoffRadius(in_wallCutoffRadius);	// in case the wall is modeled by the Tersoff potential
	}
	//@brief: obsolete!
	/*else if (in_wall == WALL_CU_LJ)
	{
		wall = WALL_CU_LJ;
		sLjWallCutoffRadius(in_wallCutoffRadius);
	}*/
	// the error message following should not be necessary currently,
	//implemented for error detection and further models implemented
	//@todo: if no error reported during test runs => removing the else-branch
	else
	{
		cout << "wall model: only Lennard-Jones TS.\n";
		//return 51;
	}
	sProfile (in_profilePhi, in_profileR, in_profile_H);
	sProfileOutputTimesteps(in_profileOutputTimesteps);
	// values set as default, if required they can be passed by main.cpp, i.e. as a user input => additional constructor
	sInitCanonical(initCanon);					// temperature raise during equilibration applied? (ask M.Horsch) => short equilibration time needed otherwise vaporisation of the drop
	sInitStatistics(500000);
	sProfileRecordingTimesteps(2);
	sOutputResWriter(40);
	sOutputXyzWriter(500000);
	// if no movie is to be made, there's no need for the visitt-writer, by default: movie = false; then also no VisItt output is generated: if(movie) wOutputVisittWriter();
	_movie = in_movie;
	if(_movie){
		sOutputVisittWriter(500);
		sOutputMmspdWriter(500);
	}
	//	cout << "\n**********************************\nConstructor of Configwriter finished\n**********************************\n";
}

ConfigWriter::~ConfigWriter()
{
	confStrm.close();
	//cout << "Config file completed.\n";
}

//@brief: superior writing method: handling the streams, calling the single write-methods
void ConfigWriter::write(){

//	cout << "\n**********************************\nwrite() method of Configwriter started\n**********************************\n";
	confFile << prefix << "_1R.cfg";					// building the file prefix.cfg
	confStrm.open(confFile.str().c_str(), ios::trunc);	// linking the file prefix.cfg with confStrm => confStrm writes in prefix.cfg

	cout << "\n**********************************\nWriting the config file \n**********************************\n\n";

	confStrm << "MDProjectConfig\n";
	wTimestepLength();
	wCutoffRadius();
	wLjCutoffRadius();
	confStrm <<"\n";
	// set the appropriate cutoff radius, the appropriate wall model, respectively
	//if(wall == WALL_CU_LJ)	wLjWallCutoffRadius();
	//else if (wall == WALL_TERSOFF) wTersoffCutoffRadius();
	wInitCanonical();
	wInitStatistics();
	wPhaseSpaceFileName();
	confStrm << "parallelization DomainDecomposition \n";
	//confStrm << "parallelization KDDecomposition2 200 3 \n";
	confStrm << "# for LinkedCells, the cellsInCutoffRadius has to be provided\n"
			 << "datastructure\tLinkedCells \t1\n";
	wOutputResWriter();
	wOutputXyzWriter();
	// only if a movie is to be made
	if(_movie){
	  //wOutputVisittWriter();
	  wOutputMmspdWriter();
	}
	wProfile();
	confStrm << "yOffset\t" << (wallLays -0.5 + 0.05)*LATTICE_CONST_WALL_LJTS << "\n"; // +0.05*latticeConstant => little offset to avoid particles beeing placed outside the simulation box
	wProfileRecordingTimesteps();
	wProfileOutputTimesteps();
	confStrm <<"profiledComponent\t1\n";
	confStrm << "SessileDrop\n";
	wProfileOutputPrefix();
	confStrm << "AlignCentre\t25\t1\n" ;
	confStrm << "ComponentForYShift\t2 3\n";
	confStrm << "nomomentum\t200\n";
//	cout << "\n**********************************\nwrite() method of Configwriter finished\n**********************************\n";

}


//@brief: implemeting the set-method s...(arg)
	void	ConfigWriter::sPrefix(char* input)
	{
		prefix = input;
	}
	void 	ConfigWriter::sTimestepLength(double input)
	{
		timestepLength = DT/input;
	}
	void	ConfigWriter::sCutoffRadius(double input)
	{
		cutoffRadius = input;
	}
	void	ConfigWriter::sLjCutoffRadius(double input)
	{
		ljCutoffRadius = input;
	}
	void	ConfigWriter::sTersoffCutoffRadius(double input)
	{
		tersoffCutoffRadius = input;
	}
	void	ConfigWriter::sLjWallCutoffRadius(double input)
		{
			ljWallCutoffRadius = input;
		}
	void	ConfigWriter::sInitCanonical(unsigned input)
	{
		initCanonical = input;
	}
	void	ConfigWriter::sInitStatistics(unsigned input)
	{
		initStatistics = input;
	}
	void	ConfigWriter::sProfileOutputTimesteps(unsigned input)
	{
		profileOutputTimesteps = input;
	}
	void	ConfigWriter::sProfileRecordingTimesteps(unsigned input)
	{
		profileRecordingTimesteps = input;
	}
	void	ConfigWriter::sProfile (unsigned profilePhi, unsigned profileR, unsigned profileH)
	{
		profile[0] = profilePhi;
		profile[1] = profileR;
		profile[2] = profileH;
	}
	void 	ConfigWriter::sHProfile(unsigned profileH){
		profile[2] = profileH;
	}
	void 	ConfigWriter::sOutputResWriter(unsigned input)
			{
				outputResWriter = input;
			}
	void	ConfigWriter::sOutputXyzWriter(unsigned input)
	{
		outputXyzWriter = input;
	}
	void	ConfigWriter::sOutputVisittWriter(unsigned input)
	{
		outputVisittWriter = input;
	}
	void 	ConfigWriter::sOutputMmspdWriter(unsigned input)
	{
		outputMmspdWriter = input;
	}



//@brief: implementing the write-methods w...()
	char*	ConfigWriter::wPrefix()
		{
		return prefix;
		}
	void 	ConfigWriter::wTimestepLength()
		{
		confStrm << "timestepLength\t"<<timestepLength<<"\n";
		}
	void 	ConfigWriter::wCutoffRadius()
		{
		confStrm << "cutoffRadius\t"<< cutoffRadius<<"\n";
		}
	void 	ConfigWriter::wLjCutoffRadius()
		{
		confStrm<<"LJCutoffRadius\t"<< ljCutoffRadius<<"\n";
		}
	void	ConfigWriter::wTersoffCutoffRadius()
		{
		confStrm<<"tersoffCutoffRadius\t"<< tersoffCutoffRadius<<"\n";
		}
	void	ConfigWriter::wLjWallCutoffRadius()
			{
			confStrm<<"ljWallCutoffRadius\t"<< ljWallCutoffRadius<<"\n";
			}
	void 	ConfigWriter::wInitCanonical()
		{
		confStrm<<"initCanonical\t"<< initCanonical<<"\n";
		}
	void 	ConfigWriter::wInitStatistics()
		{
		confStrm<<"initStatistics\t"<< initStatistics<<"\n";
		}
	void 	ConfigWriter::wProfileOutputTimesteps()
		{
		confStrm<<"profileOutputTimesteps\t"<< profileOutputTimesteps<<"\n";
		}
	void 	ConfigWriter::wProfileRecordingTimesteps()
		{
		confStrm<<"profileRecordingTimesteps\t"<< profileRecordingTimesteps<<"\n";
		}
	void 	ConfigWriter::wProfile ()
		{
		confStrm<<"profile\t"<< profile[0]<<"\t"<<profile[1]<<"\t"<<profile[2]<<"\n";
		}
	void	ConfigWriter::wPhaseSpaceFileName()
		{
		confStrm<<"phaseSpaceFile\t"<< "OldStyle\t"<< prefix<<".inp\n";
		}
	void	ConfigWriter::wOutputResWriter()
		{
		confStrm<<"output\tResultWriter\t"<< outputResWriter<<"\t"<< wPrefix()<<"\n";
		}
	void	ConfigWriter::wOutputXyzWriter()
		{
		confStrm<<"output\tXyzWriter\t"<< outputXyzWriter<< "\t"<< wPrefix()<<".buxyz\n";
		}
	void	ConfigWriter::wOutputVisittWriter()
		{
		confStrm<<"output\tVisittWriter\t"<< outputVisittWriter<<"\t"<<wPrefix()<<"\n";
		}
	void	ConfigWriter::wOutputMmspdWriter()
		{
		confStrm<<"output\tMmspdWriter\t"<< outputMmspdWriter<<"\t"<<wPrefix()<<"\n";
		}
	void	ConfigWriter::wProfileOutputPrefix()
		{
		confStrm<<"profileOutputPrefix\t"<< wPrefix()<<"\n";
		}



/*
// implementing the merge(args)-functions
	char* merge(char* s1, char* s2, char* s3)
	{
		char* merged;
		merged = strcat(s1, s2);
		merged = strcat(merged, s3);
		return merged;
	}
	char* merge(char* s1, char* s2, char* s3, char*s4)
	{
		char* merged;
		merged = strcat(s1, s2);
		merged = strcat(merged, s3);
		merged = strcat(merged, s4);
	}
	char* merge(char* s1, char* s2, char* s3, char*s4, char* s5)
	{
		char* merged;
		merged = strcat(s1, s2);
		merged = strcat(merged, s3);
		merged = strcat(merged, s4);
		merged = strcat(merged, s5);
	}
*/

