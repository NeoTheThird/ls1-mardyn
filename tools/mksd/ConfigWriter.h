/*
 * ConfigWriter.h
 *
 *  Created on: 05.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 */

#ifndef CONFIGWRITER_H_
#define CONFIGWRITER_H_


#include <cstring>
#include <iostream>
#include <fstream> // Schreiben in / lesen aus Datei
#include <sstream>


using namespace std;

class ConfigWriter{

public:
	ConfigWriter(
			char* prefix, string in_wall, int in_wallLays, double in_sigFluid,
			double in_refTime, double in_cutoffRadius, double in_ljCutoffRadius, double in_wallCutoffRadius,
			unsigned in_profilePhi, unsigned in_profileR, unsigned in_profile_H,
			unsigned in_profileOutputTimesteps, bool movie
			);
	~ConfigWriter();

	//void buildString();
	void write();

	// set and write methods for the private attributes
	void	sPrefix(char* input);
	void 	sTimestepLength(double input);
	void	sCutoffRadius(double input);
	void	sLjCutoffRadius(double input);
	void	sTersoffCutoffRadius(double input);
	void	sLjWallCutoffRadius(double input);
	void	sInitCanonical(unsigned input);
	void	sInitStatistics(unsigned input);
	void	sProfileOutputTimesteps(unsigned input);
	void	sProfileRecordingTimesteps(unsigned input);
	void	sProfile (unsigned profilePhi, unsigned profileR, unsigned profileH);
	void	sHProfile(unsigned profileH);
	void 	sOutputResWriter(unsigned input);
	void	sOutputXyzWriter(unsigned input);
	void	sOutputVisittWriter(unsigned input);
	void	sOutputMmspdWriter(unsigned input);

	// @brief: each writing-method immedeately writes one line
	// of the config file in the output-stream confStrm
	char*	wPrefix();
	void 	wTimestepLength();
	void 	wCutoffRadius();
	void 	wLjCutoffRadius();
	void	wTersoffCutoffRadius();
	void	wLjWallCutoffRadius();
	void 	wInitCanonical();
	void 	wInitStatistics();
	void 	wProfileOutputTimesteps();
	void 	wProfileRecordingTimesteps();
	void 	wProfile ();
	void	wPhaseSpaceFileName();
	void	wOutputResWriter();
	void	wOutputXyzWriter();
	void	wOutputVisittWriter();
	void 	wOutputMmspdWriter();
	void	wProfileOutputPrefix();

	/*
	// methods to concatenate multiple C-strings (char*):
	//based on the repeated application of the <cstring> function char* strcat(char* s1,const char* s2)
	char* merge(char* s1, char* s2, char* s3);
	char* merge(char* s1, char* s2, char* s3, char*s4);
	char* merge(char* s1, char* s2, char* s3, char*s4, char* s5);
	*/

	//timestepLength(double refTime)
	//cutoffLJ(double rcLjFaktor);

private:
	double timestepLength;
	double cutoffRadius;
	double ljCutoffRadius;
	double ljWallCutoffRadius;
	double tersoffCutoffRadius;
	double sigFluid;

	string wall; // depicts the wall model

	unsigned initCanonical;
	unsigned initStatistics;
	unsigned profileOutputTimesteps;
	unsigned profileRecordingTimesteps;
	unsigned outputResWriter;
	unsigned outputXyzWriter;
	unsigned outputVisittWriter;
	unsigned outputMmspdWriter;

	unsigned profile[3]; // number of profile units in the order phi-r-h => density profile in cylindrical coordinates

	int wallLays;

	char* prefix;				// ending _1R added
//	char* phaseSpacePrefix;		// prefix as passed (without _1R)

	bool _movie;

	// declaration of the output streams
	ofstream confStrm;			// Ausgabe-Stream für programminterne Nutzung: wohin C++ primär schreibt
	stringstream confFile;		// "physikalisch" vorhandene Datei prefix.cfg
	//char* outputString;	// used in case the entire data are written as a single string in the config file
};

#endif /* CONFIGWRITER_H_ */
