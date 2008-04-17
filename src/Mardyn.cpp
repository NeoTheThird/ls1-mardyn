#include <iostream>
#include <sstream>
#include <string>

#include "TestCaseRunner.h"
#include "Simulation.h"
#include "utils/Log.h"
//#include "parallel/DomainDecompBase.h"
using namespace std;

//! @page main
//! In this project, software for the simulation of Molecular Dynamics
//! with short-range forces is developed. The aim is to have a parallel code (MPI and OpenMP) 
//! for multi-centered molecules.
//!
//! The role of the main function is to run Tests for all classes
//! and to instantiate an Object of the Simulation class which
//! then is responsible for the simulation
//!
int main(int argc, char** argv){
  // Logging interface
  static utils::Log _log("mardyn");
  stringstream logtext;
  
#ifndef NDEBUG
  //Test all modules
  // TestCaseRunner runner;
  // runner.run();
  cout << " <> ";
#endif
  Simulation simulation(&argc, &argv);
  
  cout << "Simulation object created.\n";
  double runtime = double(clock())/CLOCKS_PER_SEC;

  simulation.simulate();

  runtime=double(clock())/CLOCKS_PER_SEC-runtime;

  logtext << "used " << runtime << " s" << endl;
  _log.info("main", logtext.str());
  logtext.str("");

}
