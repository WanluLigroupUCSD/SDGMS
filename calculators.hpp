#ifndef calculators_hpp
#define calculators_hpp

#include <string>
#include "chemistry.hpp"

#ifndef DBL_MAX
#define DBL_MAX 1.79769e+308;
#endif

int getBasisSetElectrons(std::string element, int complexity);
std::string getBasisSet(std::string element, int complexity);

void writeToCP2K(structure& s, std::string task, std::string filename, std::string basis_file, std::string potential_file, double a = 0, double b = 0, double c = 0);
void writeToGuassian(structure& s, std::string task, std::string filename, std::string basis_name, std::string method, int charge = 0, int state = 1, int maxCycles = 0, int maxSCF = 129, std::string customBasis = "empty");
void writeToADF(structure& s, std::string task, std::string filename, std::string basis_name, std::string method, int maxSCF = 300, int charge = 0, int state = 1, std::string converge = "1.0e-6", bool unrestricted = false);
structure* optimizationADF(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, std::string AMSHOME, int amscalculation = -1, int maxSCF = 300, std::string converge = "1.0e-6", double timeoutminutes = -1, double timeLeftMinutes = -1, bool keepFiles = false, bool unrestricted = false);
double energyADF(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, std::string AMSHOME, int adfcalculation = -1, std::string converge = "1.0e-6", int timeoutminutes = -1, int timeLeftMinutes = -1, bool keepFiles = false, bool unrestricted = false);
double basinHoppingEnergy(structure& s, int ittt, std::string task, double a, double b, double c);
structure* optimizationBasinHopping(structure& s, int ittt, std::string task, int maxCycles, double x, double y, double z);
double energyCP2K(structure& s, int ittt, std::string task, double a, double b, double c);
double energyGuassian(structure& s, int ittt, std::string task, std::string basis, std::string method, int charge = 0, int state = 0, int maxSCF = 129, bool keepFiles = false, std::string customBasis = "empty");
structure* optimizationGaussian(structure& s, int ittt, std::string task, std::string basis, std::string method, int maxCycles, int charge = 0, int state = 0, int maxSCF = 129, bool keepFiles = false, std::string customBasis = "empty");
structure* newOptimizationGaussian(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, int amscalculation = -1, int maxSCF = 300, int maxCycles = 100, double timeoutminutes = -1, double timeLeftMinutes = -1, bool keepFiles = false, std::string customBasis = "empty");
double newEnergyGaussian(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, int amscalculation = -1, int maxSCF = 300, int maxCycles = 100, double timeoutminutes = -1, double timeLeftMinutes = -1, bool keepFiles = false, std::string customBasis = "empty");

std::vector<double> direction(structure& a, structure& b);



//double energyCalculation(structure* s, std::string basis, std::string method, int charge, int state, char calculator, std::string taskName, int energyCall, bool keepFiles, std::string converge = "1.0e-6", double timeLeftMinutes = -1, double timeoutMinutes = -1, bool unrestricted = false, std::string customBasis = "empty", std::string AMSHOME = "/home/jpburkhardt/scm/ams2023.104");
//structure* optimizationCalculation(char calculator, structure* Xnew, int step, std::string taskName, int optimizationIterations, int maxSCF, double x, double y, double z, std::vector<double>& optimizationMovement, std::vector<structure*>& structures, std::string basis, std::string method, int charge, int state, bool keepFiles, std::string converge = "1.0e-6", double timeLeftMinutes = -1, double timeoutMinutes = -1, bool unrestricted = false, std::string customBasis = "empty", std::string AMSHOME = "/home/jpburkhardt/scm/ams2023.104");

/*
THese use premade DFT files
*/
structure* optimization(std::string calculator, structure* x, std::string dftFileName,std::string cancelFile, double timeoutminutes, double timeLeftMinutes, std::string xyzfilename, int calculationNo);
double energyDFT(std::string calculator, structure* x, std::string dftFileName,std::string cancelFile, double timeoutminutes, double timeLeftMinutes, std::string xyzfilename, int calculationNo);

structure* readOptFile(std::string calculator, std::string dftFileName, int calculationNo, std::vector<int> composition, std::vector<std::string> elements);
#endif