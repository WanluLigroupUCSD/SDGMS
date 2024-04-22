// ModifiedBasinHopping.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Written by Jordan Burkhardt for the Wan-Lu Li lab at UCSD

#include <iostream>
#include "chemistry.hpp"
#include "DistanceHeuristics.hpp"
#include <functional>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdlib.h>

const int DEBUG = 0;
const int STDOUT = 1;

double temperature(int step,int averageSteps,double eVLimit)
{
	/*
	at a difference of - 1.5 eV //from Wan-Lu Li paper
	at 25000 Kelvins the probability becomes 0.5
	we should reach this exponentially, but after how many steps?

	*/
	/*
	eVLimit: the difference between seed and new energy where PES signal-tonoise ratios become poor
	averageSteps: the number of steps required to reach a probability of 50% acceptance
	step: current step
	*/
	//determine the equation:
	double halfValue = -std::log(0.5);
	double temp = eVLimit / (halfValue * boltzman);
	double base = std::pow(temp, 1.0 / averageSteps);
	return 0.8 + std::pow(base, step);
}

std::vector<std::set<int>> generateAcceptedCoordinationNumbers(std::vector<structure*> s)
{
	std::vector < std::set<int>> returnValue;
	//one element at a time
	for (int e = 0; e < s[0]->composition.size(); e++)
	{
		//for each element
		std::set<int> acceptedCoordinationNumbers;
		//for each structure
		for (std::vector<structure*>::iterator structure = s.begin(); structure != s.end(); structure++)
		{
			//for each coordination number in the strucutre of the current element
			for (int a = 0; a < (*structure)->coordinationNumbers[e].size(); a++)
			{
				acceptedCoordinationNumbers.insert((*structure)->coordinationNumbers[e][a]);
			}
		}
		returnValue.push_back(acceptedCoordinationNumbers);

	}
	return returnValue;

}

bool coordinationNumber(structure& s, std::vector<std::set<int>> coordinationNumbers)
{
	for (int i = 0; i < s.coordinationNumbers.size(); i++)
	{
		//for each element
		for (int a = 0; a < s.coordinationNumbers[i].size(); a++)
		{
			//for each atom
			bool accept = false;
			for (int b = 0; b < coordinationNumbers[i].size(); b++)
			{
				//for each accetable coordination number for element i
				if (a == b) accept = true;
			}
			if (!accept) return false;
		}
	}
	return true;
}

int getBasisSetElectrons(std::string element, int complexity)
{
		if (element == "H") {
			return 1;
		}
		else if (element == "He") {
			return 2;
		}
		else if (element == "Li") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 3;
		}
		else if (element == "Be") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 4;
		}
		else if (element == "B") {
			return 3;
		}
		else if (element == "C") {
			return 4;
		}
		else if (element == "N") {
			return 5;
		}
		else if (element == "O") {
			return 6;
		}
		else if (element == "F") {
			return 7;
		}
		else if (element == "Ne") {
			return 8;
		}
		else if (element == "Na") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 9;
		}
		else if (element == "Mg") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 10;
		}
		else if (element == "Al") {
			return 3;
		}
		else if (element == "Si") {
			return 4;
		}
		else if (element == "P") {
			return 5;
		}
		else if (element == "S") {
			return 6;
		}
		else if (element == "Cl") {
			return 7;
		}
		else if (element == "Ar") {
			return 8;
		}
		else if (element == "K") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 9;
		}
		else if (element == "Ca") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 10;
		}
		else if (element == "Sc") {
			if (complexity == 1) return 3;
			else if (complexity == 2) return 11;
		}
		else if (element == "Ti") {
			if (complexity == 1) return 4;
			else if (complexity == 2) return 12;
		}
		else if (element == "V") {
			if (complexity == 1) return 5;
			else if (complexity == 2) return 13;
		}
		else if (element == "Cr") {
			if (complexity == 1) return 6;
			else if (complexity == 2) return 14;
		}
		else if (element == "Mn") {
			if (complexity == 1) return 7;
			else if (complexity == 2) return 15;
		}
		else if (element == "Fe") {
			if (complexity == 1) return 8;
			else if (complexity == 2) return 16;
		}
		else if (element == "Co") {
			if (complexity == 1) return 9;
			else if (complexity == 2) return 17;
		}
		else if (element == "Ni") {
			if (complexity == 1) return 10;
			else if (complexity == 2) return 18;
		}
		else if (element == "Cu") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 11;
			else if (complexity == 3) return 19;
		}
		else if (element == "Zn") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 12;
			else if (complexity == 3) return 20;
		}
		else if (element == "Ga") {
			if (complexity == 1) return 3;
			else if (complexity == 2) return 13;
		}
		else if (element == "Ge") {
			return 4;
		}
		else if (element == "As") {
			return 5;
		}
		else if (element == "Se") {
			return 6;
		}
		else if (element == "Br") {
			return 7;
		}
		else if (element == "Kr") {
			return 8;
		}
		else if (element == "Rb") {
			return 1;
		}
		else if (element == "Sr") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 10;
		}
		else if (element == "Y") {
			if (complexity == 1) return 3;
			else if (complexity == 2) return 11;
		}
		else if (element == "Zr") {
			if (complexity == 1) return 4;
			else if (complexity == 2) return 12;
		}
		else if (element == "Nb") {
			if (complexity == 1) return 5;
			else if (complexity == 2) return 13;
		}
		else if (element == "Mo") {
			if (complexity == 1) return 6;
			else if (complexity == 2) return 14;
		}
		else if (element == "Tc") {
			if (complexity == 1) return 7;
			else if (complexity == 2) return 15;
		}
		else if (element == "Ru") {
			if (complexity == 1) return 8;
			else if (complexity == 2) return 16;
		}
		else if (element == "Rh") {
			if (complexity == 1) return 9;
			else if (complexity == 2) return 17;
		}
		else if (element == "Pd") {
			if (complexity == 1) return 10;
			else if (complexity == 2) return 18;
		}
		else if (element == "Ag") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 11;
			else if (complexity == 3) return 19;
		}
		else if (element == "Cd") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 12;
		}
		else if (element == "In") {
			if (complexity == 1) return 3;
			else if (complexity == 2) return 13;
			else if (complexity == 3) return 146;
		}
		else if (element == "Sn") {
			return 4;
		}
		else if (element == "Sb") {
			return 5;
		}
		else if (element == "Te") {
			return 6;
		}
		else if (element == "I") {
			if (complexity == 1) return 7;
			else if (complexity == 2) return 129;
			else if (complexity == 3) return 125;
		}
		else if (element == "Xe") {
			return 8;
		}
		else if (element == "Cs") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 9;
		}
		else if (element == "Ba") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 10;
		}
		else if (element == "La") {
			return 11;
		}
		else if (element == "Ce") {
			return 12;
		}
		else if (element == "Pr") {
			return 13;
		}
		else if (element == "Nd") {
			return 14;
		}
		else if (element == "Pm") {
			return 15;
		}
		else if (element == "Sm") {
			return 16;
		}
		else if (element == "Eu") {
			return 17;
		}
		else if (element == "Gd") {
			return 18;
		}
		else if (element == "Tb") {
			return 19;
		}
		else if (element == "Dy") {
			return 20;
		}
		else if (element == "Ho") {
			return 21;
		}
		else if (element == "Er") {
			return 22;
		}
		else if (element == "Tm") {
			return 23;
		}
		else if (element == "Yd") {
			return 24;
		}
		else if (element == "Lu") {
			return 25;
		}
		else if (element == "Hf") {
			return 12;
		}
		else if (element == "Ta") {
			if (complexity == 1) return 5;
			else if (complexity == 2) return 13;
		}
		else if (element == "W") {
			if (complexity == 1) return 6;
			else if (complexity == 2) return 14;
		}
		else if (element == "Re") {
			if (complexity == 1) return 7;
			else if (complexity == 2) return 15;
		}
		else if (element == "Os") {
			if (complexity == 1) return 8;
			else if (complexity == 2) return 16;
		}
		else if (element == "Ir") {
			if (complexity == 1) return 9;
			else if (complexity == 2) return 17;
		}
		else if (element == "Pt") {
			if (complexity == 1) return 10;
			else if (complexity == 2) return 18;
		}
		else if (element == "Au") {
			if (complexity == 1) return 1;
			else if (complexity == 2) return 11;
			else if (complexity == 3) return 19;
		}
		else if (element == "Hg") {
			if (complexity == 1) return 2;
			else if (complexity == 2) return 12;
		}
		else if (element == "Tl") {
			if (complexity == 1) return 3;
			else if (complexity == 2) return 13;
		}
		else if (element == "Pb") {
			return 4;
		}
		else if (element == "Bi") {
			return 5;
		}
		else if (element == "Po") {
			return 6;
		}
		else if (element == "At") {
			return 7;
		}
		else if (element == "Rn") {
			return 8;
		}
		else if (element == "Fr") {
			//no basis?
		}
		else if (element == "Ra") {
			//no basis?
		}
		else if (element == "Ac") {
			if (complexity == 1) return 11;
			else if (complexity == 2) return 29;
		}
		else if (element == "Th") {
			if (complexity == 1) return 12;
			else if (complexity == 2) return 30;
		}
		else if (element == "Pa") {
			if (complexity == 1) return 13;
			else if (complexity == 2) return 31;
		}
		else if (element == "U") {
			if (complexity == 1) return 14;
			else if (complexity == 2) return 32;
		}
		else if (element == "Np") {
			if (complexity == 1) return 15;
			else if (complexity == 2) return 33;
		}
		else if (element == "Pu") {
			if (complexity == 1) return 16;
			else if (complexity == 2) return 34;
		}
		else if (element == "Am") {
			if (complexity == 1) return 17;
			else if (complexity == 2) return 35;
		}
		else if (element == "Cm") {
			if (complexity == 1) return 18;
			else if (complexity == 2) return 36;
		}
		else if (element == "Bk") {
			if (complexity == 1) return 20;
			else if (complexity == 2) return 38;
		}
		else if (element == "Cf") {
			//no basis?
		}
		else if (element == "Es") {
			if (complexity == 1) return 21;
			else if (complexity == 2) return 39;
		}
		else if (element == "Fm") {
			if (complexity == 1) return 22;
			else if (complexity == 2) return 40;
		}
		else if (element == "Md") {
			if (complexity == 1) return 23;
			else if (complexity == 2) return 41;
		}
		else if (element == "No") {
			if (complexity == 1) return 24;
			else if (complexity == 2) return 42;
		}
		else if (element == "Lr") {
			if (complexity == 1) return 25;
			else if (complexity == 2) return 43;
		}
		else if (element == "Rf") {
			//no basis
		}
		else if (element == "Db") {
			//no basis
		}
		else if (element == "Sg") {
		//no basis
		}
		else if (element == "Bh") {
		//no basis
		}
		else if (element == "Hs") {
		//no basis
		}
		else if (element == "Mt") {
		//no basis
		}
		else if (element == "Ds") {
		//no basis
		}
		else if (element == "Rg") {
		//no basis
		}
		else if (element == "112") {
		//no basis
		}
		else if (element == "113") {
		//no basis
		}
		else if (element == "114") {
		//no basis
		}
		else if (element == "115") {
		//no basis
		}
		else if (element == "116") {
		//no basis
		}
		else if (element == "117") {
		//no basis
		}
		else if (element == "118") {
		//no basis
		}
		return 0;
}

std::string getBasisSet(std::string element, int complexity)
{
	return "GTH-PADE-q" + std::to_string(getBasisSetElectrons(element, complexity));
}

void writeToCP2K(structure& s, std::string task, std::string filename, std::string basis_file, std::string potential_file, double a = 0, double b = 0, double c = 0)
{
	std::string output;
	std::ofstream file(filename);

	std::string global = "&GLOBAL\n\tPROJECT " + task + "\n\tRUN_TYPE ENERGY_FORCE\n\tPRINT_LEVEL LOW\n&END GLOBAL";

	std::string DFT = "\n\t&DFT"
		"\n\t\tBASIS_SET_FILE_NAME " + basis_file;
		"\n\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS " + potential_file;
		"\n\t\t&QS"
		"\n\t\t\tEPS_DEFAULT 1.0E-10"
		"\n\t\t&END QS"
		"\n\t\t&MGRID"
		"\n\t\t\tNGRIDS 4"
		"\n\t\t\tCUTOFF 300"
		"\n\t\t\tREL_CUTOFF 60"
		"\n\t\t&END MGRID"
		"\n\t\t&XC"
		"\n\t\t\t&XC_FUNCTIONAL PADE"
		"\n\t\t\t&END XC_FUNCTIONAL"
		"\n\t\t&END XC"
		"\n\t\t&SCF"
		"\n\t\t\tSCF_GUESS ATOMIC"
		"\n\t\t\tEPS_SCF 1.0E-7"
		"\n\t\t\tMAX_SCF 300"
		"\n\t\t\t&DIAGONALIZATION ON"
		"\n\t\t\t\tALGORITHM STANDARD"
		"\n\t\t\t&END DIAGONALIZATION"
		"\n\t\t\t&MIXING T"
		"\n\t\t\t\tMETHOD BROYDEN_MIXING"
		"\n\t\t\t\tALPHA 0.4"
		"\n\t\t\t\tNBROYDEN 8"
		"\n\t\t\t&END MIXING"
		"\n\t\t&END SCF"
		"\n\t&END DFT";

	std::string subsys = "\n\t&SUBSYS";
	for (int e = 0; e < s.elements.size(); e++)
	{
		subsys += "\n\t\t&KIND " + s.elements[e];
		subsys += "\n\t\t\tELEMENT\t" + s.elements[e];
		subsys += "\n\t\t\tBASIS_SET DZVP-GTH-PADE"
			"\n\t\t\tPOTENTIAL " + getBasisSet(s.elements[e], 1) +
			"\n\t\t&END KIND";
	}

	/*
	We need to know the cell shape. From a materials science perspective it makes sense for this to be the smallest shape that will fit all the atoms in side. Determine this.
	*/
	if (a == 0 && b == 0 && c == 0)
	{
		double maxa = 0;
		double maxb = 0;
		double maxc = 0;
		double mina = DBL_MAX;
		double minb = DBL_MAX;
		double minc = DBL_MAX;

		for (int e = 0; e < s.set.size(); e++)
		{
			for (int atom = 0; atom < s.set[e].size(); atom++)
			{
				if (s.set[e][atom].x > maxa) {
					maxa = s.set[e][atom].x;
				}
				else if (s.set[e][atom].x < mina)
				{
					mina = s.set[e][atom].x;
				}

				if (s.set[e][atom].y > maxb) {
					maxb = s.set[e][atom].y;
				}
				else if (s.set[e][atom].y < minb)
				{
					minb = s.set[e][atom].y;
				}

				if (s.set[e][atom].z > maxc) {
					maxc = s.set[e][atom].z;
				}
				else if (s.set[e][atom].z < minc)
				{
					minc = s.set[e][atom].z;
				}

			}
		}
		//double the cell shape to account for the negative numbers?
		a = maxa - mina;
		b = maxb - minb;
		c = maxc - minc;
	}

	subsys += "\n\t\t&CELL"
		"\n\t\t\tABC " + std::to_string(a) + " " + std::to_string(b) + " " + std::to_string(c) +
		"\n\t\t&END CELL"
		//"\n\t\t# " + std::to_string(numMolecules) + " " + formula + "(TIP5P," + std::to_string(pressureBars) + "bar, " + std::to_string(temperature) + "K) a = 9.8528"
		"\n\t\t&COORD";
	for (int e = 0; e < s.set.size(); e++)
	{
		for (int atom = 0; atom < s.set[e].size(); atom++)
		{
			subsys += "\n\t" + s.elements[e] + " " + std::to_string(s.set[e][atom].x) + " " + std::to_string(s.set[e][atom].y) + " " + std::to_string(s.set[e][atom].z);
		}
	}
	subsys += "\n\t\t&END COORD"
		"\n\t&END SUBSYS";

	std::string print = "\n\t&PRINT"
		"\n\t\t&FORCES ON"
		"\n\t\t&END FORCES"
		"\n\t&END PRINT";

	std::string forceeval = "&FORCE_EVAL\n\tMETHOD Quickstep" + DFT + subsys + print + "\n&END FORCE_EVAL";
	file << global << "\n" << forceeval;
	file.close();
}

void writeToGuassian(structure& s, std::string task, std::string filename, std::string basis_name,std::string method, int charge = 0, int state = 1, int maxCycles = 0, int maxSCF = 129)
{
	//max cycles = 0 implies this is just energy and no optimization
	//a,b,c are useless for guassian? so is potential name. the variables are left here for compatibility.
	/*
	
	# HF/6-31G(d)	Route section
 
water energy	Title section
 
0   1	Molecule specification
O  -0.464   0.177   0.0	 
H  -0.464   1.137   0.0	 
H   0.441  -0.143   0.0
	*/
	//route section
	std::string outString = "# " + method + "/" + basis_name;
	if(maxCycles > 0) outString += " Opt=(MaxCycles=" + std::to_string(maxCycles) + ")";
	if (maxSCF > 0) outString += " SCF = (MaxCycle = " + std::to_string(maxSCF) + ")";
	outString += "\n\n";
	//title section
	outString += task + "\n\n";
	//molecule specification, charge, state singlet doublet triplet
	outString += std::to_string(charge) + " " + std::to_string(state);
	//molecules
	for (int e = 0; e < s.set.size(); e++)
	{
		for (int atom = 0; atom < s.set[e].size(); atom++)
		{
			outString += "\n" + s.elements[e] + " " + std::to_string(s.set[e][atom].x) + " " + std::to_string(s.set[e][atom].y) + " " + std::to_string(s.set[e][atom].z);
		}
	}
	outString += "\n\n";//blank line at the end of the structure is required
	std::ofstream file(filename);
	file << outString;
	file.close();
	


}
void writeToADF(structure& s, std::string task, std::string filename, std::string basis_name, std::string method,int maxSCF = 300, int charge = 0, int state = 1, std::string converge = "1.0e-6")
{

	//task is either GeometryOptimization or Energy??
	std::string out = "#!/bin/sh\n\n\"$AMSBIN/ams\" << eor\n\n";
	out += "Task " + task + "\n";
	out += "System\n";
	out += "\tCharge " + std::to_string(charge) + "\n";
	out += "\tAtoms";


	for (int e = 0; e < s.set.size(); e++)
	{
		for (int atom = 0; atom < s.set[e].size(); atom++)
		{
			out += "\n\t\t" + s.elements[e];
			if (s.set[e][atom].x > 0) out += "  " + std::to_string(s.set[e][atom].x).substr(0, 7);
			else out += " " + std::to_string(s.set[e][atom].x).substr(0, 8);
			if (s.set[e][atom].y > 0) out += "  " + std::to_string(s.set[e][atom].y).substr(0, 7);
			else out += " " + std::to_string(s.set[e][atom].y).substr(0, 8);
			if (s.set[e][atom].z > 0) out += "  " + std::to_string(s.set[e][atom].z).substr(0, 7);
			else out += " " + std::to_string(s.set[e][atom].z).substr(0, 8);
		}
	}
	out += "\n\tEnd\nEnd\n";
	/*
	if (task == "GeometryOptimization")
	{
		out += "GeometryOptimization\n";
		out += "\toptim all\n";
		out += "\titerations 300\n";
		out += "\tstep rad=0.15 angle=10.0\n";
		out += "\thessupd BFGS\n";
		out += "\tConvergence Gradient=1.0e-4\n";
		out += "\tConvergence e=1.0e-4 grad=1.0e-4 rad=1.0e-2 angle=0.5\n";
		out += "End\n";
		out += "\t\n";
		out += "SCF\n";
		out += "\titerations 300\n";
		out += "\tconverge 1.0e-6 1.0e-6\n";
		out += "\tmixing 0.2\n";
		out += "\tlshift 0.0\n";
		out += "\tDIIS\n\t\tn 500\n\t\tok 0.00001\n\t\tcyc 500\n\tEnd\n";// cx = 5.0 cxx = 10.0\n";
		out += "END\n";
	}*/

	//Engine block
	out += "Engine ADF\n";//should be ADF
	out += "\tUnrestricted Yes\n";
	out += "\tspinpolarization " + std::to_string(state) + "\n";
	out += "\tBasis\n\t\tType " + basis_name + "\n\t\tCore Large\n";
	out += "\tEnd\n";
	out += "\trelativity\n";
	out += "\t\tlevel scalar\n";
	out += "\t\tformalism ZORA\n";
	out += "\tEnd\n";
	out += "\tSCF\n";
	out += "\t\titerations " + std::to_string(maxSCF) + "\n";
	out += "\t\tconverge " + converge + " " + converge + "\n";
	out += "\t\tmixing 0.2\n";
	out += "\t\tlshift 0.0\n";
	out += "\t\tdiis\n\t\t\tn 500\n\t\t\tok 0.00001\n\t\t\tcyc 500\n\t\tEnd\n";// cx = 5.0 cxx = 10.0\n";
	out += "\tEND\n";

	out += "\txc\n";
	out += "\t\tgga " + method + "\n";
	out += "\tEnd\n";
	out += "EndEngine\n";
	out += "eor";
	std::ofstream file(filename + ".run");
	file << out;
	file.close();



}

structure* optimizationADF(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, std::string AMSHOME,int amscalculation = -1,int maxSCF = 300,std::string converge="1.0e-6", double timeoutminutes = -1, double timeLeftMinutes = -1, bool keepFiles = false)
{
	
	if (amscalculation != -1)
	{
		//we are going to make separate folders for each ams

	}
	else
	{
		char* ac;
		std::string deleteFiles1 = "rm -r ams.results\n";
		std::string deleteFiles2 = "rm " + filename + ".run\n";
		std::string deleteFiles3 = "rm " + filename + ".out\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


		ac = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			ac[i] = deleteFiles2[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

		ac = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			ac[i] = deleteFiles3[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

	}
	/*
	
	*/
	std::string dirname = "calc" + std::to_string(amscalculation);
	char* ac;
	if (amscalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}

	writeToADF(s, "GeometryOptimization", dirname + "/" + filename, basis_name, method,maxSCF, charge, state,converge);
	
	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + filename + ".run");

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + filename + ".run");
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	//commands += "export AMSHOME=" + AMSHOME + "\n";//AMSHOME = /home/jpburkhardt/scm/ams2023.104
	//commands += "export AMSBIN=$AMSHOME/bin\nexport AMSRESOURCES=$AMSHOME/atomicdata\nexport SCMLICENSE=$AMSHOME/license.txt\n";
	//commands += "AMS_CALCULATION=" + filename + "\n" + "dos2unix $AMS_CALCULATION.run\n" + "chmod u+x $AMS_CALCULATION.run\n";
	//commands += "./$AMS_CALCULATION.run > " + filename + ".out\n";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + filename + ".run\n" + "chmod u+x " + filename + ".run\n";
	commands += "nohup ./" + filename + ".run > " + filename + ".out &\n";

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;


	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;


	//reading file input
	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + filename + ".out");


	//stall for cp2k to finish
	std::cout << "checking for guassian completion" << std::endl;



	while (!getline(*aOutFile, checker))
	{
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");
	}
	std::cout << "ADF started" << std::endl;

	//the file is writing complete?
	bool finishedOutput_energy = false;
	bool finishedOutput_geometry = false;
	std::vector<std::string> optStrings = {};
	std::string outLine;
	std::string energyString;




	while (!finishedOutput_energy || !finishedOutput_geometry)
	{
		if (timeLeftMinutes != -1)
		{
			
			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes*60 - 5*60)
			{
				//kill the task so we dont waste cpu resource
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return nullptr;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return nullptr;
			}

		}

		optStrings = {};
		//reopen file to get new updates
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");

		//look for the energy line
		
		while (getline(*aOutFile, outLine)) {
			
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = "Energy (hartree)";
			std::string match2 = "Optimized geometry:";
			std::string match3 = "Geometry optimization failed";
			std::string match4 = "BAD CORE INTEGRAL";
			std::string match5 = "Process received SIGTERM";
			if (outLine.size() >= match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) {
					finishedOutput_energy = true;//cannot be made false
					energyString = outLine;
				}
			}

			bool getopt = false;

			if (outLine.size() >= match2.size())
			{
				bool matching = true;
				for (int m = 0; m < match2.size(); m++)
				{
					if (outLine[m] != match2[m]) matching = false;
				}
				if (matching) {
					getopt = true;//cannot be made false
					finishedOutput_geometry = true;//geometry needs its own bool now because there are cases where energy calculates but geometry does not and crash line geometry optimization failed hasnt been printed yet
				}
			}
			if (getopt)
			{
				bool stop = false;
				for (int i = 0; i < 7; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
				for (int e = 0; e < s.elements.size(); e++)  for (int a = 0; a < s.set[e].size(); a++) {
					if(!stop) if (!getline(*aOutFile, outLine)) stop = true;
					optStrings.push_back(outLine);
				}
			}
			//when the AMS calculation is successful, energy is always printed after the geometry, and therefore its fine to get the geometry here.

			bool failed = false;
			if (outLine.size() >= match3.size())
			{
				bool matching = true;
				for (int m = 0; m < match3.size(); m++)
				{
					if (outLine[m] != match3[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (outLine.size() >= match4.size())
			{
				bool matching = true;
				for (int m = 0; m < match4.size(); m++)
				{
					if (outLine[m] != match4[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (outLine.size() >= match5.size())
			{
				bool matching = true;
				for (int m = 0; m < match5.size(); m++)
				{
					if (outLine[m] != match5[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (failed)
			{
				return nullptr;
			}
		}

		

	}
	std::cout << "finished output: " << finishedOutput_energy << std::endl;



	double seconds_since_start = difftime(time(0), start);

	std::string getout = "cd ../\n";
	if (amscalculation != -1) {
		ac = new char[getout.length() + 1];
		for (int i = 0; i < getout.length(); i++)
		{
			ac[i] = getout[i];
		}

		std::system(ac);
		std::cout << "command: " << getout << std::endl;
		delete[] ac;
	}
	std::cout << "ADF time: " << seconds_since_start << std::endl;
	
	structure* st = nullptr;

	

	std::vector<std::string> elements = {};
	std::vector<std::vector<atom>> set = {};

	for (int i = 0; i < optStrings.size(); i++)
	{
		std::cout << optStrings[i] << std::endl;

		std::string cline = optStrings[i];

		std::string element = "";
		if (cline[12] == ' ') element = cline.substr(13, 1);//it goes the other way ....
		else element = cline.substr(12, 2);
		//x
		double x;
		double y;
		double z;
		bool sign = false;
		if (cline[18] == '-') sign = true;
		if (sign) x = std::stod(cline.substr(18, 11));
		else x = std::stod(cline.substr(19, 10));
		if (cline[33] == '-') sign = true;
		else sign = false;
		if (sign) y = std::stod(cline.substr(33, 11));
		else y = std::stod(cline.substr(34, 10));
		if (cline[48] == '-') sign = true;
		else sign = false;
		if (sign) z = std::stod(cline.substr(48, 11));
		else z = std::stod(cline.substr(49, 10));

		

		//figure out which set this belongs to.
		int index = -1;
		for (int i = 0; i < elements.size(); i++)
		{
			if (elements[i] == element)
			{
				index = i;
			}
		}
		//if (buggo) std::cout << "making atom for elemental index " << index << " with coords " << x << " " << y << " " << z << std::endl;
		if (index == -1)
		{
			//if (buggo) std::cout << "not yet indexed" << std::endl;
			//atom type not yet indexed
			elements.push_back(element);
			std::vector < atom> alist = { atom(x,y,z) };
			set.push_back(alist);
		}
		else {
			//if (buggo) std::cout << "pushing back" << std::endl;
			set[index].push_back(atom(x, y, z));
		}
	}
				//if (buggo) std::cout << "this structure is not the last one starting over with next" << std::endl;

	if (st != nullptr) delete st;
	st = new structure(set, elements);

	
	std::cout << "found energy" << std::endl;
	std::string matchString = "Energy (hartree)";
	std::string num = "";
	for (int i = matchString.length(); i < energyString.length(); i++)
	{
		if (energyString[i] != ' ') num += energyString[i];
	}
	double energy = std::stod(num);
	st->energy = energy;
	// Close the file
	//if there was no convergence then nullptr will be returned
	//if (buggo) std::cout << "closing" << std::endl;
	aOutFile->close();
	delete aOutFile;


	//remove files
	if (!keepFiles)
	{
		if (amscalculation == -1)
		{

			std::string deleteFiles1 = "rm -r ams.results\n";
			std::string deleteFiles2 = "rm " + filename + ".run\n";
			std::string deleteFiles3 = "rm " + filename + ".out\n";

			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;


			ac = new char[deleteFiles2.length() + 1];
			for (int i = 0; i < deleteFiles2.length(); i++)
			{
				ac[i] = deleteFiles2[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;

			ac = new char[deleteFiles3.length() + 1];
			for (int i = 0; i < deleteFiles3.length(); i++)
			{
				ac[i] = deleteFiles3[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;
		}
		else {
			std::string deleteFiles1 = "rm -r " + dirname + "\n";
			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;

		}

	}
	
	return st;
	//return energyHartree;
}

double energyADF(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, std::string AMSHOME,int adfcalculation = -1,std::string converge = "1.0e-6",int timeoutminutes = -1,int timeLeftMinutes = -1, bool keepFiles = false)
{
	if (adfcalculation == -1)//always delete with ams, make sure its gone to prevent crashes
	{
		char* ac;
		std::string deleteFiles1 = "rm -r ams.results\n";
		std::string deleteFiles2 = "rm " + filename + ".run\n";
		std::string deleteFiles3 = "rm " + filename + ".out\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


		ac = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			ac[i] = deleteFiles2[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

		ac = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			ac[i] = deleteFiles3[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

	}
	std::string dirname = "calc" + std::to_string(adfcalculation);
	char* ac;
	if (adfcalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\n";//cd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;

	}

	writeToADF(s, "SinglePoint", dirname + "/" + filename, basis_name, method,300, charge, state,converge);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + filename + ".run");

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + filename + ".run");
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	//commands += "export AMSHOME=" + AMSHOME + "\n";//AMSHOME = /home/jpburkhardt/scm/ams2023.104
	//commands += "export AMSBIN=$AMSHOME/bin\nexport AMSRESOURCES=$AMSHOME/atomicdata\nexport SCMLICENSE=$AMSHOME/license.txt\n";
	//commands += "AMS_CALCULATION=" + filename + "\n";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + filename + ".run\n" + "chmod u+x " + filename + ".run\n";
	commands += "nohup ./" + filename + ".run > " + filename + ".out &\n";

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;


	//reading file input
	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + filename + ".out");


	//stall for cp2k to finish
	std::cout << "checking for adf completion" << std::endl;




	while (!getline(*aOutFile, checker))
	{
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");
	}
	std::cout << "ADF started" << std::endl;

	//the file is writing complete?
	bool finishedOutput = false;
	std::string outLine;
	std::string energyString;
	while (!finishedOutput)
	{
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				//kill the task so we dont waste cpu resource
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;
			}

		}


		//reopen file to get new updates
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");

		//look for the energy line

		while (getline(*aOutFile, outLine)) {
			//come back here plz
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = "Energy (hartree)";
			std::string match4 = "BAD CORE INTEGRAL";
			std::string match5 = "Process received SIGTERM";
			if (outLine.size() >= match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) {
					finishedOutput = true;//cannot be made false
					energyString = outLine;
				}
			}
			bool failed = false;
			if (outLine.size() >= match4.size())
			{
				bool matching = true;
				for (int m = 0; m < match4.size(); m++)
				{
					if (outLine[m] != match4[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (outLine.size() >= match5.size())
			{
				bool matching = true;
				for (int m = 0; m < match5.size(); m++)
				{
					if (outLine[m] != match5[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (failed)
			{
				return DBL_MAX;
			}
		}



	}
	std::cout << "finished output: " << finishedOutput << std::endl;


	double seconds_since_start = difftime(time(0), start);
	if (adfcalculation != -1) {
		std::string getout = "cd ../\n";

		ac = new char[getout.length() + 1];
		for (int i = 0; i < getout.length(); i++)
		{
			ac[i] = getout[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}
	std::cout << "ADF time: " << seconds_since_start << std::endl;


	std::cout << "found energy" << std::endl;
	std::string matchString = "Energy (hartree)";
	std::string num = "";
	for (int i = matchString.length(); i < energyString.length(); i++)
	{
		if (energyString[i] != ' ') num += energyString[i];
	}
	double energy = std::stod(num);
	// Close the file
	//if there was no convergence then nullptr will be returned
	//if (buggo) std::cout << "closing" << std::endl;
	aOutFile->close();
	delete aOutFile;


	//remove files
	if (!keepFiles)
	{
		if (adfcalculation == -1)
		{

			std::string deleteFiles1 = "rm -r ams.results\n";
			std::string deleteFiles2 = "rm " + filename + ".run\n";
			std::string deleteFiles3 = "rm " + filename + ".out\n";

			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;


			ac = new char[deleteFiles2.length() + 1];
			for (int i = 0; i < deleteFiles2.length(); i++)
			{
				ac[i] = deleteFiles2[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;

			ac = new char[deleteFiles3.length() + 1];
			for (int i = 0; i < deleteFiles3.length(); i++)
			{
				ac[i] = deleteFiles3[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;
		}
		else {
			std::string deleteFiles1 = "rm -r " + dirname + "\n";
			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;

		}

	}

	return energy;
	//return energyHartree;
}

double basinHoppingEnergy(structure& s, int ittt, std::string task, double a, double b, double c)
{
	if (DEBUG > 1)
	{
		//scannign structure
		for (std::vector<std::vector<atom>>::iterator it = s.set.begin(); it != s.set.end(); it++)
		{
			std::cout << "element set" << std::endl;
			for (int i = 0; i < it->size(); i++)
			{
				std::cout << "atom: " << (*it)[i].x << " " << (*it)[i].y << " " << (*it)[i].z << std::endl;
			}
		}
		std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
		std::string inFileName = fileName + "_in.txt";
		writeToCP2K(s, task, inFileName, "/home/jpburkhardt/GMS/data/BASIS_SET", "/home/jpburkhardt/GMS/data/GTH_POTENTIALS", a, b, c);

	}

	/*
	New basin hopping energy iterator
	*/
	double energy = 0;
	double sigma = 1;
	double epsilon = 1;
	if (DEBUG > 1) std::cout << "starting energy calculation" << std::endl;

	//the iterators i made for structure only work with do while loops
	structure::Iterator it = s.begin();
	do {
		if (DEBUG > 3) std::cout << "new jt loop about to begin" << std::endl;
		structure::Iterator jt = s.begin();
		do {
			if (DEBUG > 3) std::cout << "comparing" << "pointer " << it.m_ptr << " to " << "pointer " << jt.m_ptr << std::endl;
			if (it != jt)
			{
				//if (DEBUG) std::cout << "doing some math" << std::endl;

				if (DEBUG > 2) std::cout << "atom " << it.m_ptr << ": " << it->x << " " << it->y << " " << it->z << " vs atom2 " << jt.m_ptr << ": " << jt->x << " " << jt->y << " " << jt->z << std::endl;
				double r = euclideanDistance((*it), (*jt));
				//if (DEBUG) std::cout << "is r 0: " << r << std::endl;
				//std::cout << "r: " << r << std::endl;
				double e = pow(sigma / r, 12) - pow(sigma / r, 6);
				if (DEBUG > 1) std::cout << it.indexI << "," << it.indexJ << " vs. " << jt.indexI << "," << jt.indexJ << " e: " << e << std::endl;
				energy += e;
			}
			else {
				if (DEBUG > 1) std::cout << "did not calculate an atom against itself" << std::endl;
			}
			if (DEBUG > 3) std::cout << "about to iterate jt" << std::endl;
			jt++;
		} while (jt != s.end());

		it++;
	} while (it != s.end());



	/*
	for (std::vector<std::vector<atom>>::iterator it = s.set.begin(); it != s.set.end(); it++)
	{
		if (DEBUG) std::cout << "iterating outer set" << std::endl;
		for (int i = 0; i < it->size(); i++)
		{
			if (DEBUG) std::cout << "inner loop 1: all inner sets" << std::endl;
			for (std::vector<std::vector<atom>>::iterator itt = s.set.begin(); itt != s.set.end(); itt++)
			{
				if (DEBUG) std::cout << "inner loop 2: inner sets size: " << itt->size() << std::endl;
				for (int j = 0; j < itt->size(); j++)
				{
					if (DEBUG) std::cout << "about to do math" << std::endl;
					if (&(*it)[i] != &(*itt)[j])
					{
						if (DEBUG) std::cout << "doing some math" << std::endl;
						if (DEBUG) std::cout << "it size:" << it->size() << " " << i << std::endl;
						if (DEBUG) std::cout << "itt size:" << itt->size() << " " << j << std::endl;
						double r = euclideanDistance((*it)[i], (*itt)[j]);
						//if (DEBUG) std::cout << "is r 0: " << r << std::endl;
						//std::cout << "r: " << r << std::endl;
						double e = pow(sigma / r, 12) - pow(sigma / r, 6);
						//std::cout << "e: " << e << std::endl;
						energy += e;
					}
					else {
						if (DEBUG) std::cout << "did not calculate duplicate" << std::endl;
					}
					if (DEBUG) std::cout << "math done" << std::endl;
				}
				if (DEBUG) std::cout << "here 2" << std::endl;
			}
			if (DEBUG) std::cout << "here 1" << std::endl;
		}
		if (DEBUG) std::cout << "outer loop complete" << std::endl;
	}startModifiedBasinHopping(init, x, y, z, percentChangeI, percentChangeF, composition, elements, charge, state, dft, method, basis, optimizationCycles, coordinationSteps, HFsteps, time, localMinimaTrappingSteps, radialCriteriaPercent, outputFileName, taskName, a, b, c, blindGradientDescent, directionExponent, proportionExponent, deltaEThr);
			
	*/
	if (DEBUG > 1) std::cout << "returning" << std::endl;
	return 0.5 * epsilon * energy;
};

structure* optimizationBasinHopping(structure& s, int ittt, std::string task, int maxCycles, double x, double y, double z)
{
	//perform max cylces random steps
	std::vector<double> movementVector = {};
	double lowestEnergy = s.energy;
	structure * st = nullptr;
	for (int i = 0; i < s.elements.size(); i++)
	{
		for (int j = 0; j < s.set[i].size(); j++)
		{
			movementVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * x));
			movementVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * y));
			movementVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * z));
		}

	}
	for (int i = 0; i < maxCycles; i++)
	{
		for (int i = 0; i < s.elements.size(); i++)
		{
			for (int j = 0; j < s.set[i].size(); j++)
			{
				movementVector[i*3] = (((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)* x));
				movementVector[i*3 + 1] = (((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * y));
				movementVector[i*3 + 2] = (((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * z));
			}

		}
		structure* sn = nullptr;
		if(st != nullptr) sn = new structure(*st, movementVector);
		else sn = new structure(s, movementVector);
		double snergy = basinHoppingEnergy(*sn, 0, "task", x, y, z);
		sn->energy = snergy;
		if (snergy < lowestEnergy)
		{
			if(st != nullptr) delete st;
			st = sn;
			lowestEnergy = snergy;
		}
	}
	return st;
}



double energyCP2K(structure& s, int ittt, std::string task, double a, double b, double c)
{

	//these a b c are for cell shape
	std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
	std::string inFileName = fileName + "_in.txt";
	std::string outFileName = fileName + "_out.txt";
	writeToCP2K(s, task, inFileName, "/home/jpburkhardt/GMS/data/BASIS_SET", "/home/jpburkhardt/GMS/data/GTH_POTENTIALS", a, b, c );

	//ensure file has been written
	std::string checker;
	std::ifstream* cp2kInput = new std::ifstream( inFileName);

	while (!getline(*cp2kInput, checker))
	{
		cp2kInput->close();
		delete cp2kInput;
		cp2kInput = new std::ifstream( inFileName);
	}
	cp2kInput->close();
	delete cp2kInput;
	std::cout << "input file written" << std::endl;

	std::string cp2kcommand = "cp2k.psmp -i " + inFileName + " - o " + outFileName;
		
//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
	//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	char* cp2kc = new char[cp2kcommand.length() + 1];
	for (int i = 0; i < cp2kcommand.length(); i++)
	{
		cp2kc[i] = cp2kcommand[i];
	}
	cp2kc[cp2kcommand.length()] = '\0';
	std::cout << "command: " << cp2kc << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "command: " << cp2kc << std::endl;
	std::system(cp2kc);
	std::cout << "command written" << std::endl;

	std::ifstream* cp2koutFile = new std::ifstream( outFileName + ".txt");


	//stall for cp2k to finish
	std::cout << "checking for cp2k completion" << std::endl;
	time_t start = time(0);

	
	while (!getline(*cp2koutFile, checker))
	{
		cp2koutFile->close();
		delete cp2koutFile;
		cp2koutFile = new std::ifstream( outFileName);
	}
	std::cout << "cp2k started" << std::endl;
	//the file is located but energy may not be calculated yet
	bool energyCalculated = false;
	std::string energyLine;
	while (!energyCalculated)
	{
		//reopen file to get new updates
		cp2koutFile->close();
		delete cp2koutFile;
		cp2koutFile = new std::ifstream(outFileName);

		//look for the energy line
		while (getline(*cp2koutFile, energyLine)) {
			std::string match = "ENERGY| Total";
			if (energyLine.size() > 6)
			{
				if (energyLine[1] == 'E' && energyLine[2] == 'N' && energyLine[3] == 'E' && energyLine[4] == 'R' && energyLine[5] == 'G' && energyLine[6] == 'Y')
				{
					energyCalculated = true;
				}
			}
		}
	}

	double seconds_since_start = difftime(time(0), start);
	std::cout << "cp2k time: " << seconds_since_start << std::endl;
	std::cout << "reading file" << std::endl;

	double value = 0.0;
	// Use a while loop together with the getline() function to read the file line by line

	double enrg = 0.0;
	int count = 0;
	int j = energyLine.size() - 1;
	while (energyLine[j - 1] != ' ')
	{
		j--;
	}
	
	std::string val;
	while (j != energyLine.size())
	{
		val += energyLine[j];
		j++;
	}
	std::cout << val << std::endl;
	value = std::stod(val);


	// Close the file
	cp2koutFile->close();
	delete cp2koutFile;
	std::cout << "energy: " << value << std::endl;
	return value;
}


double energyGuassian(structure& s, int ittt, std::string task,std::string basis,std::string method,int charge = 0, int state = 0, int maxSCF = 129,bool keepFiles = false)
{
	//outfile name should end in log, not sure what the input should end in

	//these a b c are for cell shape
	std::string guassianTaskName = task + " energy";

	std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
	std::string inFileName = fileName + ".com";
	std::string outFileName = fileName + ".log";

	writeToGuassian(s, guassianTaskName, inFileName, basis, method,charge,state,0,maxSCF);

	//ensure file has been written
	std::string checker;
	std::ifstream* guassianInput = new std::ifstream(inFileName);

	while (!getline(*guassianInput, checker))
	{
		guassianInput->close();
		delete guassianInput;
		guassianInput = new std::ifstream(inFileName);
	}
	guassianInput->close();
	delete guassianInput;
	std::cout << "input file written" << std::endl;

	std::string guassiancommand = "g16 <" + inFileName + " >" + outFileName;

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
	//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	char* gc = new char[guassiancommand.length() + 1];
	for (int i = 0; i < guassiancommand.length(); i++)
	{
		gc[i] = guassiancommand[i];
	}
	gc[guassiancommand.length()] = '\0';
	std::cout << "command: " << gc << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "command: " << gc << std::endl;
	time_t start = time(0);
	std::system(gc);
	delete[]gc;
	//std::cout << "command written" << std::endl;

	std::ifstream* gOutFile = new std::ifstream(outFileName);


	//stall for cp2k to finish
	std::cout << "checking for guassian completion" << std::endl;
	


	while (!getline(*gOutFile, checker))
	{
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);
	}
	std::cout << "guassian started" << std::endl;
	//the file is located but energy may not be calculated yet
	bool finishedOutput = false;

	std::string outLine;
	while (!finishedOutput)
	{
		//reopen file to get new updates
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);

		//look for the energy line
		while (getline(*gOutFile, outLine)) {
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = " Job cpu time:";

			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
			}
		}

	}
	double seconds_since_start = difftime(time(0), start);
	std::cout << "gaussian time: " << seconds_since_start << std::endl;

	std::string energyLine;
	int lastLines = 0;
	int clines = 0;
	//reopen file to get new updates
	gOutFile->close();
	delete gOutFile;
	gOutFile = new std::ifstream(outFileName);
	bool energyCalculated = false;
		//look for the energy line
	while (!energyCalculated) {
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line
		if (getline(*gOutFile, energyLine))
		{
			std::string match = " SCF Done:  E";

			if (energyLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (energyLine[m] != match[m]) matching = false;
				}
				energyCalculated = matching;
			}
		}
		else {
			break;
		}
			
	}
	if (!energyCalculated) return DBL_MAX;//this might mess with gradient descent

	
	std::cout << "reading file" << std::endl;

	//find the start of the energy
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if(eFin){}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;
		
	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));

	


	// Close the file
	gOutFile->close();
	delete gOutFile;
	std::cout << "energy: " << energyHartree << std::endl;

	//remove files
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm " + inFileName + "\n";
		std::string deleteFiles2 = "rm " + outFileName + "\n";
		std::string deleteFiles3 = "rm -r core.*\n";

		gc = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			gc[i] = deleteFiles1[i];
		}
		std::cout << "command: " << gc << std::endl;
		std::system(gc);
		delete[] gc;


		gc = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			gc[i] = deleteFiles2[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

		gc = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			gc[i] = deleteFiles3[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

	}

	return energyHartree;
}

structure* optimizationGaussian(structure& s, int ittt, std::string task, std::string basis, std::string method, int maxCycles, int charge = 0, int state = 0, int maxSCF = 129, bool keepFiles = false)
{
	//the method should be HF instead of DFT for better convergence
	//outfile name should end in log, not sure what the input should end in


	std::string guassianTaskName = task + " opt";

	std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
	std::string inFileName = fileName + ".com";
	std::string outFileName = fileName + ".log";

	writeToGuassian(s, guassianTaskName, inFileName, basis, method, charge, state, maxCycles,maxSCF);

	//ensure file has been written
	std::string checker;
	std::ifstream* guassianInput = new std::ifstream(inFileName);

	while (!getline(*guassianInput, checker))
	{
		guassianInput->close();
		delete guassianInput;
		guassianInput = new std::ifstream(inFileName);
	}
	guassianInput->close();
	delete guassianInput;
	std::cout << "input file written" << std::endl;

	std::string guassiancommand = "g16 <" + inFileName + " >" + outFileName;

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	char* gc = new char[guassiancommand.length() + 1];
	for (int i = 0; i < guassiancommand.length(); i++)
	{
		gc[i] = guassiancommand[i];
	}
	gc[guassiancommand.length()] = '\0';
	std::cout << "command: " << gc << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "command: " << gc << std::endl;
	time_t start = time(0);
	std::system(gc);
	delete[] gc;
	//std::cout << "command written" << std::endl;

	std::ifstream* gOutFile = new std::ifstream(outFileName);


	//stall for cp2k to finish
	std::cout << "checking for guassian completion" << std::endl;



	while (!getline(*gOutFile, checker))
	{
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);
	}
	std::cout << "guassian started" << std::endl;

	//the file is writing complete?
	bool finishedOutput = false;

	std::string outLine;
	while (!finishedOutput)
	{
		//reopen file to get new updates
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);

		//look for the energy line
		while (getline(*gOutFile, outLine)) {
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = " Job cpu time:";

			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if(matching) finishedOutput = true;//cannot be made false
			}
		}

	}
	

	double seconds_since_start = difftime(time(0), start);
	std::cout << "gaussian time: " << seconds_since_start << std::endl;
	std::cout << "reading file" << std::endl;

	//find the start of the optimization
	//this will take the orientation everytime but iguess thats ok
	//restart the file just incase
	gOutFile->close();
	delete gOutFile;
	gOutFile = new std::ifstream(outFileName);
	std::string cline;
	structure *st = nullptr;

	//bool buggo = true;
	//if (buggo) std::cout << "opening line one" << std::endl;
	std::string matchString = "Input orientation:";
	std::string eMatch = " SCF Done:  E";
	std::string energyLine = "";



	std::cout << "reading file" << std::endl;

	//find the start of the energy
	
	while (getline(*gOutFile, cline)) {
		//if (buggo) std::cout << "line one opened" << std::endl;
		if (cline.length() > matchString.length())
		{

			bool matching = false;
			int matchIt = 0;
			//if (buggo) std::cout << "checking matches" << std::endl;

			for (int l = 0; l < cline.length() - matchString.length(); l++)
			{
				//if (buggo) std::cout << "l=" << l << " cline length: " << cline.length() << " , length - match: " << cline.length() - matchString.length() << " match string: " << matchString << std::endl;
				if (cline[l] == matchString[matchIt])
				{
					//if (buggo) std::cout << "match" << std::endl;
					matchIt++;
					matching = true;
				}
				else if (matchIt < 18) {
					//if (buggo) std::cout << "no match" << std::endl;
					//if (matchIt > 10 && buggo) std::cout << "closed early on: " << cline << " with a matchit of " << matchIt << " and a l of " << l << " and a cline length of " << cline.length() << " and a matchit length of " << matchString.length() << std::endl;
					matchIt = 0;
					matching = false;
				}
				else {
					//the strings are already matched let the match go through
				}
			}
			if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
			if (matching && matchIt >= matchString.length())
			{
				//if (buggo) std::cout << "matched" << std::endl;
				//a geometry has been identified

				//create set of atoms for those inside
				std::vector<std::vector<atom>> set;
				std::vector<std::string> elements;
				//the next 4 lines do not contain atom info
				//if (buggo) std::cout << "skipping 4 lines" << std::endl;
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				char space = ' ';
				
				
				//if (buggo) std::cout << "reading an atom" << std::endl;

				while (getline(*gOutFile, cline) && cline[1] != '-')
				{
					std::string cnumber = "";

					//if (buggo) std::cout << cline << std::endl;
					int num = 0;
					std::string element;
					double x = 0;
					double y = 0;
					double z = 0;
					for (int i = 0; i < cline.length(); i++)
					{
						if (cline[i] == space)
						{
							if (cnumber.length() > 0)
							{
								//if (buggo) std::cout << "cnumber larger than one. it is " << cnumber << std::endl;
								num++;
								switch (num)
								{
								case 1:
									//center number, not sure what this is for. do nothing
									break;
								case 2:
									//atomic number
									element = atomicNumber(std::stoi(cnumber));
									break;
								case 3:
									//atomic type, not sure what this is. do nothing.
									break;
								case 4:
									//x coordinate
									x = std::stof(cnumber);
									break;
								case 5:
									//y
									y = std::stof(cnumber);
									break;
								case 6:
									//z
									z = std::stof(cnumber);
									break;
								default:
									break;

								}
								cnumber = "";
							}
							//else {
								//if (buggo) std::cout << "current char not large than one" << std::endl;
							//}

						}
						else {
							//std::cout << "adding " << std::to_string(cline[i]) << " or " << cline[i] << "? to cnumber" << std::endl;
							cnumber += cline[i];
						}
					}
					//figure out which set this belongs to.
					int index = -1;
					for (int i = 0; i < elements.size(); i++)
					{
						if (elements[i] == element)
						{
							index = i;
						}
					}
					//if (buggo) std::cout << "making atom for elemental index " << index << " with coords " << x << " " << y << " " << z << std::endl;
					if (index == -1)
					{
						//if (buggo) std::cout << "not yet indexed" << std::endl;
						//atom type not yet indexed
						elements.push_back(element);
						std::vector < atom> alist = { atom(x,y,z) };
						set.push_back(alist);
					}
					else {
						//if (buggo) std::cout << "pushing back" << std::endl;
						set[index].push_back(atom(x, y, z));
					}
				}
				//if (buggo) std::cout << "this structure is not the last one starting over with next" << std::endl;

				if (st != nullptr) delete st;
				st = new structure(set, elements);


			}
			//else if (buggo) std::cout << "no match" << std::endl;

		}
		if (cline.size() > eMatch.size())
		{
			bool matching = true;
			for (int m = 0; m < eMatch.size(); m++)
			{
				if (cline[m] != eMatch[m]) matching = false;
			}
			if (matching) energyLine = cline;
		}
			//structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool sort = false)
		
	}
	// Close the file
	//if there was no convergence then nullptr will be returned
	//if (buggo) std::cout << "closing" << std::endl;
	gOutFile->close();
	delete gOutFile;
	////
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));
	st->energy = energyHartree;
	/////
	//remove files
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm " + inFileName + "\n";
		std::string deleteFiles2 = "rm "+ outFileName + "\n";
		std::string deleteFiles3 = "rm -r core.*\n";

		gc = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			gc[i] = deleteFiles1[i];
		}
		std::cout << "command: " << gc << std::endl;
		std::system(gc);
		delete[] gc;


		gc = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			gc[i] = deleteFiles2[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

		gc = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			gc[i] = deleteFiles3[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;
		
	}

	return st;
}

structure* newOptimizationGaussian(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, int amscalculation = -1, int maxSCF = 300,int maxCycles = 100, double timeoutminutes = -1, double timeLeftMinutes = -1, bool keepFiles = false)
{
	std::string dirname = "calc" + std::to_string(amscalculation);
	char* ac;
	if (amscalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}

	std::string guassianTaskName = filename + " opt";

	std::string inFileName =  filename + ".com";
	std::string outFileName =  filename + ".log";

	writeToGuassian(s, guassianTaskName, dirname + "/" + inFileName, basis_name, method, charge, state, maxCycles, maxSCF);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + inFileName);

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + inFileName);
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + inFileName + "\n";
	commands += "nohup g16 <" + inFileName + " >" + outFileName + " &\n";
	std::cout << "printing input" << std::endl;


	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;

	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + outFileName);


	std::cout << "checking for gaussian completion" << std::endl;



	while (!getline(*aOutFile, checker))
	{
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);
	}
	std::cout << "Gaussian started" << std::endl;

	bool finishedOutput = false;
	bool finishedEnergy = false;
	bool finishedStructure = false;
	std::vector<std::string> optStrings = {};
	std::string outLine;
	std::string energyString;



	
	while (!finishedOutput || !finishedEnergy || !finishedStructure)
	{
		if (finishedOutput && !(finishedEnergy && finishedStructure))
		{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				std::system(ac);
				delete[] ac;
				return nullptr;

		}
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return nullptr;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return nullptr;
			}

		}

		optStrings = {};
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);

		std::string match = " Job cpu time:";
		std::string matchString = "Input orientation:";
		std::string eMatch = " SCF Done:  E";
		while (getline(*aOutFile, outLine)) {
			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
			}
			if (outLine.size() > matchString.size())
			{
				bool matching = true;
				for (int m = 0; m < matchString.size(); m++)
				{
					if (outLine[m] != matchString[m]) matching = false;
				}
				if (matching) finishedStructure = true;//cannot be made false
			}
			if (outLine.size() > eMatch.size())
			{
				bool matching = true;
				for (int m = 0; m < eMatch.size(); m++)
				{
					if (outLine[m] != eMatch[m]) matching = false;
				}
				if (matching) finishedEnergy = true;//cannot be made false
			}
		}

		
		
		

	}
	std::cout << "finished output: " << finishedOutput << std::endl;
	double seconds_since_start = difftime(time(0), start);


	aOutFile->close();
	delete aOutFile;
	std::ifstream* gOutFile = new std::ifstream(dirname + "/" + outFileName);
	std::string cline;
	structure* st = nullptr;

	std::string matchString = "Input orientation:";
	std::string eMatch = " SCF Done:  E";
	std::string energyLine = "";

	std::cout << "reading file" << std::endl;

	while (getline(*gOutFile, cline)) {
		if (cline.length() > matchString.length())
		{

			bool matching = false;
			int matchIt = 0;
			for (int l = 0; l < cline.length() - matchString.length(); l++)
			{
				if (cline[l] == matchString[matchIt])
				{
					matchIt++;
					matching = true;
				}
				else if (matchIt < 18) {
					matchIt = 0;
					matching = false;
				}
				else {
					//the strings are already matched let the match go through
				}
			}
			if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
			if (matching && matchIt >= matchString.length())
			{
				std::vector<std::vector<atom>> set;
				std::vector<std::string> elements;
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				char space = ' ';
				while (getline(*gOutFile, cline) && cline[1] != '-')
				{
					std::string cnumber = "";
					int num = 0;
					std::string element;
					double x = 0;
					double y = 0;
					double z = 0;
					for (int i = 0; i < cline.length(); i++)
					{
						if (cline[i] == space)
						{
							if (cnumber.length() > 0)
							{
								num++;
								switch (num)
								{
								case 1:
									//center number, not sure what this is for. do nothing
									break;
								case 2:
									//atomic number
									element = atomicNumber(std::stoi(cnumber));
									break;
								case 3:
									//atomic type, not sure what this is. do nothing.
									break;
								case 4:
									//x coordinate
									x = std::stof(cnumber);
									break;
								case 5:
									//y
									y = std::stof(cnumber);
									break;
								case 6:
									//z
									z = std::stof(cnumber);
									break;
								default:
									break;

								}
								cnumber = "";
							}
						}
						else {
							cnumber += cline[i];
						}
					}
					int index = -1;
					for (int i = 0; i < elements.size(); i++)
					{
						if (elements[i] == element)
						{
							index = i;
						}
					}
					if (index == -1)
					{
						elements.push_back(element);
						std::vector < atom> alist = { atom(x,y,z) };
						set.push_back(alist);
					}
					else {
						set[index].push_back(atom(x, y, z));
					}
				}
				if (st != nullptr) delete st;
				st = new structure(set, elements);


			}
		}
		if (cline.size() > eMatch.size())
		{
			bool matching = true;
			for (int m = 0; m < eMatch.size(); m++)
			{
				if (cline[m] != eMatch[m]) matching = false;
			}
			if (matching) energyLine = cline;
		}
	}
	gOutFile->close();
	delete gOutFile;

	//find the start of the energy
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));
	st->energy = energyHartree;

	std::cout << "Gaussian time: " << seconds_since_start << std::endl;
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm -r " + dirname + "\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


	}
	return st;
}
double newEnergyGaussian(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, int amscalculation = -1, int maxSCF = 300, int maxCycles = 100, double timeoutminutes = -1, double timeLeftMinutes = -1, bool keepFiles = false)
{
	std::string dirname = "calc" + std::to_string(amscalculation);
	char* ac;
	if (amscalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}

	std::string guassianTaskName = filename + " energy";

	std::string inFileName = filename + ".com";
	std::string outFileName = filename + ".log";

	writeToGuassian(s, guassianTaskName, dirname + "/" + inFileName, basis_name, method, charge, state,0, maxSCF);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + inFileName);

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + inFileName);
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + inFileName + "\n";
	commands += "nohup g16 <" + inFileName + " >" + outFileName + " &\n";
	std::cout << "printing input" << std::endl;


	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;

	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + outFileName);


	std::cout << "checking for gaussian completion" << std::endl;



	while (!getline(*aOutFile, checker))
	{
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);
	}
	std::cout << "Gaussian started" << std::endl;

	bool finishedOutput = false;
	bool finishedEnergy = false;
	std::vector<std::string> optStrings = {};
	std::string outLine;
	std::string energyString;



	std::string energyLine = "";
	while (!finishedOutput || !finishedEnergy)
	{
		if (finishedOutput && !(finishedEnergy))
		{

				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

		}
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;
			}

		}

		optStrings = {};
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);

		std::string match = " Job cpu time:";
		std::string eMatch = " SCF Done:  E";
		while (getline(*aOutFile, outLine)) {
			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
				else {
					//std::cout << match << "\n vs \n" << outLine << std::endl;
				}
			}
			if (outLine.size() > eMatch.size())
			{
				bool matching = true;
				for (int m = 0; m < eMatch.size(); m++)
				{
					if (outLine[m] != eMatch[m]) matching = false;
				}
				if (matching) {
					finishedEnergy = true;//cannot be made false
					energyLine = outLine;
				}
			}
		}


	}
	std::cout << "finished output: " << finishedOutput << std::endl;
	double seconds_since_start = difftime(time(0), start);


	aOutFile->close();
	delete aOutFile;
	
	//find the start of the energy
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));

	std::cout << "Gaussian time: " << seconds_since_start << std::endl;
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm -r " + dirname + "\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


	}
	return energyHartree;
}


std::vector<double> direction(structure& a, structure& b)
{
	//returns the direction between two structures
	std::vector<double> dir;
	for (int i = 0; i < a.elements.size(); i++)
	{
		for (int j = 0; j < a.set[i].size(); j++)
		{
			dir.push_back(b.set[i][j].x - a.set[i][j].x);
			dir.push_back(b.set[i][j].y - a.set[i][j].y);
			dir.push_back(b.set[i][j].z - a.set[i][j].z);
		}
	}

	return dir;
}

structure* getAcceptedStructure(double x, double y, double z, std::vector<std::string> elements, std::vector<int> composition, std::vector<structure*>& structures, std::vector<structure*>& localMinima, int coordinationSteps, std::vector<std::set<int>>& coordinationNumbers, int radialCriteriaPercent, double RIDthreshold )
{
	structure* retValue = nullptr;
	bool seedAccepted = false;
	if (DEBUG) std::cout << "looking for a new seed " << std::endl;
	while (!seedAccepted)
	{
		//cant delete s because it doesn't exist
		//delete s;
		retValue = new structure(x, y, z, composition, elements);
		std::vector<structure*>::iterator it = structures.begin();

		if (localMinima.size() >= coordinationSteps)
		{
			while (!coordinationNumber(*retValue, coordinationNumbers) || !radialCriteria(*retValue, radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria ";
				if (DEBUG) {
					if (coordinationNumber(*retValue, coordinationNumbers)) std::cout << " radial" << std::endl;
else std::cout << " coordination numbers" << std::endl;
				}
				delete retValue;
				retValue = new structure(x, y, z, composition, elements);
			}
		}
		else {
			while (!radialCriteria(*retValue, radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria (radial) in seed creation" << std::endl;
				delete retValue;
				retValue = new structure(x, y, z, composition, elements);
			}
		}

		//check that the new structure does not violate radius criteria
		if (DEBUG) std::cout << "checking for uniqueness " << std::endl;


		bool unique = true;

		while (it != structures.end() && unique)
		{
			double dist = RID(*retValue, *(*it));
			if (dist < RIDthreshold)
			{
				unique = false;
			}
			it++;
		}
		seedAccepted = unique;
		if (!seedAccepted)
		{
			delete retValue;
			//delete s here because we know we arent accepting it and we cant delete it at the top of the loop
		}
		else {
			if (DEBUG) std::cout << "seed structure accepted" << std::endl;
		}

		


		
	}
	return retValue;
}

void createSeeds()
{
	//this function will use symmetry to create seeds

}

double proportion(double energy,double eThr)
{
	/*
	TODO, should be a sigmoid function
	*/
	//double eThr = 0;
	//eThr is the threshold for energy change that should be considered a complete movement
	//ranges from 1 to -1 based on x
	//the change in energy is negative if good, but we want this scalar to always be positive, and eThr/e normalization input is positive.
	//therefore
	double eScalar = 0;
	if (energy > 0)	eScalar = energy / eThr;
	else eScalar = -energy / eThr;

	/*
	MAJOR CHANGE to fix an error hopefully

	if the change in energy is negative this is good so we want the proportion to be positive.
	if the change in energy is positive this is bad so we want the proportion to be negative.

	*/
	eScalar = -energy / eThr;
	//a negative sign was placed here as well infront of escalar
	double sigmoid = 2 / (1 + exp((-eScalar)*2*exp(1))) - 1;
	
	return sigmoid;
}


std::pair<structure*,std::vector<double>> createXnew(structure* s, double percentChangeI, double percentChangeF, int step, double x, double y, double z, std::vector<double>& directionVector, bool blindGradient = false, int exponent1 = 0, int exponent2 = 0,  double deltaE = 0,double eThr = 1, int hybridSteps = 0, int ccSteps = 1)
{
	//direction vector must be size s.natoms*3 
	std::vector<double> newDirVector = {};
	

	if (blindGradient && ccSteps >= hybridSteps)
	{
		//new blind gradient descent method

		/*
		Negative last step: (previous direction)*(proportionality scalar) + abs(previous direction)*(random vector 1 to -1)*(proportionality scalar)^n
		Positive last step: (previous direction)*(proportionality scalar) + abs(previous direction)*(random vector 1 to -1)*(1 - proportionality scalar)^n
		*/
		double p = proportion(deltaE,eThr);
		//if (DEBUG) std::cout << "VECTORS: before and after vectors with proportion: " << p << " from deltaE: " << deltaE << std::endl;
		//if (DEBUG) for (int i = 0; i < directionVector.size(); i++) std::cout << directionVector[i] << " ";
		//if(DEBUG) std::cout << std::endl;

		/*
				With the MAJOR CHANGE in proportion, in the case of positive energy in the past it was ok for (1-p) to be positive as the signs would randomly change anyways, but now we want it to be 1+p as p should be negative
				in proportion now in this case p is negative so we want for the previous 1-p term to be negativ
		*/

		/*//old way
		if (deltaE > 0)
		{
			for (int i = 0; i < directionVector.size(); i++)
			{
				
				if (directionVector[i] > 0) newDirVector.push_back(directionVector[i] * p + pow(directionVector[i], exponent1) * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * pow((1-p), exponent2));
				else newDirVector.push_back(directionVector[i] * p + pow(-directionVector[i], exponent1) * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * pow((1 - p), exponent2));
			}
		}
		else {
			for (int i = 0; i < directionVector.size(); i++)
			{
				if (directionVector[i] > 0) newDirVector.push_back( directionVector[i] * p + pow(directionVector[i], exponent1) * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * pow((p), exponent2));
				else newDirVector.push_back( directionVector[i] * p + pow(-directionVector[i], exponent1) * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * pow((p), exponent2));
			}
		}
		*/

		if (deltaE > 0)
		{
			for (int i = 0; i < directionVector.size(); i++)
			{
				newDirVector.push_back(directionVector[i] * p + pow(directionVector[i], exponent1) * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * pow((1 + p), exponent2));
			}
		}
		else {
			for (int i = 0; i < directionVector.size(); i++)
			{
				newDirVector.push_back(directionVector[i] * p + pow(directionVector[i], exponent1) * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * pow((p), exponent2));
			}
		}
		//if (DEBUG) for (int i = 0; i < directionVector.size(); i++) std::cout << directionVector[i] << " ";
		//if (DEBUG) std::cout << std::endl;
		
	}
	else {
		//old basin hopping method
		//create the movement vector
		//get number of atoms first
		//movement parameters are based on sigmoid variation from step on x,y,z maximum change.
		//maybe later this should be changed to based the slope of the gradient which can be determiend from the proportion of change in energy to last movement.
		//std::vector<double> movement;
		double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));


		for (int i = 0; i < s->elements.size(); i++)
		{
			for (int j = 0; j < s->set[i].size(); j++)
			{
				//should this be randomized between -1 to 1?
				//x, y , z are the maximum movement, defining those parameters.
				//movement.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * x * sigmoidPercentageFromStep / 100);
				//movement.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * y * sigmoidPercentageFromStep / 100);
				//movement.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * z * sigmoidPercentageFromStep / 100);
				newDirVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * x * sigmoidPercentageFromStep / 100));
				newDirVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * y * sigmoidPercentageFromStep / 100));
				newDirVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * z * sigmoidPercentageFromStep / 100));
			}

		}
		//direction vec

	}
	return std::pair<structure*, std::vector<double>>(new structure(*s, newDirVector), newDirVector);
}


std::pair< std::vector<int>, std::vector<std::string>> getComposition(std::string comp)
{
	bool bugging = false;
	std::vector<int> composition;
	std::vector<std::string> elements;
	bool instring = true;//assume we start in a string
	std::string current = "";
	int currint = 0;
	for (int i = 0; i < comp.size(); i++)
	{
		if (bugging) std::cout << comp[i] << std::endl;
		if (comp[i] >= 48 && comp[i] <= 57)
		{
			//numeric
			if (i > 0) {
				if (instring) {
					if (bugging) std::cout << "added " << current << " to elements" << std::endl;
					//switch
					instring = false;
					elements.push_back(current);
					current = "";
					
				}
			}
			currint = currint * 10 + (comp[i] - 48);
		}
		if (comp[i] >= 65 && comp[i] <= 90)
		{
			//uppercase
			if (i > 0) {
				if (!instring) {
					if (bugging) std::cout << "added " << currint << " to composition" << std::endl;
					//switch
					instring = true;
					composition.push_back(currint);
					currint = 0;
					
				}
				else {
					if (bugging) std::cout << "added " << current << " to elements and " << 1 << " to composition" << std::endl;
					//elements begin with capitals so the last one is over
					elements.push_back(current);
					composition.push_back(1);//only one of the last element
					current = "";
					
				}
			}
			current = current + comp[i];
		}
		if (comp[i] >= 97 && comp[i] <= 122)
		{
			//lowercase
			if (i > 0) {
				if (!instring) {
					if (bugging) std::cout << "added " << currint << " to composition" << std::endl;
					//switch
					instring = true;
					composition.push_back(currint);
					currint = 0;
					
				}
			}
			current = current + comp[i];
		}
	}
	//dont forget who was added last
	if (!instring) {
		if (bugging) std::cout << "added " << currint << " to composition" << std::endl;
		composition.push_back(currint);
		currint = 0;
		
	}
	else {
		if (bugging) std::cout << "added " << current << " to elements and " << 1 << " to composition" << std::endl;
		//elements begin with capitals so the last one is over
		elements.push_back(current);
		composition.push_back(1);//omly one of the last element
		current = "";
		
	}
	return std::pair< std::vector<int>, std::vector<std::string>>(composition, elements);
}
//(structure* initializationStructure, double x, double y, double z, double percentChangeI, double percentChangeF, double RIDpercent, std::vector<int> composition, std::vector<std::string> elements, int charge, int state, char calculator, std::string method, std::string basis, int optimizationIterations, int coordinationSteps, int HFsteps, int time, int localMinimaTrapping, int radialCriteriaPercent, int radialCriteriaTrappingIteration, int radialCriteriaTrappingFreedom, std::string outputFileName, std::string taskName, double a, double b, double c, bool blindGradientDescent = false, int directionExponent = 0, int proportionExponent = 1, double deltaEThr = 0.01, int hybridCycles = 0)

structure* optimize(structure& s)
{
	/*
	Todo
	*/
	//thsi causes a crash dont use!
	return nullptr;
}

double energyCalculation(structure* s,std::string basis, std::string method,int charge, int state, char calculator, std::string taskName, int energyCall, bool keepFiles, std::string converge = "1.0e-6", double timeLeftMinutes = -1, double timeoutMinutes = -1, std::string AMSHOME = "/home/jpburkhardt/scm/ams2023.104") {
	double retValue = 0;
	switch (calculator)
	{
	case 'a':
		//adf/scm energy
		retValue = energyADF(*s, taskName, basis, method, charge, state, AMSHOME,energyCall,converge, timeoutMinutes, timeLeftMinutes, keepFiles);
		break;
	case 'b':
		//basin hopping energy
		retValue = basinHoppingEnergy(*s, energyCall++, taskName, 0, 0, 0);
		break;
	case 'g':
		//gaussain energy
		retValue = newEnergyGaussian(*s, taskName, basis, method, charge, state, energyCall, 300, 0, timeoutMinutes, timeLeftMinutes, keepFiles);
		//retValue = energyGuassian(*s, energyCall++, taskName, basis, method, charge, state, keepFiles);
		break;
	default:
		retValue = energyGuassian(*s, energyCall++, taskName, basis, method, charge, state, keepFiles);
		break;
	}
	return retValue;
}

structure* optimizationCalculation(char calculator, structure* Xnew, int step, std::string taskName, int optimizationIterations,int maxSCF, double x, double y, double z, std::vector<double> & optimizationMovement,std::vector<structure*>& structures,std::string basis, std::string method, int charge, int state, bool keepFiles, std::string converge = "1.0e-6",double timeLeftMinutes = -1,double timeoutMinutes = -1, std::string AMSHOME = "/home/jpburkhardt/scm/ams2023.104")
{
	structure* optimizedStructure;
	switch (calculator)
	{
		case 'a':
		{
			optimizedStructure = optimizationADF(*Xnew, taskName, basis, method, charge, state, AMSHOME,step,maxSCF,converge,timeoutMinutes,timeLeftMinutes, keepFiles);
			break;
		}
		case 'b':
		{
			optimizedStructure = optimizationBasinHopping(*Xnew, step, taskName, optimizationIterations, x, y, z);
			break;
		}
		case 'g':
		{
			optimizedStructure = newOptimizationGaussian(*Xnew, taskName, basis, method, charge, state, step, maxSCF, optimizationIterations, timeoutMinutes, timeLeftMinutes, keepFiles);
			//optimizedStructure = optimizationGaussian(*Xnew, step, taskName, basis, method, optimizationIterations, charge, state, maxSCF, keepFiles);
			break;
		}
		default:
		{
			optimizedStructure = optimizationGaussian(*Xnew, step, taskName, basis, method, optimizationIterations, charge, state,maxSCF, keepFiles);
			break;
		}
	}
	if (optimizedStructure == nullptr)
	{
		if (DEBUG) std::cout << "failed to optimize" << std::endl;
		//optimizedStructure = Xnew;//give it back Xnew with a terrible energy
		//optimizedStructure->energy = DBL_MAX;
	}
	else {
		optimizationMovement = direction(*Xnew, *optimizedStructure);
	}
	return optimizedStructure;
}

void groupTheoryGradientDescent(std::string previousRun, std::vector<int> composition, std::vector<std::string> elements, double x, double y, double z, double rlimit, char seedMode, double percentChangeI, double percentChangeF, double RIDpercent, int charge, int state, char calculator, std::string method, std::string basis, int optimizationIterations, int SCFsteps, int HFsteps, double time, int localMinimaTrapping, int uphillstopping, int misstepTrapping, int radialCriteriaPercent, std::string outputFileName, std::string taskName, int directionExponent = 0, int proportionExponent = 1, double deltaEThr = 1, bool checkEnergy = false, bool optimizePositive = true, bool keepFiles = false, bool filterOff = false, std::string converge = "1.0e-6", int timeOutMinutes = -1, int timeOutMinutesSP = -1, bool covalentCriteriaMode = false, double covalentCriteriaPercent = 0, int covalentCriteriaLimit = 1000, int bondRequirement = 0, int criteriaIterations = 10000, bool planar = false, bool partialSymmetry = false, bool extraHoles = false, bool alignHoles = true,bool variantSymmetry = false,double nMultiplier = 1, double temperature = 293.15)
{
	auto start = std::chrono::high_resolution_clock::now();

	if (planar)
	{
		std::cout << "planar functionality not yet implemented" << std::endl;
	}

	int energyCall = 0;
	double RIDthreshold = 0;
	std::vector<structure*> structures = {};
	//std::vector<std::set<int>> coordinationNumbers;
	structure* gmeS = nullptr;//the global minima structure known
	double globalMinimumEnergy = 0;//energy corresponding to gmeS
	double globalMinimumTime = 0;
	
	int totalStructuresSearchedAtGME = 0;

	structure* ls = nullptr; //current/last structure explored
	double lastEnergy = 0;//energy corresponding to ls

	structure* cMinima;//current local minima since the last seed reset
	double latestLowestEnergy;//energy corresponding to cminima
	std::vector<double> lowMovementVecor;//corresponding to cminima
	double cMinDeltaE;//delta E at cmin

	int stepsSinceNew = 0;
	int recyclesSinceNew = 0;
	int upHillSteps = 0;
	int optimizationsInCurrentSeed = 0;

	std::vector<structure*> localMinima = {};
	int totalStructuresSearched = 0;

	//not used in continue
	structure* s = nullptr;
	double seedEnergy = 0;
	int ccSteps = 0;//current chain steps

	int calculations = 0;
	int step = 0;

	std::vector<std::pair<int, bool>> symmetries = {};//used in seed mode i & r
	std::vector<structure*> seeds = {};//used in seed mode o
	if(seedMode != 'o')  symmetries = possibleSymmetry(composition, elements, rlimit, radialCriteriaPercent, 1, criteriaIterations);
	int seedNo = 0;

	if (previousRun != "")
	{
		if (STDOUT) std::cout << "continuing from previous file_: " << previousRun << std::endl;
		energyCall = 1;
		std::ifstream startFile(previousRun);
		std::string inLine;
		std::vector<int> composition;
		std::vector<std::string> elements;
		int stepsSinceNew = 0;
		std::getline(startFile, inLine);// "Composition:\n";
		std::getline(startFile, inLine);
		std::pair< std::vector<int>, std::vector<std::string>> comp = getComposition(inLine);
		int debugo = 0;
		composition = comp.first;
		elements = comp.second;
		std::getline(startFile, inLine); //"threshold:\n";
		std::getline(startFile, inLine);
		RIDthreshold = std::stoi(inLine);
		std::getline(startFile, inLine);// best structure
		std::getline(startFile, inLine);
		gmeS = new structure(inLine, elements, composition);
		globalMinimumEnergy = gmeS->energy;//a variable for monitoring the global minimum found
		globalMinimumTime = 0;
		std::getline(startFile, inLine);//# minutes searched prior
		int priorEnd = 0;
		for (int i = 0; i < inLine.size(); i++)
		{
			if (inLine[i] == ' ' && priorEnd == 0) priorEnd = i;
		}
		globalMinimumTime = std::stof(inLine.substr(0, priorEnd));
		std::getline(startFile, inLine);//# minutes searched prior
		priorEnd = 0;
		for (int i = 0; i < inLine.size(); i++)
		{
			if (inLine[i] == ' ' && priorEnd == 0) priorEnd = i;
		}
		totalStructuresSearchedAtGME = std::stof(inLine.substr(0, priorEnd));
		std::cout << "Current global minimum structure found after " << totalStructuresSearchedAtGME << " structures with energy " << globalMinimumEnergy << ": " << std::endl;
		gmeS->print();
		std::getline(startFile, inLine);//energy of best structure

		for (int i = 0; i < composition.size(); i++) for (int j = 0; j < composition[i]; j++) std::getline(startFile, inLine);//element followed by coordinates
		std::getline(startFile, inLine);// "last energy,structure:\n";
		std::getline(startFile, inLine);
		ls = new structure(inLine, elements, composition);
		lastEnergy = ls->energy;
		std::getline(startFile, inLine);//steps since new
		std::getline(startFile, inLine);
		stepsSinceNew = std::stoi(inLine);
		std::getline(startFile, inLine);//recycles since new
		std::getline(startFile, inLine);
		recyclesSinceNew = std::stoi(inLine);
		std::getline(startFile, inLine);//uphill steps
		std::getline(startFile, inLine);
		upHillSteps = std::stoi(inLine);
		std::getline(startFile, inLine);//uphill steps
		std::getline(startFile, inLine);
		optimizationsInCurrentSeed = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		step = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		ccSteps = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		seedNo = std::stoi(inLine);
		std::getline(startFile, inLine);
		seeds = {};
		std::getline(startFile, inLine);
		while (inLine[0] != 'l')
		{
			structure* next = new structure(inLine, elements, composition);
			seeds.push_back(next);
			std::getline(startFile, inLine);
		}
		/*
		Todo set steps since new
		*/
		std::getline(startFile, inLine);// "local minima energy,structure:\n";
		localMinima = {};
		cMinima = nullptr;
		std::getline(startFile, inLine);
		while (inLine[0] != 'a')
		{
			structure* next = new structure(inLine, elements, composition);
			localMinima.push_back(next);
			cMinima = next;//constantly reset until we get the actual last one

			//load the next line. the reason this isn't done in the conditional is because we need to stop at a text line
			std::getline(startFile, inLine);
		}
		//std::getline(startFile, line);// "all structures energy,structure:\n";// this should already be ignored
		//get all the structures now
		latestLowestEnergy = cMinima->energy;
		while (std::getline(startFile, inLine))
		{
			structure* next = new structure(inLine, elements, composition);
			structures.push_back(next);
		}
		totalStructuresSearched = structures.size();
		// Close the file
		startFile.close();
	}
	else {
		RIDthreshold = threshold(RID, composition, elements, rlimit, rlimit, rlimit, 1, x, y, z, RIDpercent);
		stepsSinceNew = 0;
		upHillSteps = 0;
		recyclesSinceNew = 0;//NEW METHOD - now steps since new will trigger returns to the last local minima, while recycles since new will trigger new seeds. monitored by the limit missteptrapping
		optimizationsInCurrentSeed = 0;
		structures = {};
		localMinima = {};
		totalStructuresSearched = structures.size();

		if (seedMode != 'o') {
			do {
				//s = seedFromGroupTheory(composition, elements, rlimit);
				s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);
				while (!radialCriteria(*s, radialCriteriaPercent))
				{
					if (DEBUG) std::cout << "did not pass criteria (radial) in seed creation" << std::endl;
					delete s;
					s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);
				}
			} while (!bondingRequirement(s, radialCriteriaPercent, bondRequirement));
			energyCall = 0;
			double timeLeftMinutes = (double)time * 60;

			seedEnergy = energyCalculation(s, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
		}
		else {
			//create all the seeds, check all their energies, order them by energy and go one by one
			std::cout << "Creating Seeds, time(m): " << (double)time * 60 << std::endl;
			seeds = getSeeds(composition, elements, rlimit, radialCriteriaPercent, bondRequirement, criteriaIterations, partialSymmetry, extraHoles, alignHoles, variantSymmetry, nMultiplier);
			std::cout << "Calculating energy of " << seeds.size() << " seeds, time(m) : " << (double)time * 60 << std::endl;
			for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
			{
				double timeLeftMinutes = (double)time * 60;
				double currentEnergy = energyCalculation(*it, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
				(*it)->energy = currentEnergy;
			}
			//sort seeds by energy
			std::cout << "Sorting seeds, time(m): " << (double)time * 60 << std::endl;
			std::sort(seeds.begin(), seeds.end(),energyCompare);
			std::ofstream outputfile(outputFileName + "_seed_energies" + ".txt");
			// Write to the file
			for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
			{
				std::string line = "energy: " + std::to_string((*it)->energy) + "\n";
				outputfile << line;
			}

			// Close the file
			outputfile.close();
			writeToXyz(seeds, outputFileName + "_seeds" + ".xyz");

			//start with seeed 1
			s = seeds[0];

		}
		gmeS = s;//the global minima structure known
		globalMinimumEnergy = seedEnergy;//energy corresponding to gmeS
		totalStructuresSearchedAtGME = 1;

		ls = s; //current/last structure explored
		lastEnergy = seedEnergy;//energy corresponding to ls

		cMinima = s;//current local minima since the last seed reset
		latestLowestEnergy = seedEnergy;//energy corresponding to cminima
		

		ccSteps = 0;//current chain steps


		step = 0;



	}
	
	
	cMinDeltaE = 0;//delta E at cmin

	/*
	Create seed
	*/

	std::vector<double> movement = {};
	int natoms = 0;
	for (int i = 0; i < ls->elements.size(); i++)
	{
		for (int j = 0; j < ls->set[i].size(); j++)
		{
			natoms++;
			movement.push_back(x * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
			movement.push_back(y * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
			movement.push_back(z * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
		}
	}
	lowMovementVecor = movement;//corresponding to cminima//data is lost in continue
	
	

	
	
	//also initialize the energy changel

	double deltaE = -0.1 * deltaEThr;
	//change 10

	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double minutes = (double)duration.count() / (double)60;
	if (DEBUG) std::cout << "minutes:  " << minutes << " ; time limit: " << double(time) * 60 - 30 << std::endl;
	int maxint_ = 0;//remove later!

	while (minutes < double(time) * 60 - 30)
	{
		ccSteps++;//current chain steps. gets reset every time a new seed is found, and is used to determine when to switch to blind gradient descent in the hybrid model. unfortunately it is distinct from all other step counters.
		maxint_++;
		if (DEBUG) std::cout << "MONITOR______________: " << " totalStructures: " << structures.size() << " total steps: " << maxint_ << std::endl;
		if (DEBUG || STDOUT) std::cout << "step: " << step << " time: " << minutes / 60 << std::endl;
		//step += 1;
		structure* Xnew = nullptr;

		//move back to cMinima
		if (stepsSinceNew > localMinimaTrapping || upHillSteps > uphillstopping)
		{
			//NEW METHOD
			//we will return to the last local minima
			ls = cMinima;
			lastEnergy = latestLowestEnergy;
			movement = lowMovementVecor;

			recyclesSinceNew++;
			stepsSinceNew = 0;

			//reset movement
			movement = {};
			for (int i = 0; i < ls->elements.size(); i++)
			{
				for (int j = 0; j < ls->set[i].size(); j++)
				{
					movement.push_back(x * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
					movement.push_back(y * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
					movement.push_back(z * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
				}
			}
			//make delta Ethr small at first
			deltaE = -0.1 * deltaEThr;
		}
		//Find a new seed
		if (recyclesSinceNew > misstepTrapping)
		{
			seedNo++;
			movement = {};
			for (int i = 0; i < ls->elements.size(); i++)
			{
				for (int j = 0; j < ls->set[i].size(); j++)
				{
					movement.push_back(x * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
					movement.push_back(y * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
					movement.push_back(z * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
				}
			}
			bool seedAccepted = false;
			if (seedMode == 'o')
			{
				if (seedNo < seeds.size()) {
					s = seeds[seedNo];
					seedEnergy = s->energy;
				}
				else {
					seedNo = 0;
					//seedMode = 'n';//we ran through all the seeds, go to no symmetry for continued running
					/*
					* TODO
					There is currently an error where all the seeds made with no symmetry are broken, while this is not fixed we should loop the old seeds.
					*/
				}
			}
			if (seedMode != 'o') {
				while (!seedAccepted) {
					std::vector<structure*>::iterator it = structures.begin();
					do {
						//s = seedFromGroupTheory(composition, elements, rlimit);
						s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);
						while (!radialCriteria(*s, radialCriteriaPercent))
						{
							if (DEBUG) std::cout << "did not pass criteria (radial) in seed creation" << std::endl;
							delete s;
							s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);
							//s = seedFromGroupTheory(composition, elements, rlimit);
						}
					} while (!bondingRequirement(s, radialCriteriaPercent, bondRequirement));
					if (DEBUG) std::cout << "checking for uniqueness " << std::endl;
					bool unique = true;
					if (!filterOff) {
						while (it != structures.end() && unique)
						{
							double dist = RID(*s, *(*it));
							if (dist < RIDthreshold)
							{
								unique = false;
							}
							it++;
						}
					}
					seedAccepted = unique;
				}
				double timeLeftMinutes = (double)time * 60 - minutes;
				seedEnergy = energyCalculation(s, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
			}
			


			deltaE = -0.1 * deltaEThr;
			recyclesSinceNew = 0;
			stepsSinceNew = 0;
			optimizationsInCurrentSeed = 0;
			

			//set the last structure to the new seed
			ls = s;
			lastEnergy = seedEnergy;

			//add the last minima to list
			localMinima.push_back(cMinima);
			

			//set the current local minima to seed as well
			latestLowestEnergy = seedEnergy;
			cMinima = s;

		}
		if (STDOUT) std::cout << "Creating a new structure" << std::endl;
		std::pair<structure*, std::vector<double>> XnewPair = createXnew(ls, percentChangeI, percentChangeF, step, x, y, z, movement, true, directionExponent, proportionExponent, deltaE, deltaEThr, 0, ccSteps);
		if (covalentCriteriaMode)
		{
			int attempts = 0;
			while ((!covalentCriteria(XnewPair.first, covalentCriteriaPercent) || !bondingRequirement(XnewPair.first, covalentCriteriaPercent, bondRequirement)) && attempts < covalentCriteriaLimit)
			{
				delete XnewPair.first;
				XnewPair = createXnew(ls, percentChangeI, percentChangeF, step, x, y, z, movement, true, directionExponent, proportionExponent, deltaE, deltaEThr, 0, ccSteps);
				attempts++;
			}
			if (attempts >= covalentCriteriaLimit)
			{
				std::cout << "BGD algorithm trapped by covalent criteria. Returning to LLM with minor inefficency" << std::endl;
				//group theory gradient is stuck, return to the last local minima
				XnewPair = std::pair<structure*, std::vector<double>>(cMinima, lowMovementVecor);
				XnewPair.first->energy = latestLowestEnergy;

			}
		}
		Xnew = XnewPair.first;
		std::vector<double> proposedDirectionVector = XnewPair.second;
		structures.push_back(Xnew);

		movement = proposedDirectionVector;
		std::vector<double> optimizationMovement;
		if (STDOUT)
		{
			std::cout << "New structure before optimization:" << std::endl;
			Xnew->print();
		}


		bool optimize = true;
		if (checkEnergy)
		{
			
			//prevent optimization on terrible structures
			double timeLeftMinutes = (double)time * 60 - minutes;

			double tempEnergy = energyCalculation(Xnew, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles,converge, timeLeftMinutes, timeOutMinutesSP);
			if (!optimizePositive) {
				if (tempEnergy < 0)
				{
					optimize = true;
				}
				else {
					optimize = false;
				}
			}
			Xnew->energy = tempEnergy;
		}
		if (optimizationIterations > 0 && optimize) {
			if (STDOUT) std::cout << "optimizing current new structure" << std::endl;
			
			std::string cmethod = "HF";
			if (optimizationsInCurrentSeed >= HFsteps) cmethod = method;

			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
			minutes = (double)duration.count();
			double timeLeftMinutes = (double)time * 60 - minutes;

			structure* tempXnew = optimizationCalculation(calculator, Xnew,++calculations, taskName, optimizationIterations, SCFsteps, x, y, z, optimizationMovement, structures, basis, cmethod, charge, state,keepFiles,converge, timeLeftMinutes,timeOutMinutes);
			if (tempXnew == nullptr)
			{
				if (STDOUT) std::cout << "optimization failed" << std::endl;
				std::cout << "optimization failed.";
				if (checkEnergy) std::cout << "Energy before optimization (used) : " << Xnew->energy;
				std::cout << std::endl;
				if(!checkEnergy) Xnew->energy = DBL_MAX;//as of now energy is not taken from a failed optimization, the failed optimization structure is also not used but in future versions this will be implemented. so the energy should be set to double max for the failure, therefore with check energy turned on the algorithm will be more proficient
			}
			else {
				Xnew = tempXnew;
				structures.push_back(Xnew);
				std::vector<double> combinedMovement;
				if (STDOUT) std::cout << "movement before optimization: ";
				//print the movement and optimization movement vectors
				if (STDOUT) for (int i = 0; i < movement.size(); i++) std::cout << " " << movement[i];
				if (STDOUT) std::cout << std::endl << "optimization movement: ";
				if (STDOUT) for (int i = 0; i < optimizationMovement.size(); i++) std::cout << " " << optimizationMovement[i];
				if (STDOUT) std::cout << std::endl << "combined: ";
				std::cout << std::endl << "movement size:" << movement.size() << " " << "opt size: " << optimizationMovement.size() << std::endl;
				for (int i = 0; i < movement.size(); i++) combinedMovement.push_back(movement[i] + optimizationMovement[i]);
				std::cout << std::endl << "movement size:" << movement.size() << " " << "opt size: " << optimizationMovement.size() << std::endl;
				if (STDOUT) for (int i = 0; i < combinedMovement.size(); i++) std::cout << " " << combinedMovement[i];
				if (STDOUT) std::cout << std::endl;
				movement = combinedMovement;
				if (STDOUT)
				{
					std::cout << "New structure after optimization with energy:"  << Xnew->energy  << std::endl;
					Xnew->print();
					
				}

			}

		}


		if (DEBUG) std::cout << "calling energy of Xnew " << std::endl;
		double newEnergy = Xnew->energy;//energy is now taken from the optimizers
			//energyCalculation(Xnew, basis, method, charge, state, calculator, taskName, energyCall, keepFiles);
		if (DEBUG) std::cout << "Xnew energy: " << newEnergy << std::endl;
		if (STDOUT) std::cout << "New structure energy: " << newEnergy << std::endl;
		if (STDOUT) Xnew->printZmatrix();

		totalStructuresSearched += 1;//only add 1 as we only checked the optimized energy, not the xnew energy in that case
		if (globalMinimumEnergy > newEnergy)
		{
			globalMinimumEnergy = newEnergy;
			globalMinimumTime = (double)duration.count();
			gmeS = Xnew;
			totalStructuresSearchedAtGME = totalStructuresSearched;
		}

		Xnew->energy = newEnergy;

		//Checking criteria
		if (newEnergy < lastEnergy)
		{
			upHillSteps = 0;
			if (STDOUT) std::cout << "BGG Mode: New structure accepted for lower energy, set new local minima " << std::endl;


			if (newEnergy < latestLowestEnergy)
			{
				latestLowestEnergy = newEnergy;
				stepsSinceNew = 0;
				
				cMinima = Xnew;
				lowMovementVecor = movement;
				cMinDeltaE = newEnergy - lastEnergy;
			}
			else {
				//not cminima, steps since new increasing anyways
				stepsSinceNew++;
			}
		}
		else
		{
			upHillSteps++;
			stepsSinceNew++;
			if (STDOUT) std::cout << "BGG Mode: New structure has higher energy, steps since new local Minima increased: " << stepsSinceNew << " with " << localMinimaTrapping << " limit" << std::endl;

		}

		//we will always be accepting
		ls = Xnew;
		deltaE = newEnergy - lastEnergy;
		lastEnergy = newEnergy;
		step++;




		if (STDOUT) {
			std::cout << "CYCLE END: Current global minima energy, time,  and structure: " << globalMinimumEnergy << " " << globalMinimumTime << std::endl;
			gmeS->print();
			gmeS->printZmatrix();
		}



		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutes = (double)duration.count();
	}

	localMinima.push_back(cMinima);



	if (DEBUG || STDOUT) std::cout << "printing results " << std::endl;
	gmeS->writeToXyz(outputFileName + ".xyz");
	//gmeS->writeToObj(outputFileName + ".obj");
	/*
	print all structures
	*/
	std::string outputString = "";
	//add composition
	outputString += "Composition:\n";
	for (int i = 0; i < composition.size(); i++)
	{
		outputString += elements[i] + std::to_string(composition[i]);
	}
	outputString += "\n";

	//rid threshold
	outputString += "RID threshold:\n";
	outputString += std::to_string(0) + "\n";
	//add the Global minima structure
	outputString += "Best structure:\n";
	std::string line = "energy: " + std::to_string((gmeS)->energy);
	for (int e = 0; e < (gmeS)->set.size(); e++)
	{
		for (int a = 0; a < (gmeS)->set[e].size(); a++)
		{
			//each atom
			line += " " + (gmeS)->elements[e] + " " + std::to_string((gmeS)->set[e][a].x) + " " + std::to_string((gmeS)->set[e][a].y) + " " + std::to_string((gmeS)->set[e][a].z);
		}
	}
	line += "\n";
	outputString += line;
	
	std::string gline = std::to_string((int)globalMinimumTime) + " minutes prior\n" + std::to_string(totalStructuresSearchedAtGME) + " structures searched prior" + "\n" + "Energy: " + std::to_string(globalMinimumEnergy) + " \n";
	int zitt = 0;
	gmeS->makeZmatrix();
	for (int e = 0; e < gmeS->set.size(); e++)
	{
		for (int a = 0; a < gmeS->set[e].size(); a++)
		{
			//each atom
			gline += " " + gmeS->elements[e] + " " + std::to_string(gmeS->zmatrix[zitt][0]) + " " + std::to_string(gmeS->zmatrix[zitt][1]) + " " + std::to_string(gmeS->zmatrix[zitt][2]) + "\n";
			zitt++;
		}
	}
	/*
	for (int e = 0; e < gmeS->set.size(); e++)
	{
		for (int a = 0; a < gmeS->set[e].size(); a++)
		{
			//each atom
			gline += " " + gmeS->elements[e] + " " + std::to_string(gmeS->set[e][a].x) + " " + std::to_string(gmeS->set[e][a].y) + " " + std::to_string(gmeS->set[e][a].z) + "\n";
		}
	}*/

	outputString += gline;
	//add the last accepted structure
	outputString += "last energy,structure:\n";
	line = std::to_string(s->energy);
	for (int e = 0; e < s->set.size(); e++)
	{
		for (int a = 0; a < s->set[e].size(); a++)
		{
			//each atom
			line += " " + s->elements[e] + " " + std::to_string(s->set[e][a].x) + " " + std::to_string(s->set[e][a].y) + " " + std::to_string(s->set[e][a].z);
		}
	}
	line += "\n";
	outputString += line;
	//add the most recent steps since new
	outputString += "steps since new:\n";
	outputString += std::to_string(stepsSinceNew) + "\n";
	outputString += "recycles since new:\n";
	outputString += std::to_string(recyclesSinceNew) + "\n";
	outputString += "uphill steps since new:\n";
	outputString += std::to_string(upHillSteps) + "\n"; 
	outputString += "optimizations in current seed:\n";
	outputString += std::to_string(optimizationsInCurrentSeed) + "\n";
	outputString += "step:\n";
	outputString += std::to_string(step) + "\n";
	outputString += "ccStep:\n";
	outputString += std::to_string(ccSteps) + "\n";
	//seedstructures with energies NEW, incorporate back into the conitnue function
	outputString += "seedNo:\n";
	outputString += std::to_string(seedNo) + "\n";
	outputString += "seeds, energy structure:\n";
	for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
	{
		std::string line = "energy: " + std::to_string((*it)->energy);
		for (int e = 0; e < (*it)->set.size(); e++)
		{
			for (int a = 0; a < (*it)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*it)->elements[e] + " " + std::to_string((*it)->set[e][a].x) + " " + std::to_string((*it)->set[e][a].y) + " " + std::to_string((*it)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	outputString += "local minima only: energy,structure:\n";
	//Old
	//int iii = 0;
	std::cout << "local minima size: " << localMinima.size() << std::endl;
	for (std::vector<structure*>::iterator c_structure = localMinima.begin(); c_structure != localMinima.end(); c_structure++)
	{

		//std::cout << "going through the localMinima " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}

	outputString += "all structures: energy,structure:\n";
	//iii = 0;
	for (std::vector<structure*>::iterator c_structure = structures.begin(); c_structure != structures.end(); c_structure++)
	{
		//for each structure
		//std::cout << "going through the all structures " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	int numStructures = structures.size();
	std::cout << "output filename: " << outputFileName + "_" + std::to_string(numStructures) + ".txt" << std::endl;
	std::ofstream outputfile(outputFileName + "_" + std::to_string(numStructures) + ".txt");

	// Write to the file
	outputfile << outputString;

	// Close the file
	outputfile.close();

	/*
	Create an xyz file with all local minima structures, and also a corresponding file with energies and proportion from partition function;
	*/
	std::ofstream localMinimaEnergyFile(outputFileName + "_local_minima_energies"  + ".txt");
	std::sort(localMinima.begin(), localMinima.end(),energyCompare);
	
	int lmi = 0;
	for (std::vector<structure*>::iterator it = localMinima.begin(); it != localMinima.end(); it++)
	{
		//e^ (gme - lme)*(hartree to joule conversion)/k_b T
		double partitionFunctionRatio = exp((globalMinimumEnergy - (*it)->energy) * (4.35974 / 1.380649) * (double)100000 / temperature);
		localMinimaEnergyFile << "structure: " << std::to_string(++lmi) << "\nenergy (Hartree): " << (*it)->energy << "\nratio at " + std::to_string(temperature) + " K: " << partitionFunctionRatio << std::endl;
	}

	localMinimaEnergyFile.close();
	writeToXyz(localMinima, outputFileName + "_" + "localMinima" + ".xyz");


	//this syntax should be removed just incase
	//delete everything except output string
	if (DEBUG) std::cout << "cleaning memory " << std::endl;
	for (int i = 0; i < structures.size(); i++)
	{
		delete structures[i];
	}
	//open a file and print the output
	






}

void basinHoppingCriteria(std::string previousRun, char seedMode, std::vector<int> composition, std::vector<std::string> elements, double x, double y, double z, double rlimit, double percentChangeI, double percentChangeF, double RIDpercent, int charge, int state, char calculator, std::string method, std::string basis, int optimizationIterations, int SCFsteps, int HFsteps, double time, int localMinimaTrapping, int radialCriteriaPercent, std::string outputFileName, std::string taskName, double deltaEThr = 1, bool checkEnergy = false, bool optimizePositive = true, bool keepFiles = false, bool filterOff = false, std::string converge = "1.0e-6", int timeOutMinutes = -1, int timeOutMinutesSP = -1, bool covalentCriteriaMode = false, double covalentCriteriaPercent = 0, int covalentCriteriaLimit = 1000, int bondRequirement = 0, bool planar = false)
{
	double RIDthreshold = 0;
	int energyCall = 0;
	std::vector<structure*> structures = {};
	//std::vector<std::set<int>> coordinationNumbers;
	structure* gmeS = nullptr;//the global minima structure known
	double globalMinimumEnergy = 0;//energy corresponding to gmeS
	double globalMinimumTime = 0;

	int totalStructuresSearchedAtGME = 0;

	structure* ls = nullptr; //current/last structure explored
	double lastEnergy = 0;//energy corresponding to ls

	structure* cMinima;//current local minima since the last seed reset
	double latestLowestEnergy;//energy corresponding to cminima
	std::vector<double> lowMovementVecor;//corresponding to cminima
	double cMinDeltaE;//delta E at cmin

	int stepsSinceNew = 0;
	
	std::vector<structure*> localMinima = {};
	int totalStructuresSearched = 0;

	//not used in continue
	structure* s = nullptr;
	double seedEnergy = 0;
	int ccSteps = 0;//current chain steps

	int calculations = 0;
	int step = 0;

	std::vector<std::pair<int, bool>> symmetries = possibleSymmetry(composition, elements, rlimit, radialCriteriaPercent, 1, 10000);
	int seedNo = 0;
	if (previousRun != "")
	{
		if (STDOUT) std::cout << "continuing from previous file_: " << previousRun << std::endl;
		energyCall = 1;
		std::ifstream startFile(previousRun);
		std::string inLine;
		std::vector<int> composition;
		std::vector<std::string> elements;
		int stepsSinceNew = 0;
		std::getline(startFile, inLine);// "Composition:\n";
		std::getline(startFile, inLine);
		std::pair< std::vector<int>, std::vector<std::string>> comp = getComposition(inLine);
		int debugo = 0;
		composition = comp.first;
		elements = comp.second;
		std::getline(startFile, inLine); //"threshold:\n";
		std::getline(startFile, inLine);
		//RIDthreshold = std::stoi(inLine);
		std::getline(startFile, inLine);// best structure
		std::getline(startFile, inLine);
		gmeS = new structure(inLine, elements, composition);
		globalMinimumEnergy = gmeS->energy;//a variable for monitoring the global minimum found
		globalMinimumTime = 0;
		std::getline(startFile, inLine);//# minutes searched prior
		int priorEnd = 0;
		for (int i = 0; i < inLine.size(); i++)
		{
			if (inLine[i] == ' ' && priorEnd == 0) priorEnd = i;
		}
		globalMinimumTime = std::stof(inLine.substr(0, priorEnd));
		std::getline(startFile, inLine);//# minutes searched prior
		priorEnd = 0;
		for (int i = 0; i < inLine.size(); i++)
		{
			if (inLine[i] == ' ' && priorEnd == 0) priorEnd = i;
		}
		totalStructuresSearchedAtGME = std::stof(inLine.substr(0, priorEnd));
		std::cout << "Current global minimum structure found after " << totalStructuresSearchedAtGME << " structures with energy " << globalMinimumEnergy << ": " << std::endl;
		gmeS->print();
		std::getline(startFile, inLine);//energy of best structure

		for (int i = 0; i < composition.size(); i++) for (int j = 0; j < composition[i]; j++) std::getline(startFile, inLine);//element followed by coordinates
		std::getline(startFile, inLine);// "last energy,structure:\n";
		std::getline(startFile, inLine);
		ls = new structure(inLine, elements, composition);
		lastEnergy = ls->energy;
		std::getline(startFile, inLine);//steps since new
		std::getline(startFile, inLine);
		stepsSinceNew = std::stoi(inLine);
		std::getline(startFile, inLine);//recycles since new
		std::getline(startFile, inLine);
		//recyclesSinceNew = std::stoi(inLine);
		std::getline(startFile, inLine);//uphill steps
		std::getline(startFile, inLine);
		//upHillSteps = std::stoi(inLine);
		std::getline(startFile, inLine);//uphill steps
		std::getline(startFile, inLine);
		//optimizationsInCurrentSeed = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		step = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		ccSteps = std::stoi(inLine);
		/*
		Todo set steps since new
		*/
		std::getline(startFile, inLine);// "local minima energy,structure:\n";
		localMinima = {};
		cMinima = nullptr;
		std::getline(startFile, inLine);
		while (inLine[0] != 'a')
		{
			structure* next = new structure(inLine, elements, composition);
			localMinima.push_back(next);
			cMinima = next;//constantly reset until we get the actual last one

			//load the next line. the reason this isn't done in the conditional is because we need to stop at a text line
			std::getline(startFile, inLine);
		}
		//std::getline(startFile, line);// "all structures energy,structure:\n";// this should already be ignored
		//get all the structures now
		latestLowestEnergy = cMinima->energy;
		while (std::getline(startFile, inLine))
		{
			structure* next = new structure(inLine, elements, composition);
			structures.push_back(next);
		}
		totalStructuresSearched = structures.size();
		// Close the file
		startFile.close();
	}
	else {
		RIDthreshold = threshold(RID, composition, elements, rlimit, rlimit, rlimit, 1, x, y, z, RIDpercent);
		stepsSinceNew = 0;
		structures = {};
		localMinima = {};
		totalStructuresSearched = structures.size();
		do {
			s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);
			while (!radialCriteria(*s, radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria (radial) in seed creation" << std::endl;
				delete s;
				s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);

			}
		} while (!bondingRequirement(s, radialCriteriaPercent, bondRequirement));
		energyCall = 0;
		double timeLeftMinutes = (double)time * 60;

		seedEnergy = energyCalculation(s, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
		gmeS = s;//the global minima structure known
		globalMinimumEnergy = seedEnergy;//energy corresponding to gmeS
		totalStructuresSearchedAtGME = 1;

		ls = s; //current/last structure explored
		lastEnergy = seedEnergy;//energy corresponding to ls

		cMinima = s;//current local minima since the last seed reset
		latestLowestEnergy = seedEnergy;//energy corresponding to cminima


		ccSteps = 0;//current chain steps


		step = 0;



	}


	cMinDeltaE = 0;//delta E at cmin

	/*
	Create seed
	*/



	//also initialize the energy changel

	double deltaE = -0.1 * deltaEThr;
	//change 10

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double minutes = (double)duration.count() / (double)60;
	if (DEBUG) std::cout << "minutes:  " << minutes << " ; time limit: " << double(time) * 60 - 30 << std::endl;
	int maxint_ = 0;//remove later!

	while (minutes < double(time) * 60 - 30)
	{
		ccSteps++;//current chain steps. gets reset every time a new seed is found, and is used to determine when to switch to blind gradient descent in the hybrid model. unfortunately it is distinct from all other step counters.
		maxint_++;
		if (DEBUG) std::cout << "MONITOR______________: " << " totalStructures: " << structures.size() << " total steps: " << maxint_ << std::endl;
		if (DEBUG || STDOUT) std::cout << "step: " << step << " time: " << minutes / 60 << std::endl;
		//step += 1;
		structure* Xnew = nullptr;

		if (STDOUT) std::cout << "Creating a new structure" << std::endl;
		//COME BACK HEREEEEEEE
		std::vector<double> movement = {};
		for (int i = 0; i < ls->elements.size(); i++)
		{
			for (int j = 0; j < ls->set[i].size(); j++)
			{
				movement.push_back(x * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
				movement.push_back(y * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
				movement.push_back(z * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
			}
		}
		std::pair<structure*, std::vector<double>> XnewPair = createXnew(ls, percentChangeI, percentChangeF, step, x, y, z, movement, false);
		if (covalentCriteriaMode)
		{
			int attempts = 0;
			while ((!covalentCriteria(XnewPair.first, covalentCriteriaPercent) || !bondingRequirement(XnewPair.first, covalentCriteriaPercent, bondRequirement)) && attempts < covalentCriteriaLimit)
			{
				if(XnewPair.first->composition.size() > 0) delete XnewPair.first;//this if statement only prevents some memory error i do not understand
				movement = {};
				for (int i = 0; i < ls->elements.size(); i++)
				{
					for (int j = 0; j < ls->set[i].size(); j++)
					{
						movement.push_back(x * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
						movement.push_back(y * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
						movement.push_back(z * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
					}
				}
				std::pair<structure*, std::vector<double>> XnewPair = createXnew(ls, percentChangeI, percentChangeF, step, x, y, z, movement, false);
				attempts++;
			}
			if (attempts >= covalentCriteriaLimit)
			{
				std::cout << "BGD algorithm trapped by covalent criteria. Returning to LLM with minor inefficency" << std::endl;
				//group theory gradient is stuck, return to the last local minima
				XnewPair = std::pair<structure*, std::vector<double>>(cMinima, lowMovementVecor);
				XnewPair.first->energy = latestLowestEnergy;

			}
		}
		Xnew = XnewPair.first;

		structures.push_back(Xnew);
		if (STDOUT)
		{
			std::cout << "New structure before optimization:" << std::endl;
			Xnew->print();
		}


		bool optimize = true;
		if (checkEnergy)
		{
			//prevent optimization on terrible structures
			double timeLeftMinutes = (double)time * 60 - minutes;

			double tempEnergy = energyCalculation(Xnew, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
			if (!optimizePositive) {
				if (tempEnergy < 0)
				{
					optimize = true;
				}
				else {
					optimize = false;
				}
			}
			Xnew->energy = tempEnergy;
		}
		std::cout << "to here? 3" << std::endl;
		if (optimizationIterations > 0 && optimize) {
			if (STDOUT) std::cout << "optimizing current new structure" << std::endl;

			std::string cmethod = method;
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
			minutes = (double)duration.count();
			double timeLeftMinutes = (double)time * 60 - minutes;

			structure* tempXnew = optimizationCalculation(calculator, Xnew, ++calculations, taskName, optimizationIterations, SCFsteps, x, y, z, movement, structures, basis, cmethod, charge, state, keepFiles, converge, timeLeftMinutes, timeOutMinutes);
			if (tempXnew == nullptr)
			{
				if (!checkEnergy) Xnew->energy = DBL_MAX;//as of now energy is not taken from a failed optimization, the failed optimization structure is also not used but in future versions this will be implemented. so the energy should be set to double max for the failure, therefore with check energy turned on the algorithm will be more proficient
			}
			else {
				Xnew = tempXnew;
				structures.push_back(Xnew);

			}

		}

		std::cout << "to here? 4" << std::endl;
		if (DEBUG) std::cout << "calling energy of Xnew " << std::endl;
		double newEnergy = Xnew->energy;//energy is now taken from the optimizers
			//energyCalculation(Xnew, basis, method, charge, state, calculator, taskName, energyCall, keepFiles);
		if (DEBUG) std::cout << "Xnew energy: " << newEnergy << std::endl;
		if (STDOUT) std::cout << "New structure energy: " << newEnergy << std::endl;
		if (STDOUT) Xnew->printZmatrix();

		totalStructuresSearched += 1;//only add 1 as we only checked the optimized energy, not the xnew energy in that case
		if (globalMinimumEnergy > newEnergy)
		{
			globalMinimumEnergy = newEnergy;
			globalMinimumTime = (double)duration.count();
			gmeS = Xnew;
			totalStructuresSearchedAtGME = totalStructuresSearched;
		}

		Xnew->energy = newEnergy;

		//Checking criteria
		if (newEnergy < lastEnergy || double(rand())/double(RAND_MAX) < exp((lastEnergy - newEnergy)/deltaEThr))
		{
			if (newEnergy < lastEnergy) {
				stepsSinceNew = 0;
				if (newEnergy < latestLowestEnergy)
				{
					latestLowestEnergy = newEnergy;

					cMinima = Xnew;
					lowMovementVecor = movement;
					cMinDeltaE = newEnergy - lastEnergy;
				}
			}
			//only accept during these criteria
			ls = Xnew;
			deltaE = newEnergy - lastEnergy;
			lastEnergy = newEnergy;
			
		}
		else
		{
			stepsSinceNew++;
		}
		
		step++;




		if (STDOUT) {
			std::cout << "CYCLE END: Current global minima energy, time,  and structure: " << globalMinimumEnergy << " " << globalMinimumTime << std::endl;
			gmeS->print();
			gmeS->printZmatrix();
		}



		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutes = (double)duration.count();
	}

	localMinima.push_back(cMinima);



	if (DEBUG || STDOUT) std::cout << "printing results " << std::endl;
	gmeS->writeToXyz(outputFileName + ".xyz");
	//gmeS->writeToObj(outputFileName + ".obj");
	/*
	print all structures
	*/
	std::string outputString = "";
	//add composition
	outputString += "Composition:\n";
	for (int i = 0; i < composition.size(); i++)
	{
		outputString += elements[i] + std::to_string(composition[i]);
	}
	outputString += "\n";

	//rid threshold
	outputString += "RID threshold:\n";
	outputString += std::to_string(0) + "\n";
	//add the Global minima structure
	outputString += "Best structure:\n";
	std::string line = "energy: " + std::to_string((gmeS)->energy);
	for (int e = 0; e < (gmeS)->set.size(); e++)
	{
		for (int a = 0; a < (gmeS)->set[e].size(); a++)
		{
			//each atom
			line += " " + (gmeS)->elements[e] + " " + std::to_string((gmeS)->set[e][a].x) + " " + std::to_string((gmeS)->set[e][a].y) + " " + std::to_string((gmeS)->set[e][a].z);
		}
	}
	line += "\n";
	outputString += line;

	std::string gline = std::to_string((int)globalMinimumTime) + " minutes prior\n" + std::to_string(totalStructuresSearchedAtGME) + " structures searched prior" + "\n" + "Energy: " + std::to_string(globalMinimumEnergy) + " \n";
	int zitt = 0;
	gmeS->makeZmatrix();
	for (int e = 0; e < gmeS->set.size(); e++)
	{
		for (int a = 0; a < gmeS->set[e].size(); a++)
		{
			//each atom
			gline += " " + gmeS->elements[e] + " " + std::to_string(gmeS->zmatrix[zitt][0]) + " " + std::to_string(gmeS->zmatrix[zitt][1]) + " " + std::to_string(gmeS->zmatrix[zitt][2]) + "\n";
			zitt++;
		}
	}
	/*
	for (int e = 0; e < gmeS->set.size(); e++)
	{
		for (int a = 0; a < gmeS->set[e].size(); a++)
		{
			//each atom
			gline += " " + gmeS->elements[e] + " " + std::to_string(gmeS->set[e][a].x) + " " + std::to_string(gmeS->set[e][a].y) + " " + std::to_string(gmeS->set[e][a].z) + "\n";
		}
	}*/

	outputString += gline;
	//add the last accepted structure
	outputString += "last energy,structure:\n";
	line = std::to_string(s->energy);
	for (int e = 0; e < s->set.size(); e++)
	{
		for (int a = 0; a < s->set[e].size(); a++)
		{
			//each atom
			line += " " + s->elements[e] + " " + std::to_string(s->set[e][a].x) + " " + std::to_string(s->set[e][a].y) + " " + std::to_string(s->set[e][a].z);
		}
	}
	line += "\n";
	outputString += line;
	//add the most recent steps since new
	outputString += "steps since new:\n";
	outputString += std::to_string(stepsSinceNew) + "\n";
	outputString += "recycles since new:\n";
	outputString += std::to_string(0) + "\n";
	outputString += "uphill steps since new:\n";
	outputString += std::to_string(0) + "\n";
	outputString += "optimizations in current seed:\n";
	outputString += std::to_string(0) + "\n";
	outputString += "step:\n";
	outputString += std::to_string(step) + "\n";
	outputString += "ccStep:\n";
	outputString += std::to_string(ccSteps) + "\n";
	//structures with energies
	outputString += "local minima only: energy,structure:\n";

	//int iii = 0;
	std::cout << "local minima size: " << localMinima.size() << std::endl;
	for (std::vector<structure*>::iterator c_structure = localMinima.begin(); c_structure != localMinima.end(); c_structure++)
	{

		//std::cout << "going through the localMinima " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}

	outputString += "all structures: energy,structure:\n";
	//iii = 0;
	for (std::vector<structure*>::iterator c_structure = structures.begin(); c_structure != structures.end(); c_structure++)
	{
		//for each structure
		//std::cout << "going through the all structures " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	int numStructures = structures.size();
	std::cout << "output filename: " << outputFileName + "_" + std::to_string(numStructures) + ".txt" << std::endl;
	std::ofstream outputfile(outputFileName + "_" + std::to_string(numStructures) + ".txt");

	// Write to the file
	outputfile << outputString;

	// Close the file
	outputfile.close();

	//this syntax should be removed just incase
	//delete everything except output string
	if (DEBUG) std::cout << "cleaning memory " << std::endl;
	for (int i = 0; i < structures.size(); i++)
	{
		delete structures[i];
	}
	//open a file and print the output







}

void basinHopping(std::string previousRun, char seedMode, std::vector<int> composition, std::vector<std::string> elements, double x, double y, double z, double rlimit, double percentChangeI, double percentChangeF, double RIDpercent, int charge, int state, char calculator, std::string method, std::string basis, int optimizationIterations, int SCFsteps, int HFsteps, double time, int localMinimaTrapping, std::string outputFileName, std::string taskName, double deltaEThr = 1, bool checkEnergy = false, bool optimizePositive = true, bool keepFiles = false, std::string converge = "1.0e-6", int timeOutMinutes = -1, int timeOutMinutesSP = -1, bool planar = false)
{
	double RIDthreshold = 0;
	int energyCall = 0;
	std::vector<structure*> structures = {};
	//std::vector<std::set<int>> coordinationNumbers;
	structure* gmeS = nullptr;//the global minima structure known
	double globalMinimumEnergy = 0;//energy corresponding to gmeS
	double globalMinimumTime = 0;

	int totalStructuresSearchedAtGME = 0;

	structure* ls = nullptr; //current/last structure explored
	double lastEnergy = 0;//energy corresponding to ls

	structure* cMinima;//current local minima since the last seed reset
	double latestLowestEnergy;//energy corresponding to cminima
	std::vector<double> lowMovementVecor;//corresponding to cminima
	double cMinDeltaE;//delta E at cmin

	int stepsSinceNew = 0;

	std::vector<structure*> localMinima = {};
	int totalStructuresSearched = 0;

	//not used in continue
	structure* s = nullptr;
	double seedEnergy = 0;
	int ccSteps = 0;//current chain steps

	int calculations = 0;
	int step = 0;

	std::vector<std::pair<int, bool>> symmetries = possibleSymmetry(composition, elements, rlimit,20, 1, 10000);
	int seedNo = 0;
	if (previousRun != "")
	{
		if (STDOUT) std::cout << "continuing from previous file_: " << previousRun << std::endl;
		energyCall = 1;
		std::ifstream startFile(previousRun);
		std::string inLine;
		std::vector<int> composition;
		std::vector<std::string> elements;
		int stepsSinceNew = 0;
		std::getline(startFile, inLine);// "Composition:\n";
		std::getline(startFile, inLine);
		std::pair< std::vector<int>, std::vector<std::string>> comp = getComposition(inLine);
		int debugo = 0;
		composition = comp.first;
		elements = comp.second;
		std::getline(startFile, inLine); //"threshold:\n";
		std::getline(startFile, inLine);
		//RIDthreshold = std::stoi(inLine);
		std::getline(startFile, inLine);// best structure
		std::getline(startFile, inLine);
		gmeS = new structure(inLine, elements, composition);
		globalMinimumEnergy = gmeS->energy;//a variable for monitoring the global minimum found
		globalMinimumTime = 0;
		std::getline(startFile, inLine);//# minutes searched prior
		int priorEnd = 0;
		for (int i = 0; i < inLine.size(); i++)
		{
			if (inLine[i] == ' ' && priorEnd == 0) priorEnd = i;
		}
		globalMinimumTime = std::stof(inLine.substr(0, priorEnd));
		std::getline(startFile, inLine);//# minutes searched prior
		priorEnd = 0;
		for (int i = 0; i < inLine.size(); i++)
		{
			if (inLine[i] == ' ' && priorEnd == 0) priorEnd = i;
		}
		totalStructuresSearchedAtGME = std::stof(inLine.substr(0, priorEnd));
		std::cout << "Current global minimum structure found after " << totalStructuresSearchedAtGME << " structures with energy " << globalMinimumEnergy << ": " << std::endl;
		gmeS->print();
		std::getline(startFile, inLine);//energy of best structure

		for (int i = 0; i < composition.size(); i++) for (int j = 0; j < composition[i]; j++) std::getline(startFile, inLine);//element followed by coordinates
		std::getline(startFile, inLine);// "last energy,structure:\n";
		std::getline(startFile, inLine);
		ls = new structure(inLine, elements, composition);
		lastEnergy = ls->energy;
		std::getline(startFile, inLine);//steps since new
		std::getline(startFile, inLine);
		stepsSinceNew = std::stoi(inLine);
		std::getline(startFile, inLine);//recycles since new
		std::getline(startFile, inLine);
		//recyclesSinceNew = std::stoi(inLine);
		std::getline(startFile, inLine);//uphill steps
		std::getline(startFile, inLine);
		//upHillSteps = std::stoi(inLine);
		std::getline(startFile, inLine);//uphill steps
		std::getline(startFile, inLine);
		//optimizationsInCurrentSeed = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		step = std::stoi(inLine);
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);
		ccSteps = std::stoi(inLine);
		/*
		Todo set steps since new
		*/
		std::getline(startFile, inLine);// "local minima energy,structure:\n";
		localMinima = {};
		cMinima = nullptr;
		std::getline(startFile, inLine);
		while (inLine[0] != 'a')
		{
			structure* next = new structure(inLine, elements, composition);
			localMinima.push_back(next);
			cMinima = next;//constantly reset until we get the actual last one

			//load the next line. the reason this isn't done in the conditional is because we need to stop at a text line
			std::getline(startFile, inLine);
		}
		//std::getline(startFile, line);// "all structures energy,structure:\n";// this should already be ignored
		//get all the structures now
		latestLowestEnergy = cMinima->energy;
		while (std::getline(startFile, inLine))
		{
			structure* next = new structure(inLine, elements, composition);
			structures.push_back(next);
		}
		totalStructuresSearched = structures.size();
		// Close the file
		startFile.close();
	}
	else {
		RIDthreshold = threshold(RID, composition, elements, rlimit, rlimit, rlimit, 1, x, y, z, RIDpercent);
		stepsSinceNew = 0;
		structures = {};
		localMinima = {};
		totalStructuresSearched = structures.size();
		s = getSeedFromSymmetryMode(seedMode, composition, elements, rlimit, seedNo, symmetries);
			
		energyCall = 0;
		double timeLeftMinutes = (double)time * 60;

		seedEnergy = energyCalculation(s, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
		gmeS = s;//the global minima structure known
		globalMinimumEnergy = seedEnergy;//energy corresponding to gmeS
		totalStructuresSearchedAtGME = 1;

		ls = s; //current/last structure explored
		lastEnergy = seedEnergy;//energy corresponding to ls

		cMinima = s;//current local minima since the last seed reset
		latestLowestEnergy = seedEnergy;//energy corresponding to cminima


		ccSteps = 0;//current chain steps


		step = 0;



	}


	cMinDeltaE = 0;//delta E at cmin

	/*
	Create seed
	*/



	//also initialize the energy changel

	double deltaE = -0.1 * deltaEThr;
	//change 10

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double minutes = (double)duration.count() / (double)60;
	if (DEBUG) std::cout << "minutes:  " << minutes << " ; time limit: " << double(time) * 60 - 30 << std::endl;
	int maxint_ = 0;//remove later!

	while (minutes < double(time) * 60 - 30)
	{
		ccSteps++;//current chain steps. gets reset every time a new seed is found, and is used to determine when to switch to blind gradient descent in the hybrid model. unfortunately it is distinct from all other step counters.
		maxint_++;
		if (DEBUG) std::cout << "MONITOR______________: " << " totalStructures: " << structures.size() << " total steps: " << maxint_ << std::endl;
		if (DEBUG || STDOUT) std::cout << "step: " << step << " time: " << minutes / 60 << std::endl;
		//step += 1;
		structure* Xnew = nullptr;

		if (STDOUT) std::cout << "Creating a new structure" << std::endl;
		//COME BACK HEREEEEEEE
		std::vector<double> movement = {};
		for (int i = 0; i < ls->elements.size(); i++)
		{
			for (int j = 0; j < ls->set[i].size(); j++)
			{
				movement.push_back(x * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
				movement.push_back(y * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
				movement.push_back(z * (2.0 * (rand() % RAND_MAX) / RAND_MAX - 1));
			}
		}
		std::pair<structure*, std::vector<double>> XnewPair = createXnew(ls, percentChangeI, percentChangeF, step, x, y, z, movement, false);
		Xnew = XnewPair.first;

		structures.push_back(Xnew);
		if (STDOUT)
		{
			std::cout << "New structure before optimization:" << std::endl;
			Xnew->print();
		}


		bool optimize = true;
		if (checkEnergy)
		{
			//prevent optimization on terrible structures
			double timeLeftMinutes = (double)time * 60 - minutes;

			double tempEnergy = energyCalculation(Xnew, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP);
			if (!optimizePositive) {
				if (tempEnergy < 0)
				{
					optimize = true;
				}
				else {
					optimize = false;
				}
			}
			Xnew->energy = tempEnergy;
		}
		if (optimizationIterations > 0 && optimize) {
			if (STDOUT) std::cout << "optimizing current new structure" << std::endl;

			std::string cmethod = method;
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
			minutes = (double)duration.count();
			double timeLeftMinutes = (double)time * 60 - minutes;

			structure* tempXnew = optimizationCalculation(calculator, Xnew, ++calculations, taskName, optimizationIterations, SCFsteps, x, y, z, movement, structures, basis, cmethod, charge, state, keepFiles, converge, timeLeftMinutes, timeOutMinutes);
			if (tempXnew == nullptr)
			{
				if (!checkEnergy) Xnew->energy = DBL_MAX;//as of now energy is not taken from a failed optimization, the failed optimization structure is also not used but in future versions this will be implemented. so the energy should be set to double max for the failure, therefore with check energy turned on the algorithm will be more proficient
			}
			else {
				Xnew = tempXnew;
				structures.push_back(Xnew);

			}

		}


		if (DEBUG) std::cout << "calling energy of Xnew " << std::endl;
		double newEnergy = Xnew->energy;//energy is now taken from the optimizers
			//energyCalculation(Xnew, basis, method, charge, state, calculator, taskName, energyCall, keepFiles);
		if (DEBUG) std::cout << "Xnew energy: " << newEnergy << std::endl;
		if (STDOUT) std::cout << "New structure energy: " << newEnergy << std::endl;
		if (STDOUT) Xnew->printZmatrix();

		totalStructuresSearched += 1;//only add 1 as we only checked the optimized energy, not the xnew energy in that case
		if (globalMinimumEnergy > newEnergy)
		{
			globalMinimumEnergy = newEnergy;
			globalMinimumTime = (double)duration.count();
			gmeS = Xnew;
			totalStructuresSearchedAtGME = totalStructuresSearched;
		}

		Xnew->energy = newEnergy;

		//Checking criteria
		if (newEnergy < lastEnergy || double(rand()) / double(RAND_MAX) < exp((lastEnergy - newEnergy) / deltaEThr))
		{
			if (newEnergy < lastEnergy) {
				stepsSinceNew = 0;
				if (newEnergy < latestLowestEnergy)
				{
					latestLowestEnergy = newEnergy;

					cMinima = Xnew;
					lowMovementVecor = movement;
					cMinDeltaE = newEnergy - lastEnergy;
				}
			}
			//only accept during these criteria
			ls = Xnew;
			deltaE = newEnergy - lastEnergy;
			lastEnergy = newEnergy;

		}
		else
		{
			stepsSinceNew++;
		}

		step++;




		if (STDOUT) {
			std::cout << "CYCLE END: Current global minima energy, time,  and structure: " << globalMinimumEnergy << " " << globalMinimumTime << std::endl;
			gmeS->print();
			gmeS->printZmatrix();
		}



		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutes = (double)duration.count();
	}

	localMinima.push_back(cMinima);



	if (DEBUG || STDOUT) std::cout << "printing results " << std::endl;
	gmeS->writeToXyz(outputFileName + ".xyz");
	//gmeS->writeToObj(outputFileName + ".obj");
	/*
	print all structures
	*/
	std::string outputString = "";
	//add composition
	outputString += "Composition:\n";
	for (int i = 0; i < composition.size(); i++)
	{
		outputString += elements[i] + std::to_string(composition[i]);
	}
	outputString += "\n";

	//rid threshold
	outputString += "RID threshold:\n";
	outputString += std::to_string(0) + "\n";
	//add the Global minima structure
	outputString += "Best structure:\n";
	std::string line = "energy: " + std::to_string((gmeS)->energy);
	for (int e = 0; e < (gmeS)->set.size(); e++)
	{
		for (int a = 0; a < (gmeS)->set[e].size(); a++)
		{
			//each atom
			line += " " + (gmeS)->elements[e] + " " + std::to_string((gmeS)->set[e][a].x) + " " + std::to_string((gmeS)->set[e][a].y) + " " + std::to_string((gmeS)->set[e][a].z);
		}
	}
	line += "\n";
	outputString += line;

	std::string gline = std::to_string((int)globalMinimumTime) + " minutes prior\n" + std::to_string(totalStructuresSearchedAtGME) + " structures searched prior" + "\n" + "Energy: " + std::to_string(globalMinimumEnergy) + " \n";
	int zitt = 0;
	gmeS->makeZmatrix();
	for (int e = 0; e < gmeS->set.size(); e++)
	{
		for (int a = 0; a < gmeS->set[e].size(); a++)
		{
			//each atom
			gline += " " + gmeS->elements[e] + " " + std::to_string(gmeS->zmatrix[zitt][0]) + " " + std::to_string(gmeS->zmatrix[zitt][1]) + " " + std::to_string(gmeS->zmatrix[zitt][2]) + "\n";
			zitt++;
		}
	}
	/*
	for (int e = 0; e < gmeS->set.size(); e++)
	{
		for (int a = 0; a < gmeS->set[e].size(); a++)
		{
			//each atom
			gline += " " + gmeS->elements[e] + " " + std::to_string(gmeS->set[e][a].x) + " " + std::to_string(gmeS->set[e][a].y) + " " + std::to_string(gmeS->set[e][a].z) + "\n";
		}
	}*/

	outputString += gline;
	//add the last accepted structure
	outputString += "last energy,structure:\n";
	line = std::to_string(s->energy);
	for (int e = 0; e < s->set.size(); e++)
	{
		for (int a = 0; a < s->set[e].size(); a++)
		{
			//each atom
			line += " " + s->elements[e] + " " + std::to_string(s->set[e][a].x) + " " + std::to_string(s->set[e][a].y) + " " + std::to_string(s->set[e][a].z);
		}
	}
	line += "\n";
	outputString += line;
	//add the most recent steps since new
	outputString += "steps since new:\n";
	outputString += std::to_string(stepsSinceNew) + "\n";
	outputString += "recycles since new:\n";
	outputString += std::to_string(0) + "\n";
	outputString += "uphill steps since new:\n";
	outputString += std::to_string(0) + "\n";
	outputString += "optimizations in current seed:\n";
	outputString += std::to_string(0) + "\n";
	outputString += "step:\n";
	outputString += std::to_string(step) + "\n";
	outputString += "ccStep:\n";
	outputString += std::to_string(ccSteps) + "\n";
	//structures with energies
	outputString += "local minima only: energy,structure:\n";

	//int iii = 0;
	std::cout << "local minima size: " << localMinima.size() << std::endl;
	for (std::vector<structure*>::iterator c_structure = localMinima.begin(); c_structure != localMinima.end(); c_structure++)
	{

		//std::cout << "going through the localMinima " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}

	outputString += "all structures: energy,structure:\n";
	//iii = 0;
	for (std::vector<structure*>::iterator c_structure = structures.begin(); c_structure != structures.end(); c_structure++)
	{
		//for each structure
		//std::cout << "going through the all structures " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	int numStructures = structures.size();
	std::cout << "output filename: " << outputFileName + "_" + std::to_string(numStructures) + ".txt" << std::endl;
	std::ofstream outputfile(outputFileName + "_" + std::to_string(numStructures) + ".txt");

	// Write to the file
	outputfile << outputString;

	// Close the file
	outputfile.close();

	//this syntax should be removed just incase
	//delete everything except output string
	if (DEBUG) std::cout << "cleaning memory " << std::endl;
	for (int i = 0; i < structures.size(); i++)
	{
		delete structures[i];
	}
	//open a file and print the output







}

//seed structure

int main(int argv, char *argc[])
{
	
	srand(time(0));
	//std::vector<structure*> sfxyz = structuresFromXYZ("testSeeds.xyz");
	//for (std::vector<structure*>::iterator it = sfxyz.begin(); it != sfxyz.end(); it++) std::cout << radialCriteria(**it, 20);
	//return 0;
	/*
	std::vector<atom> set;
	set.push_back(atom(0.690242, 3.18115, -0.17103));
	set.push_back(atom(-2.91749, -1.44376, -0.17103));
	set.push_back(atom(-2.00214, -2.56664, -0.17103));
	set.push_back(atom(-0.690242, -3.18115, -0.17103));
	set.push_back(atom(0.758364, -3.16561, -0.17103));
	set.push_back(atom(2.05677, -2.52307, -0.17103));
	set.push_back(atom(2.9478, -1.38081, -0.17103));

	structure pen = structure({ set }, { "B" });
	radialCriteria(pen, 20);
	*/
	//std::vector<atom> set = partialSymmetrySet(8, 5, 0, 5, 1);
	
	/*
	std::vector<bool*> ce = alignedEnumeration(7, 3);
	for (int i = 0; i < ce.size(); i++) {
		for (int j = 0; j < 7; j++) std::cout << ce[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::vector<bool*> ce2 = enumeration(7, 3);
	for (int i = 0; i < ce2.size(); i++) {
		for (int j = 0; j < 7; j++) std::cout << ce2[i][j] << " ";
		std::cout << std::endl;
	}*/

	///*
	//CHECK SEED GENERATION
	std::vector<structure*> sps = getSeeds({ 7 }, { "B" }, 2.5, 20, 1, 10000,false, false, true, false, 1);
	//std::vector<structure*> sps = seedsWithPartialSymmetry(false,true,false, 1.3, { 1,9 }, {"Ag","B"}, 4.5, 20, 1, 10000);
	writeToXyz(sps, "testSeeds.xyz");
	std::cout << sps.size() << std::endl;
	return 0;

	//*/




	/*
	std::vector<std::pair<int, bool>> as = possibleSymmetry({ 7 }, { "B" }, 3.5, 20, 2, 100,10000);
	for (std::vector<std::pair<int, bool>>::iterator it = as.begin(); it != as.end(); it++) std::cout << it->first << "," << it->second << std::endl;
	return 0;
	*/

	/*
	std::vector<std::string> elementsEE = { "La","B" };
	std::vector<int> compositionCC = { 2,8 };

	
	srand(time(0));
	for (int i = 0; i < 1; i++)
	{

		double rad = 3.0;
		int rcp = 20;
		int bondNum = 2;

		std::cout << "attempt " << i << std::endl;
		structure* seed;
		do {
			seed = seedWithSymmetry(8,1, compositionCC, elementsEE, rad);

			while (!radialCriteria(*seed, rcp))
			{
				delete seed;
				seed = seedWithSymmetry(8, 1, compositionCC, elementsEE, rad);
			}
		} while (!bondingRequirement(seed, rcp, bondNum));
		std::cout << "accepted" << std::endl;
		std::cout << bondingRequirement(seed, rcp, bondNum) << std::endl;
		seed->writeToXyz("checkLa2B8_" + std::to_string(i) + ".xyz");
	}
	return 9;


	std::vector<std::string> elementsE = { "B" };
	std::vector<int> compositionC = { 7 };

	double rad = 3.5;
	int rcp = 20;
	int bondNum = 2;
	srand(time(0));
	for (int i = 0; i < 10; i++)
	{
		std::cout << "attempt " << i << std::endl;
		structure* seed;
		do {
			seed = proceduralSeedFromGroupTheory(i,compositionC, elementsE, rad);

			while (!radialCriteria(*seed, rcp))
			{
				delete seed;
				seed = proceduralSeedFromGroupTheory(i, compositionC, elementsE, rad);
			}
			std::cout << "symmetry: " << seed->nsym << " " << seed->hsym << std::endl;
		}while(!bondingRequirement(seed, rcp, bondNum));
		std::cout << "accepted" << std::endl;
		std::cout << bondingRequirement(seed, rcp, bondNum) << std::endl;
		seed->writeToXyz("b7_" + std::to_string(i) + ".xyz");
	}
	return 0; 

	std::cout << "generating 10 seeds for 5 radial requirements and checking bond requirements" << std::endl;
	std::vector<double> radii = { 3.5, 3.0, 2.7, 2.5, 2.0,1.5 };

	for (int r = 0; r < radii.size(); r++)
	{
		std::cout << "radius: " << radii[r] << std::endl;
		structure* seed;
		for (int i = 0; i < 1; i++)
		{
			seed = seedFromGroupTheory(compositionC, elementsE, radii[r]);
			do {
				while (!radialCriteria(*seed, rcp))
				{
					delete seed;
					seed = seedFromGroupTheory(compositionC, elementsE, radii[r]);
				}
			} while (!bondingRequirement(seed, rcp, bondNum));
			std::cout << "seed passed" << std::endl;
		}
	}

	return 0;
	//*/
	
	/*
	-b basis
	-m  method
	-l previous file name
	-i structure initialization file
	-f output file
	-n taskName
	-comp composition
	-pr radial criteria percent
	-prti radial criteria trapping iteration - at least 1. default 5
	-prtf radial criteria trapping freedom - greater than 0. default 10
	-pi percent change i
	-pf percent change f
	-ps RID percent
	-charge charge
	-s state
	-o optimization cycles
	-t time
	-stepsC coordination steps
	-stepsL local minima steps
	-stepsS scf steps
	-xyz
	-ABC
	-g DirectionExponent proportionExponent Enormalization
	-h hybrid switching step
	-keepFiles
	-


	alternate modes:

	-vflag input output :  creates .obj file for  visualization of the structure
	-vgflag input output : creates .obj file for visualization of gaussian optimization output
	-cgo input output -b basis -m method -o optimizationcycles -stepsS SCFsteps -n taskName -charge charge -s state : takes a structure and creates a gaussian optimization input file for it
	*/
	
	std::string basis;
	std::string method;
	std::string previousFilename = "";
	std::string outputFileName;
	std::string initializationFile;
	bool initialized = false;
	std::string taskName;
	std::vector<std::string> elements;
	std::vector<int> composition;
	std::string compositionString;
	//std::string potential;
	double radialCriteriaPercent;
	double radialCriteriaTrappingIteration = 5;
	double radialCriteriaTrappingFreedom = 10;
	double percentChangeI;
	double percentChangeF;
	double percentRID = 2;
	bool keepFiles = false;
	int timeOutMinutes = -1;
	int timeOutMinutesSP = -1;

	double deltaEThr = 0;
	bool blindGradientDescent = false;
	int proportionExponent = 0;
	int directionExponent = 0;
	int hybridSteps = 0;


	int charge;
	int state;
	int optimizationCycles;
	double time;
	int coordinationSteps;
	int localMinimaTrappingSteps;
	int misstepTrapping;
	int HFsteps;
	int SCFsteps;
	std::string converge = "1.0e-6";
	char dft;//dft calculator for energy and optimization

	bool continuePastRun = false;
	double x, y, z;
	double a = 0; 
	double b = 0;
	double c = 0;
	double rlimit = 0;
	int uphillstopping = 1;


	bool fflag = false;
	bool bflag = false;
	bool pflag = false;
	bool cflag = false;
	int xyzflag = 0;
	bool mflag = false;
	bool tflag = false;
	bool sflag = false;

	bool filterOff = false;

	std::string structureName = "";
	std::string objName = "";

	//alternate modes
	bool vflag = false;
	bool vgflag = false;
	bool vlflag = false;
	double lmEthr = 1;
	bool cgo = false;
	bool groupTheoryFlag = false;

	bool refine = false;
	bool gtgd = false;
	bool checkEnergy = false;
	bool optimizePositive = false;

	bool planarLimitation = false;
	bool covalentCriteria = false;
	bool covalentCriteriaPercent = 0;
	int covalentCriteriaSteps = 0;
	int bondRequirement = 0;
	int criteriaIterations = 10000;

	bool BHC = false;//run basin hopping algorrithm with criteira
	bool BH = false;//run basin hopping algorithm
	bool BHS = false;

	char symmetryMode = 'n';

	bool partialSymmetry = false;
	bool extraHoles = false;
	bool alignHoles = true;
	bool variantSymmetry = false;
	double nMultiplier = 1;
	double temp = 293.15;
	/*
	/*
	//p is procedural, each symmetry allowed will be tried one by one
	//r is random, from allowed symmetries a random one will be tried
	//s this generates a seed with random symmetry, symmetries with more reasonable structures will be visited more
	//n no symmetry, a random seed
	*/

	//ModifiedBasinHopping.exe -c CH4 -xyz 5 5 5 -abc 20 20 20 -t 1 -s 100 -l 30 -n isItWorking testOutOne -r 50
	for (int i = 1; i < argv; i++)
	{
		if (i + 1 < argv) std::cout << "current argument: " << argc[i] << " next argument: " << argc[i + 1] << std::endl;
		else std::cout << "current argument: " << argc[i] << std::endl;
		if (i + 1 > argv)
		{
			std::cout << "flag out of bounds" << std::endl;
		}
		else {
		if ((std::string)argc[i] == "-l") {
			previousFilename = argc[i + 1];
			if (DEBUG) std::cout << "previous filename: " << previousFilename << std::endl;
			continuePastRun = true;
			i++;
		}else if ((std::string)argc[i] == "-refine") {
			previousFilename = argc[i + 1];
			if (DEBUG) std::cout << "previous filename: " << previousFilename << std::endl;
			refine = true;
			i++;
		}
		else if ((std::string)argc[i] == "-gt") {
			gtgd = true;
			proportionExponent = std::atoi(argc[i + 2]);
			directionExponent = std::atoi(argc[i + 1]);
			deltaEThr = std::stof(argc[i + 3]);
			i += 3;
			if (DEBUG) std::cout << "group theory gradient descent " << std::endl;
		}
		else if ((std::string)argc[i] == "-BH") {
			BH = true;
			deltaEThr = std::stof(argc[i + 1]);
			i ++;
			if (DEBUG) std::cout << "basin hopping" << std::endl;
		}
		else if ((std::string)argc[i] == "-BHC") {
			BHC = true;
			deltaEThr = std::stof(argc[i + 1]);
			i++;
			if (DEBUG) std::cout << "basin hopping with criteria" << std::endl;
		}
		else if ((std::string)argc[i] == "-BHS") {
			BHS = true;
			if (DEBUG) std::cout << "basin hopping with symmetry" << std::endl;
		}
		else if ((std::string)argc[i] == "-r") {
			rlimit = std::atof(argc[i + 1]);
			if (DEBUG) std::cout << "radial limi: " << rlimit << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-stepsU") {
			uphillstopping = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "uphillstopping: " <<uphillstopping << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-d")
		{
			std::string dftS = argc[i + 1];
			dft = dftS[0];
			if (DEBUG) std::cout << "dftS: " << dft << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-b")
		{
			basis = argc[i + 1];
			if (DEBUG) std::cout << "basis: " << basis << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-m")
		{
			method = argc[i + 1];
			if (DEBUG) std::cout << "method: " << method << std::endl;
			i++;
		}
		else if (std::string(argc[i]) == "-comp")
		{
			compositionString = std::string( argc[i + 1]);
			if (DEBUG) std::cout << "composition string: " << compositionString << std::endl;
			std::pair< std::vector<int>, std::vector<std::string>> comppair = getComposition(compositionString);
			composition = comppair.first;
			elements = comppair.second;
			
			for (std::vector<int>::iterator it = composition.begin(); it != composition.end(); it++)
			{
				if (DEBUG) std::cout << *it << " ";
			}
			for (std::vector<std::string>::iterator it = elements.begin(); it != elements.end(); it++)
			{
				if (DEBUG) std::cout << *it << " ";
			}
			if (DEBUG) std::cout << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-xyz" || (std::string)argc[i] == "-XYZ")
		{
			x = std::atof(argc[i + 1]);
			y = std::atof(argc[i + 2]);
			z = std::atof(argc[i + 3]);
			if (DEBUG) std::cout << "xyz: " << x << " " << y << " " << z << std::endl;
			i += 3;
		}
		else if ((std::string)argc[i] == "-abc" || (std::string)argc[i] == "-ABC")
		{
			a = std::atof(argc[i + 1]);
			b = std::atof(argc[i + 2]);
			c = std::atof(argc[i + 3]);
			if (DEBUG) std::cout << "abc: " << a << " " << b << " " << c << std::endl;
			i += 3;
		}
		else if ((std::string)argc[i] == "-t")
		{
			time = std::atof(argc[i + 1]);
			if (DEBUG) std::cout << "time allowed: " << time << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-stepsC")
		{
			coordinationSteps = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "coordination steps: " << coordinationSteps << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-stepsL")
		{
			localMinimaTrappingSteps = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "local minima trapping steps: " << localMinimaTrappingSteps << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-stepsM")
		{
			misstepTrapping = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "misstep trapping steps: " << misstepTrapping << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-stepsH")
		{
			HFsteps = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "Hartree fock steps: " << HFsteps << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-stepsS")
		{
			SCFsteps = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "scf steps: " << SCFsteps << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-converge")
		{
		converge = argc[i + 1];
		if (DEBUG) std::cout << "scf converge: " << converge << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-n")
		{
			taskName = argc[i + 1];
			if (DEBUG) std::cout << "task name: " << taskName << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-f")
		{
			outputFileName = argc[i + 1];
			if (DEBUG) std::cout << "output filename: " << outputFileName << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-i")
		{
		initializationFile = argc[i + 1];
		initialized = true;
		if (DEBUG) std::cout << "initializationfilename: " << initializationFile << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-pr")
		{
			radialCriteriaPercent = std::stof(argc[i + 1]);
			if (DEBUG) std::cout << "radial criteria percent: " << radialCriteriaPercent << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-cc")
		{
		covalentCriteria = true;
		covalentCriteriaPercent = std::stof(argc[i + 1]);
		covalentCriteriaSteps = std::stof(argc[i + 2]);
		if (DEBUG) std::cout << "covalent criteria percent: " << covalentCriteriaPercent << " , covalent criteria steps: " << covalentCriteriaSteps << std::endl;
		i+=2;
		}
		else if ((std::string)argc[i] == "-br")
		{
		bondRequirement = std::stof(argc[i + 1]);
		if (DEBUG) std::cout << "bond requirement: " << bondRequirement << std::endl;
		i ++;
		}
		else if ((std::string)argc[i] == "-cI")
		{
		criteriaIterations = std::stoi(argc[i + 1]);
		if (DEBUG) std::cout << "criteria iterations: " << criteriaIterations << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-prti")
		{
		radialCriteriaTrappingIteration = std::stof(argc[i + 1]);
		if (DEBUG) std::cout << "radial criteria trapping iteration: " << radialCriteriaTrappingIteration << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-prtf")
		{
		radialCriteriaTrappingFreedom = std::stof(argc[i + 1]);
		if (DEBUG) std::cout << "radial criteria trapping freedom: " << radialCriteriaTrappingFreedom << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-pi")
		{
			percentChangeI = std::stof(argc[i + 1]);
			if (DEBUG) std::cout << "percent change I: " << percentChangeI << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-pf")
		{
			percentChangeF = std::stof(argc[i + 1]);
			if (DEBUG) std::cout << "percent change F: " << percentChangeF << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-ps")
		{
		percentRID = std::stof(argc[i + 1]);
		if (DEBUG) std::cout << "RID percent: " << percentRID << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-charge")
		{
			charge = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "charge: " << charge << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-s")
		{
			state = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "state: " << state << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-o")
		{
		optimizationCycles = std::atoi(argc[i + 1]);
		if (DEBUG) std::cout << "optimization cycles: " << optimizationCycles << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-g")
		{
		blindGradientDescent = true;
		proportionExponent = std::atoi(argc[i + 2]);
		directionExponent = std::atoi(argc[i + 1]);
		deltaEThr = std::stof(argc[i + 3]);
		if (DEBUG) std::cout << "blind gradient descent: " << proportionExponent << " " << directionExponent << " " << deltaEThr << std::endl;
		i+=3;
		}
		else if ((std::string)argc[i] == "-h")
		{
		hybridSteps = std::atoi(argc[i + 1]);
		if (DEBUG) std::cout << "hybrid steps: " << hybridSteps << std::endl;
		i++;
		}
		else if ((std::string)argc[i] == "-to")
		{
		timeOutMinutes = std::atoi(argc[i + 1]);
		timeOutMinutesSP = std::atoi(argc[i + 2]);
		if (DEBUG) std::cout << "timeoutminutes: " << timeOutMinutes << std::endl;
		i+=2;
		}
		else if ((std::string)argc[i] == "-keepFiles" || (std::string)argc[i] == "-keep" || (std::string)argc[i] == "-keepfiles")
		{
			std::cout << "keep files on" << std::endl;
			keepFiles = true;
		}
		else if ((std::string)argc[i] == "-filterOff" || (std::string)argc[i] == "-filteroff" || (std::string)argc[i] == "-nofilter")
		{
		std::cout << "similarity filter turned off" << std::endl;
		filterOff = true;
		}
		else if ((std::string)argc[i] == "-v")
		{
		structureName = argc[i + 1];
		objName = argc[i + 2];
		vflag = true;
		i+=2;
		}
		else if ((std::string)argc[i] == "-vg")
		{
		structureName = argc[i + 1];
		objName = argc[i + 2];
		vgflag = true;
		i += 2;
		}
		else if ((std::string)argc[i] == "-vl")
		{
		structureName = argc[i + 1];
		objName = argc[i + 2];
		lmEthr = std::stod(argc[i + 3]);
		vlflag = true;
		i += 3;
		}
		else if ((std::string)argc[i] == "-groupTheorySeeds")
		{
		//use composition string as well
		objName = argc[i + 1];
		groupTheoryFlag = true;
		i += 1;
		}
		else if ((std::string)argc[i] == "-seed")
		{
		//use composition string as well
		symmetryMode = argc[i + 1][0];

		/*
		case 'i'://i for iterate through the possible symmetries
		case 'r'://r for randomly choose from the possible symmetries
		case 'g'://symetrical seed with random symmetry. g for group theory seed
		case 'n'://no symmetry
		case 'o'://either symmetrical only, or also partial symmetry. create all seeds first, then order by energy and begin blind gradient descent.

		bool partialSymmetry = false;
		bool extraHoles = false;
		bool alignHoles = true;
		bool variantSymmetry = false;
		double nMultiplier = 1;
		double temp = 293.15;
		*/
			if (symmetryMode == 'o')
			{
				partialSymmetry = std::stoi(argc[i + 2]);
				extraHoles = std::stoi(argc[i + 3]);
				alignHoles = std::stoi(argc[i + 4]);
				variantSymmetry = std::stoi(argc[i + 5]);
				nMultiplier = std::stod(argc[i + 6]);
				i += 6;
			}
			else {
				i++;
			}

		
		}
		else if ((std::string)argc[i] == "-temp")
		{
		temp = std::stod(argc[i + 1]);
		i ++;
		}
		else if ((std::string)argc[i] == "-cgo")
		{
		structureName = argc[i + 1];
		cgo = true;
		objName = argc[i + 2];
		i += 2;
		}
		
		else if ((std::string)argc[i] == "-check" || (std::string)argc[i] == "-checkEnergy")
		{
			checkEnergy = true;
		}
		else if ((std::string)argc[i] == "-op" || (std::string)argc[i] == "-optimizePositive")
		{
			optimizePositive = true;
			checkEnergy = true;
		}
		else if ((std::string)argc[i] == "-planar")
		{
			planarLimitation = true;
		}
		else {
			std::cout << argc[i] << " not identified as a flag" << std::endl;
		}
		}



		
	}
	if (DEBUG) std::cout << "opening if else " << std::endl;
	if (vflag)
	{

		std::vector<int> compositionV;
		std::vector<std::string> elementsV;
	
		

		std::ifstream startFile(structureName);
		std::string inLine;
		std::getline(startFile, inLine);// "Composition:\n";
		std::getline(startFile, inLine);//composition
		std::pair< std::vector<int>, std::vector<std::string>> compV = getComposition(inLine);
		compositionV = compV.first;
		elementsV = compV.second;
		std::getline(startFile, inLine);//;"Cell size:\n";
		std::getline(startFile, inLine); //std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
		std::getline(startFile, inLine); //"threshold:\n";
		std::getline(startFile, inLine);//threshold
		std::getline(startFile, inLine);// best structure
		std::getline(startFile, inLine);
		std::cout << inLine << std::endl;
		structure readStructure(inLine, elementsV, compositionV);
		readStructure.writeToObj(objName);

	}
	else if (vgflag)
	{
		//same as vflag except its a gaussian outwput file instead
		std::vector<int> compositionV;
		std::vector<std::string> elementsV;

		std::ifstream gOutFile(structureName);

		std::string cline;
		structure* st = nullptr;

		std::string matchString = "Input orientation:";
		while (std::getline(gOutFile, cline)) {
			if (cline.length() > matchString.length())
			{

				bool matching = false;
				int matchIt = 0;
				for (int l = 0; l < cline.length() - matchString.length(); l++)
				{
					if (cline[l] == matchString[matchIt])
					{
						matchIt++;
						matching = true;
					}
					else if (matchIt < 18) {
						matchIt = 0;
						matching = false;
					}
					else {
					}
				}
				if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
				if (matching && matchIt >= matchString.length())
				{
					std::vector<std::vector<atom>> set;
					std::vector<std::string> elements;
					getline(gOutFile, cline);
					getline(gOutFile, cline);
					getline(gOutFile, cline);
					getline(gOutFile, cline);
					char space = ' ';


					while (std::getline(gOutFile, cline) && cline[1] != '-')
					{
						std::string cnumber = "";

						int num = 0;
						std::string element;
						double x = 0;
						double y = 0;
						double z = 0;
						for (int i = 0; i < cline.length(); i++)
						{
							if (cline[i] == space)
							{
								if (cnumber.length() > 0)
								{
									num++;
									switch (num)
									{
									case 1:
										break;
									case 2:
										element = atomicNumber(std::stoi(cnumber));
										break;
									case 3:
										break;
									case 4:
										x = std::stof(cnumber);
										break;
									case 5:
										y = std::stof(cnumber);
										break;
									case 6:
										z = std::stof(cnumber);
										break;
									default:
										break;
									}
									cnumber = "";
								}
							}
							else {
								cnumber += cline[i];
							}
						}
						int index = -1;
						for (int i = 0; i < elements.size(); i++)
						{
							if (elements[i] == element)
							{
								index = i;
							}
						}
						if (index == -1)
						{
							elements.push_back(element);
							std::vector < atom> alist = { atom(x,y,z) };
							set.push_back(alist);
						}
						else {
							set[index].push_back(atom(x, y, z));
						}
					}
					if (st != nullptr) delete st;
					st = new structure(set, elements);
				}
			}
		}
		gOutFile.close();

		st->writeToObj(objName);

		delete st;
	}
	else if (vlflag) {
		//obj name should not end in .obj
		std::vector<int> compositionV;
		std::vector<std::string> elementsV;



		std::ifstream startFile(structureName);
		std::string inLine;
		std::getline(startFile, inLine);// "Composition:\n";
		std::getline(startFile, inLine);//composition
		std::pair< std::vector<int>, std::vector<std::string>> compV = getComposition(inLine);
		compositionV = compV.first;
		elementsV = compV.second;
		std::getline(startFile, inLine);//;"rlimit:\n";
		std::getline(startFile, inLine); //std::to_string(r) 
		std::getline(startFile, inLine); //"threshold:\n";
		std::getline(startFile, inLine);//threshold
		std::getline(startFile, inLine);// best structure
		std::getline(startFile, inLine);
		structure bestStructure(inLine, elementsV, compositionV);
		double bestEnergy = bestStructure.energy;
		//readStructure.writeToObj(objName + "_Best.obj" );
		std::getline(startFile, inLine);//structures searched prior
		std::getline(startFile, inLine);//energy
		for(int i = 0; i < compositionV.size(); i++) for(int e = 0; e < compositionV[i];e++) std::getline(startFile, inLine);//zmatrix lines
		std::getline(startFile, inLine);//last energy structure
		std::getline(startFile, inLine);//energy
		std::getline(startFile, inLine);//steps since new
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);//recycles since new
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);//uphill
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);//optimizations in current seed
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);//step
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);//ccStep
		std::getline(startFile, inLine);
		std::getline(startFile, inLine);//local minima only
		std::vector<structure*> desiredMinima = {};//minima within the threshold from best
		std::getline(startFile, inLine);
		do {
			//std::cout << inLine << std::endl;
			structure* cStructure = new structure(inLine, elementsV, compositionV);
			//std::cout << "succeded to make structure" << std::endl;
			if (cStructure->energy < bestEnergy + lmEthr) desiredMinima.push_back(cStructure);
			else delete cStructure;
			//std::cout << "failed to get enrgy" << std::endl;
			std::getline(startFile, inLine);
		} while (inLine[0] != 'a');
		std::cout << "debug 3" << std::endl;
		for (int s = 0; s < desiredMinima.size(); s++)
		{
			desiredMinima[s]->writeToObj(objName + "_" + std::to_string(s) + ".obj");
			delete desiredMinima[s];
		}
		std::cout << "debug 4" << std::endl;

	}
	else if (cgo)
	{
		std::vector<int> compositionV;
		std::vector<std::string> elementsV;

		std::ifstream startFile(structureName);
		std::string inLine;
		std::getline(startFile, inLine);// "Composition:\n";
		std::getline(startFile, inLine);//composition
		std::pair< std::vector<int>, std::vector<std::string>> compV = getComposition(inLine);
		compositionV = compV.first;
		elementsV = compV.second;
		std::getline(startFile, inLine);//;"Cell size:\n";
		std::getline(startFile, inLine); //std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
		std::getline(startFile, inLine); //"threshold:\n";
		std::getline(startFile, inLine);//threshold
		std::getline(startFile, inLine);// best structure
		std::getline(startFile, inLine);
		std::cout << inLine << std::endl;
		structure readStructure(inLine, elementsV, compositionV);
		
		writeToGuassian(readStructure, taskName, objName, basis, method, charge, state, optimizationCycles, SCFsteps);
		

	}
	else if (groupTheoryFlag)
	{
		//required other flags:
		//composition
		//radial criteria percent
		//ModifiedBasinHopping.exe -groupTheorySeeds b8seeds/b8 -comp B8 -pr 60 -r 2
		std::vector<structure*> structures;
		for (int i = 0; i < 100; i++)
		{
			structure* s = seedFromGroupTheory(composition, elements, rlimit);
			while (!radialCriteria(*s, radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria (radial) in seed creation" << std::endl;
				delete s;
				s = seedFromGroupTheory(composition, elements, rlimit);
			}
			std::string objFileName = objName + "_C" + std::to_string(s->nsym);
			if (s->hsym != 0) objFileName += "h";
			objFileName += "_" + std::to_string(i) + ".obj";

			s->writeToObj(objFileName);
			delete s;
		}
		
		
	}
	else {
		if (gtgd)
		{
			groupTheoryGradientDescent(previousFilename,composition, elements, x, y, z, rlimit,symmetryMode, percentChangeI, percentChangeF, percentRID, charge, state, dft, method, basis, optimizationCycles, SCFsteps, HFsteps, time, localMinimaTrappingSteps, uphillstopping, misstepTrapping, radialCriteriaPercent,  outputFileName, taskName, directionExponent, proportionExponent, deltaEThr,checkEnergy,optimizePositive, keepFiles,filterOff,converge,timeOutMinutes,timeOutMinutesSP,covalentCriteria,covalentCriteriaPercent,covalentCriteriaSteps,bondRequirement,criteriaIterations,planarLimitation ,partialSymmetry,extraHoles,alignHoles,variantSymmetry,nMultiplier , temp);
			/*
			ModifiedBasinHopping.exe -gt 1 1 0.1 -m PBEPBE -b '6-31+g(d,p)' -d b -f latest -n gtb -comp B6 -pr 60 -pi 10 -pf 100 -charge -1 -s 1 -o 30 -t 0.505 -stepsL 30 -stepsM 5 -stepsU 10 -stepsH 10 -stepsS 129 -xyz 0.2 0.2 0.2 -r 2 -ps 2
			*/
		}else if (BH)
		{
			basinHopping(previousFilename,symmetryMode,composition,elements,x,y,z,rlimit,percentChangeI,percentChangeF,percentRID,charge,state,dft,method,basis,optimizationCycles,SCFsteps,HFsteps,time,localMinimaTrappingSteps,outputFileName,taskName,deltaEThr,checkEnergy,optimizePositive,keepFiles,converge,timeOutMinutes,timeOutMinutesSP);
		}
		else if (BHC)
		{
			basinHoppingCriteria(previousFilename,symmetryMode,composition,elements,x,y,z,rlimit,percentChangeI,percentChangeF,percentRID,charge,state,dft,method,basis,optimizationCycles,SCFsteps,HFsteps,time,localMinimaTrappingSteps,radialCriteriaPercent,outputFileName,taskName,deltaEThr,checkEnergy,optimizePositive,keepFiles, filterOff, converge,timeOutMinutes,timeOutMinutesSP,covalentCriteria,covalentCriteriaPercent,covalentCriteriaSteps,bondRequirement);
		}
	}
	
}

