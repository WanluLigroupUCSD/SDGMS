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
#include <fstream>
#include <stdlib.h>

const int DEBUG = 1;

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
	double halfValue = std::log(0.5);
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
	}
	*/
	if (DEBUG > 1) std::cout << "returning" << std::endl;
	return 0.5 * epsilon * energy;
};

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
				if (DEBUG) std::cout << "did not pass criteria (radial) " << std::endl;
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
			if (DEBUG) std::cout << "structure accepted" << std::endl;
		}

		


		
	}
	return retValue;
}

void startModifiedBasinHopping(double x, double y, double z,double percentChangeI, double percentChangeF, std::vector<int> composition, std::vector<std::string> elements,std::function<double(structure&, int, std::string, double, double, double)> energy, std::function<structure*(structure&)> optimize,std::function<bool(structure&,int)> radialCriteria,std::function<bool(structure&, std::vector<std::set<int>>)> coordinationNumber,int coordinationSteps, int time, int localMinimaTrapping,int radialCriteriaPercent, std::string outputFileName, std::string taskName,double a, double b, double c)
{
	std::vector<structure*> structures = { };
	std::vector<structure*> localMinima = {};
	std::vector<std::set<int>> coordinationNumbers;
	double RIDthreshold = threshold(RID,composition, x, y, z, 1, 0);
	structure* s = getAcceptedStructure(x, y, z, elements, composition, structures,localMinima, coordinationSteps, coordinationNumbers, radialCriteriaPercent, RIDthreshold);
	

		//new structure(x, y, z, composition, elements);
	int energyCall = 1;
	double seedEnergy = energy(*s,energyCall++,taskName,a,b,c);
	if(DEBUG) std::cout << "seed structure energy: " << seedEnergy << std::endl;
	//check it passes radial criteria?
	s->energy = seedEnergy;

	structures.push_back(s);
	structure* lastLocalMinima = nullptr;
	
	if (DEBUG) std::cout << "RID threshold: " << RIDthreshold << std::endl;

	int step = 0;
	int stepsSinceNew = 0;

	

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double hours = (double)duration.count() / (double)60;
	if (DEBUG) std::cout << "chours:  " << hours << " ; time limit: " << double(time) - 0.5 <<  std::endl;
	//int maxint_ = 0;//remove later!
	while (hours < double(time) - 0.5/* && 50 > maxint_++*/)
	{
		if (DEBUG) std::cout << "step: " << step << " time: " << double(time) << std::endl;
		step += 1;
		structure* Xnew = nullptr;

		//to avoid local minima trapping, if for a while a new structure cant be accepted then move to a new one
		if (stepsSinceNew > localMinimaTrapping)
		{
			if (DEBUG) std::cout << "steps since new passed "  << std::endl;
			//should we restart the step size here for the purpose of deviation from seed?
			step = 1;
			//to avoid local minima trapping find a new (unique) seed
			
			bool seedAccepted = false;
			if (DEBUG) std::cout << "looking for a new seed "  << std::endl;
			/*
			while (!seedAccepted)
			{
				//cant delete s because it might still be the last good structure
				//delete s;
				s = new structure(x, y, z, composition, elements);
				std::vector<structure*>::iterator it = structures.begin();

				if (localMinima.size() >= coordinationSteps)
				{
					while (!coordinationNumber(*s,coordinationNumbers) || !radialCriteria(*s,radialCriteriaPercent))
					{
						delete s;
						s = new structure(x, y, z, composition, elements);
					}
				}
				else {
					while (!radialCriteria(*s,radialCriteriaPercent))
					{
						delete s;
						s = new structure(x, y, z, composition, elements);
					}
				}
				bool unique = true;

				while (it != structures.end() && unique)
				{
					double dist = RID(*s, *(*it));
					if (dist < RIDthreshold)
					{
						unique = false;
					}
					it++;
				}
				seedAccepted = unique;
				if (!seedAccepted)
				{
					delete s;
					//delete s here because we know we arent accepting it and we cant delete it at the top of the loop
				}
			}
			*/
			s = getAcceptedStructure(x, y, z, elements, composition, structures, localMinima, coordinationSteps, coordinationNumbers, radialCriteriaPercent, RIDthreshold);
			// a new seed to start with in s that is unique and passes criteria
			double seedEnergy = energy(*s,energyCall++, taskName, a, b, c);
			if (DEBUG) std::cout << "new seed found: " << seedEnergy << std::endl;
			structures.push_back(s);
			stepsSinceNew = 0;
		}
		//std::cout << "remove later, doing a completely random structure" << std::endl;
		//energy(structure(*s, step), energyCall++, taskName, a, b, c);

		

		if (DEBUG) std::cout << "creating a new structure " << std::endl;

		//determine variation based on step
		double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));
		double variationX = x* sigmoidPercentageFromStep /100;
		double variationY = y * sigmoidPercentageFromStep / 100;
		double variationZ = z * sigmoidPercentageFromStep / 100;
		Xnew = new structure(*s, variationX, variationY, variationZ, a,b,c);
		//tests
		/*
		//std::cout << "calculating energy right now!" << std::endl;
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << " basic call worked" << std::endl;
		radialCriteria(*Xnew, radialCriteriaPercent);
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << " radial criteria worked" << std::endl;
		scoreAtoms(*Xnew);
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << " the sub task worked" << std::endl;
		RID(*Xnew, *(*structures.begin()));
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << "rid worked" << std::endl;
		*/
		//end of tests

		//as the commented out energy calls  worked above, something after this point causes the energy check to fail.
		
		if (localMinima.size() == coordinationSteps)
		{
			if (DEBUG) std::cout << "generating coordination numbers " << std::endl;
			coordinationNumbers = generateAcceptedCoordinationNumbers(localMinima);
		}
		if (localMinima.size() >= coordinationSteps)
		{
			while (!coordinationNumber(*Xnew,coordinationNumbers) || !radialCriteria(*Xnew,radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria ";
				if (DEBUG) {
					if (coordinationNumber(*Xnew, coordinationNumbers)) std::cout << " radial" << std::endl;
					else std::cout << " coordination numbers" << std::endl;
				}
				delete Xnew;
				double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));
				double variationX = x * sigmoidPercentageFromStep / 100;
				double variationY = y * sigmoidPercentageFromStep / 100;
				double variationZ = z * sigmoidPercentageFromStep / 100;
				Xnew = new structure(*s, variationX, variationY, variationZ, a, b, c);
			}
		}
		else {
			while (!radialCriteria(*Xnew,radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria (radial) " << std::endl;
				delete Xnew;
				double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));
				double variationX = x * sigmoidPercentageFromStep / 100;
				double variationY = y * sigmoidPercentageFromStep / 100;
				double variationZ = z * sigmoidPercentageFromStep / 100;
				Xnew = new structure(*s, variationX, variationY, variationZ, a, b, c);
			}
		}


		//check that the new structure does not violate radius criteria
		if (DEBUG) std::cout << "checking for uniqueness " << std::endl;

		//check for uniqueness
		std::vector<structure*>::iterator it = structures.begin();
		bool unique = true;
		while (it != structures.end() && unique)
		{
			double dist = RID(*Xnew, *(*it));
			if(dist < RIDthreshold)
			{
				unique = false;
			}
			it++;
		}


		if (unique)
		{
			if (DEBUG) std::cout << "unique " << std::endl;
			// add to structures
			structures.push_back(Xnew);
			/*
			when do you add a structure to local minima?
			TODO
			*/
			//optimize Xnew
			//Xnew = optimize(*Xnew);
			// this is causing a crash because there is no optimize
			//should we add the optimized structure to structures? check if it already is in structures?
			/*
			TODO
			*/
			if (DEBUG) std::cout << "calling energy of Xnew " << std::endl;
			double newEnergy = energy(*Xnew, energyCall++, taskName, a, b, c);
			if (DEBUG) std::cout << "Xnew energy: " << newEnergy << std::endl;
			Xnew->energy = newEnergy;
			bool accept = false;

			if (newEnergy < seedEnergy)
			{
				if (DEBUG) std::cout << "accepted for lower, set new local minima " << std::endl;
				accept = true;
				stepsSinceNew = 0;
				lastLocalMinima = Xnew;
			}
			else
			{
				/*
				current seed is a local minima. add to local minima?
				*/


				double p = rand() / RAND_MAX;

				double pThreshold = std::exp((seedEnergy - newEnergy)/(boltzman*temperature(stepsSinceNew,100,1.5)));//e^(Es-Ex)/KbT
				if (p < pThreshold)
				{
					if (DEBUG) std::cout << "accepted by statistics, pushed last local minima " << std::endl;
					accept = true;

					//somehow its possible to get to this line of code when lastLocalMinima has not yet been set.
					//this only occurs at the beginning, so maybe we should push back the seed in this case, however I am just going to not push anything in that case for now although likely it should be the seed.
					if (lastLocalMinima != nullptr) {
						std::cout << "pushing the last local minima" << lastLocalMinima->energy << std::endl;
						localMinima.push_back(lastLocalMinima);
						std::cout << " its energy is " << localMinima[localMinima.size() - 1]->energy << std::endl;
					}

						
				}
				
			}
			if (accept)
			{
				stepsSinceNew = 0;
				//localMinima.push_back(Xnew);
				//now only adding to local minima when moving to a higher energy or starting over
				s = Xnew;
				seedEnergy = newEnergy;
				step++;
				//need some way to save all energies?, maybe a variable in structure
			}else {
				if (DEBUG) std::cout << "rejected " << std::endl;
				stepsSinceNew++;
			}


		}
		else
		{
			if (DEBUG) std::cout << "not unique" << std::endl;
			stepsSinceNew++;
		}

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		hours = (double)duration.count() / (double)60;
	}
	if (DEBUG) std::cout << "prinitng results " << std::endl;
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

	//add cell size x y z
	outputString += "Cell size:\n";
	outputString += std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";

	//rid threshold
	outputString += "threshold:\n";
	outputString += std::to_string(RIDthreshold)  + "\n";
	//add the last accepted structure
	outputString += "last energy,structure:\n";
	std::string line = std::to_string(s->energy);
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
	//coordination numbers
	outputString += "coordination numbers:\n";
	if (coordinationSteps <= localMinima.size())
	{
		coordinationNumbers = generateAcceptedCoordinationNumbers(localMinima);
	}

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
	/*
	//the below code had a memory error which i dont understand
	for (int i = 0; i < localMinima.size(); i++)
	{
		//for each structure
		std::string line = "energy: " + std::to_string(localMinima[i]->energy);
		for (int e = 0; e < localMinima[i]->set.size(); e++)
		{
			for (int a = 0; a < localMinima[i]->set[e].size(); a++)
			{
				//each atom
				line += " " + localMinima[i]->elements[e] + " " + std::to_string(localMinima[i]->set[e][a].x) + " " + std::to_string(localMinima[i]->set[e][a].y) + " " + std::to_string(localMinima[i]->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	*/
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
	//this syntax should be removed just incase
	/*
	for (int i = 0; i < structures.size(); i++)
	{
		//for each structure
		std::string line = "energy: " + std::to_string(structures[i]->energy);
		for (int e = 0; e < structures[i]->set.size(); e++)
		{
			for (int a = 0; a < structures[i]->set[e].size(); a++)
			{
				//each atom
				line += " " + structures[i]->elements[e] + " " + std::to_string(structures[i]->set[e][a].x) + " " + std::to_string(structures[i]->set[e][a].y) + " " + std::to_string(structures[i]->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	*/
	int numStructures = structures.size();
	//delete everything except output string
	if (DEBUG) std::cout << "cleaning memory " << std::endl;
	for (int i = 0; i < structures.size(); i++)
	{
		delete structures[i];
	}
	//open a file and print the output
	std::ofstream outputfile(outputFileName + "_" + std::to_string(numStructures) + ".txt");
	
	// Write to the file
	outputfile << outputString;

	// Close the file
	outputfile.close();

}

std::pair< std::vector<int>, std::vector<std::string>> getComposition(std::string comp)
{
	std::vector<int> composition;
	std::vector<std::string> elements;
	bool instring = true;//assume we start in a string
	std::string current = "";
	int currint = 0;
	for (int i = 0; i < comp.size(); i++)
	{
		if (comp[i] > 48 && comp[i] < 57)
		{
			//numeric
			if (i > 0) {
				if (instring) {
					if (DEBUG) std::cout << "added " << current << " to elements" << std::endl;
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
					if (DEBUG) std::cout << "added " << currint << " to composition" << std::endl;
					//switch
					instring = true;
					composition.push_back(currint);
					currint = 0;
					
				}
				else {
					if (DEBUG) std::cout << "added " << current << " to elements and " << 1 << " to composition" << std::endl;
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
					if (DEBUG) std::cout << "added " << currint << " to composition" << std::endl;
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
		if (DEBUG) std::cout << "added " << currint << " to composition" << std::endl;
		composition.push_back(currint);
		currint = 0;
		
	}
	else {
		if (DEBUG) std::cout << "added " << current << " to elements and " << 1 << " to composition" << std::endl;
		//elements begin with capitals so the last one is over
		elements.push_back(current);
		composition.push_back(1);//omly one of the last element
		current = "";
		
	}
	return std::pair< std::vector<int>, std::vector<std::string>>(composition, elements);
}

void continueModifiedBasinHopping(std::string startFileName, std::function<double(structure&, int, std::string, double, double, double)>  energy, std::function<structure* (structure&)> optimize, std::function<bool(structure&,int)> radialCriteria, std::function<bool(structure&, std::vector<std::set<int>>)> coordinationNumber, double percentChangeI, double percentChangeF, int coordinationSteps, int time, int localMinimaTrapping, int radialCriteriaPercent, std::string outputFileName,std::string taskName, double a, double b, double c)
{
	int energyCall = 1;
	std::ifstream startFile(startFileName);
	std::string line;

	double RIDthreshold;
	std::vector<structure*> structures;
	double x;
	double y;
	double z;
	std::vector<int> composition;
	std::vector<std::string> elements;
	structure* s;
	int stepsSinceNew;

	std::getline(startFile, line);// "Composition:\n";
	std::getline(startFile, line);
	std::pair< std::vector<int>, std::vector<std::string>> comp = getComposition(line);
	int debugo = 0;
	composition = comp.first;
	elements = comp.second;
	
	std::getline(startFile, line);//;"Cell size:\n";
	std::getline(startFile, line); //std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
	int xstart = 0;
	int ystart = 0;
	int zstart = 0;
	for (int i = 0; i < line.size(); i++)
	{
		if (line[i] == ' ')
		{
			if (ystart == 0) ystart = i;
			else zstart = i;
		}
	}
	x = std::stof(line.substr(xstart, ystart-xstart));
	y = std::stof(line.substr(ystart, zstart - ystart));
	z = std::stof(line.substr(zstart, zstart-line.size()));

	std::getline(startFile, line); //"threshold:\n";
	std::getline(startFile, line);
	RIDthreshold = std::stoi(line);

	std::getline(startFile, line);// "last energy,structure:\n";
	std::getline(startFile, line);
	s = new structure(line, elements, composition);
	/*
	Todo get the last structure
	
	*/
	std::getline(startFile, line);// steps since new
	std::getline(startFile, line);
	stepsSinceNew = std::stoi(line);
	/*
	Todo set steps since new
	*/
	std::getline(startFile, line);// "coordination numbers:\n"; ?? what should i do with this?? - we recalculate later?
	std::getline(startFile, line);// "local minima energy,structure:\n";
	std::vector<structure*> previousLocalMinima;
	std::vector<structure*> allLocalMinima = {};
	structure* lastLocalMinima = nullptr;
	std::getline(startFile, line);
	while (line[0] != 'a')
	{
		structure* next = new structure(line, elements, composition);
		previousLocalMinima.push_back(next);
		allLocalMinima.push_back(next);
		lastLocalMinima = next;//constantly reset until we get the actual last one
		
		//load the next line. the reason this isn't done in the conditional is because we need to stop at a text line
		std::getline(startFile, line);
	}
	//std::getline(startFile, line);// "all structures energy,structure:\n";// this should already be ignored

	int totalStructures = 0;
	while (std::getline(startFile, line))
	{
		totalStructures++;
	}
	

	// Close the file
	startFile.close();



	
	double seedEnergy = s->energy;
	structures.push_back(s);
	//check it passes radial criteria?


	//std::vector<structure*> structures = { s };
	
	
	//double RIDthreshold = threshold(RID, s->composition, x, y, z, 1, 0);

	int step = 0;
	//int stepsSinceNew = 0;
	std::vector<std::set<int>> coordinationNumbers = generateAcceptedCoordinationNumbers(previousLocalMinima);

	auto start = std::chrono::high_resolution_clock::now();
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double hours = (double)duration.count() / (double)60;

	//int maxint_ = 0;//remove later!
	while (hours < double(time) - 0.5/* && 50 > maxint_++*/)
	{
		if (DEBUG) std::cout << "step: " << step << " time: " << double(time) << std::endl;
		step += 1;
		structure* Xnew = nullptr;

		//to avoid local minima trapping, if for a while a new structure cant be accepted then move to a new one
		if (stepsSinceNew > localMinimaTrapping)
		{
			if (DEBUG) std::cout << "steps since new passed "  << std::endl;
			//should we restart the step size here for the purpose of deviation from seed?
			step = 1;
			//to avoid local minima trapping find a new (unique) seed
			
			bool seedAccepted = false;
			if (DEBUG) std::cout << "looking for a new seed "  << std::endl;
			/*
			while (!seedAccepted)
			{
				//cant delete s because it might still be the last good structure
				//delete s;
				s = new structure(x, y, z, composition, elements);
				std::vector<structure*>::iterator it = structures.begin();

				if (localMinima.size() >= coordinationSteps)
				{
					while (!coordinationNumber(*s,coordinationNumbers) || !radialCriteria(*s,radialCriteriaPercent))
					{
						delete s;
						s = new structure(x, y, z, composition, elements);
					}
				}
				else {
					while (!radialCriteria(*s,radialCriteriaPercent))
					{
						delete s;
						s = new structure(x, y, z, composition, elements);
					}
				}
				bool unique = true;

				while (it != structures.end() && unique)
				{
					double dist = RID(*s, *(*it));
					if (dist < RIDthreshold)
					{
						unique = false;
					}
					it++;
				}
				seedAccepted = unique;
				if (!seedAccepted)
				{
					delete s;
					//delete s here because we know we arent accepting it and we cant delete it at the top of the loop
				}
			}
			*/
			s = getAcceptedStructure(x, y, z, elements, composition, structures, allLocalMinima, coordinationSteps, coordinationNumbers, radialCriteriaPercent, RIDthreshold);
			// a new seed to start with in s that is unique and passes criteria
			double seedEnergy = energy(*s,energyCall++, taskName, a, b, c);
			if (DEBUG) std::cout << "new seed found: " << seedEnergy << std::endl;
			structures.push_back(s);
			stepsSinceNew = 0;
		}
		//std::cout << "remove later, doing a completely random structure" << std::endl;
		//energy(structure(*s, step), energyCall++, taskName, a, b, c);

		if (DEBUG) std::cout << "creating a new structure " << std::endl;

		
		double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));
		double variationX = x * sigmoidPercentageFromStep / 100;
		double variationY = y * sigmoidPercentageFromStep / 100;
		double variationZ = z * sigmoidPercentageFromStep / 100;
		Xnew = new structure(*s, variationX, variationY, variationZ, a, b, c);
		//tests
		//
		/*
		//std::cout << "calculating energy right now!" << std::endl;
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << " basic call worked" << std::endl;
		radialCriteria(*Xnew, radialCriteriaPercent);
		std::cout << "testing" << std::endl;
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << " radial criteria worked" << std::endl;
		scoreAtoms(*Xnew);
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << " the sub task worked" << std::endl;
		RID(*Xnew, *(*structures.begin()));
		std::cout << " the sub task worked" << std::endl;
		energy(*Xnew, energyCall++, taskName, a, b, c);
		std::cout << "rid worked" << std::endl;
		*/
		//end of tests

		//as the commented out energy calls  worked above, something after this point causes the energy check to fail.
		
		if (allLocalMinima.size() == coordinationSteps)
		{
			if (DEBUG) std::cout << "generating coordination numbers " << std::endl;
			coordinationNumbers = generateAcceptedCoordinationNumbers(allLocalMinima);
		}
		if (allLocalMinima.size() >= coordinationSteps)
		{
			while (!coordinationNumber(*Xnew,coordinationNumbers) || !radialCriteria(*Xnew,radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria ";
				if (DEBUG) {
					if (coordinationNumber(*Xnew, coordinationNumbers)) std::cout << " radial" << std::endl;
					else std::cout << " coordination numbers" << std::endl;
				}
				delete Xnew;
				double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));
				double variationX = x * sigmoidPercentageFromStep / 100;
				double variationY = y * sigmoidPercentageFromStep / 100;
				double variationZ = z * sigmoidPercentageFromStep / 100;
				Xnew = new structure(*s, variationX, variationY, variationZ, a, b, c);
			}
		}
		else {
			while (!radialCriteria(*Xnew,radialCriteriaPercent))
			{
				if (DEBUG) std::cout << "did not pass criteria (radial) " << std::endl;
				delete Xnew;
				double sigmoidPercentageFromStep = percentChangeI + (percentChangeF - percentChangeI) * (1 / (1 + exp(-step + 4)));
				double variationX = x * sigmoidPercentageFromStep / 100;
				double variationY = y * sigmoidPercentageFromStep / 100;
				double variationZ = z * sigmoidPercentageFromStep / 100;
				Xnew = new structure(*s, variationX, variationY, variationZ, a, b, c);
			}
		}


		//check that the new structure does not violate radius criteria
		if (DEBUG) std::cout << "checking for uniqueness " << std::endl;

		//check for uniqueness
		std::vector<structure*>::iterator it = structures.begin();
		bool unique = true;
		while (it != structures.end() && unique)
		{
			double dist = RID(*Xnew, *(*it));
			if(dist < RIDthreshold)
			{
				unique = false;
			}
			it++;
		}


		if (unique)
		{
			if (DEBUG) std::cout << "unique " << std::endl;
			// add to structures
			structures.push_back(Xnew);
			/*
			when do you add a structure to local minima?
			TODO
			*/
			//optimize Xnew
			//Xnew = optimize(*Xnew);
			// this is causing a crash because there is no optimize
			//should we add the optimized structure to structures? check if it already is in structures?
			/*
			TODO
			*/
			if (DEBUG) std::cout << "calling energy of Xnew " << std::endl;
			double newEnergy = energy(*Xnew, energyCall++, taskName, a, b, c);
			if (DEBUG) std::cout << "Xnew energy: " << newEnergy << std::endl;
			Xnew->energy = newEnergy;
			bool accept = false;

			if (newEnergy < seedEnergy)
			{
				if (DEBUG) std::cout << "last energy: " << seedEnergy << " new energy: " << newEnergy << std::endl;
				if (DEBUG) std::cout << "accepted for lower, set new local minima " << std::endl;
				accept = true;
				stepsSinceNew = 0;
				lastLocalMinima = Xnew;
			}
			else
			{
				if (DEBUG) std::cout << "last energy: " << seedEnergy << " new energy: " << newEnergy << std::endl;
				/*
				current seed is a local minima. add to local minima?
				*/


				double p = rand() / RAND_MAX;

				double pThreshold = std::exp((seedEnergy - newEnergy)/(boltzman*temperature(stepsSinceNew,100,1.5)));//e^(Es-Ex)/KbT
				if (p < pThreshold)
				{
					if (DEBUG) std::cout << "accepted by statistics, pushed last local minima " << std::endl;
					accept = true;

					//somehow its possible to get to this line of code when lastLocalMinima has not yet been set.
					//this only occurs at the beginning, so maybe we should push back the seed in this case, however I am just going to not push anything in that case for now although likely it should be the seed.
					if (lastLocalMinima != nullptr) {
						std::cout << "pushing the last local minima" << lastLocalMinima->energy << std::endl;
						allLocalMinima.push_back(lastLocalMinima);
						std::cout << " its energy is " << allLocalMinima[allLocalMinima.size() - 1]->energy << std::endl;
					}

						
				}
				else {
					if (DEBUG) std::cout << "rejected by statistics" << std::endl;
				}
				
			}
			if (accept)
			{
				stepsSinceNew = 0;
				//localMinima.push_back(Xnew);
				//now only adding to local minima when moving to a higher energy or starting over
				s = Xnew;
				seedEnergy = newEnergy;
				step++;
				//need some way to save all energies?, maybe a variable in structure
			}else {
				if (DEBUG) std::cout << "rejected " << std::endl;
				stepsSinceNew++;
			}


		}
		else
		{
			if (DEBUG) std::cout << "not unique" << std::endl;
			stepsSinceNew++;
		}

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		hours = (double)duration.count() / (double)60;
	}
	if (DEBUG) std::cout << "prinitng results " << std::endl;
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

	//add cell size x y z
	outputString += "Cell size:\n";
	outputString += std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";

	//rid threshold
	outputString += "threshold:\n";
	outputString += std::to_string(RIDthreshold)  + "\n";
	//add the last accepted structure
	outputString += "last energy,structure:\n";
	std::string line2 = std::to_string(s->energy);
	for (int e = 0; e < s->set.size(); e++)
	{
		for (int a = 0; a < s->set[e].size(); a++)
		{
			//each atom
			line2 += " " + s->elements[e] + " " + std::to_string(s->set[e][a].x) + " " + std::to_string(s->set[e][a].y) + " " + std::to_string(s->set[e][a].z);
		}
	}
	line2 += "\n";
	outputString += line2;
	//add the most recent steps since new
	outputString += "steps since new:\n";
	outputString += std::to_string(stepsSinceNew) + "\n";
	//coordination numbers
	outputString += "coordination numbers:\n";
	if (coordinationSteps <= allLocalMinima.size())
	{
		coordinationNumbers = generateAcceptedCoordinationNumbers(allLocalMinima);
	}

	//structures with energies
	outputString += "local minima only: energy,structure:\n";

	//int iii = 0;
	std::cout << "local minima size: " << allLocalMinima.size() << std::endl;
	for (std::vector<structure*>::iterator c_structure = allLocalMinima.begin(); c_structure != allLocalMinima.end(); c_structure++)
	{
		
		//std::cout << "going through the localMinima " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line3 = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line3 += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line3 += "\n";
		outputString += line3;
	}
	/*
	//the below code had a memory error which i dont understand
	for (int i = 0; i < localMinima.size(); i++)
	{
		//for each structure
		std::string line = "energy: " + std::to_string(localMinima[i]->energy);
		for (int e = 0; e < localMinima[i]->set.size(); e++)
		{
			for (int a = 0; a < localMinima[i]->set[e].size(); a++)
			{
				//each atom
				line += " " + localMinima[i]->elements[e] + " " + std::to_string(localMinima[i]->set[e][a].x) + " " + std::to_string(localMinima[i]->set[e][a].y) + " " + std::to_string(localMinima[i]->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	*/
	outputString += "all structures: energy,structure:\n";
	//iii = 0;
	for (std::vector<structure*>::iterator c_structure = structures.begin(); c_structure != structures.end(); c_structure++)
	{
		//for each structure
		//std::cout << "going through the all structures " << iii++ << " e: " << (*c_structure)->energy << std::endl;//remove later
		std::string line3 = "energy: " + std::to_string((*c_structure)->energy);
		for (int e = 0; e < (*c_structure)->set.size(); e++)
		{
			for (int a = 0; a < (*c_structure)->set[e].size(); a++)
			{
				//each atom
				line3 += " " + (*c_structure)->elements[e] + " " + std::to_string((*c_structure)->set[e][a].x) + " " + std::to_string((*c_structure)->set[e][a].y) + " " + std::to_string((*c_structure)->set[e][a].z);
			}
		}
		line3 += "\n";
		outputString += line3;
	}
	//this syntax should be removed just incase
	/*
	for (int i = 0; i < structures.size(); i++)
	{
		//for each structure
		std::string line = "energy: " + std::to_string(structures[i]->energy);
		for (int e = 0; e < structures[i]->set.size(); e++)
		{
			for (int a = 0; a < structures[i]->set[e].size(); a++)
			{
				//each atom
				line += " " + structures[i]->elements[e] + " " + std::to_string(structures[i]->set[e][a].x) + " " + std::to_string(structures[i]->set[e][a].y) + " " + std::to_string(structures[i]->set[e][a].z);
			}
		}
		line += "\n";
		outputString += line;
	}
	*/
	int numStructures = structures.size();
	//delete everything except output string
	if (DEBUG) std::cout << "cleaning memory " << std::endl;
	for (int i = 0; i < structures.size(); i++)
	{
		delete structures[i];
	}
	//open a file and print the output
	//num structures is current and total structures is from previous runs. i believe there is an error with this
	std::ofstream outputfile(outputFileName + "_" + std::to_string(numStructures + totalStructures) + ".txt");
	
	// Write to the file
	outputfile << outputString;

	// Close the file
	outputfile.close();















}

structure* optimize(structure& s)
{
	/*
	Todo
	*/
	//thsi causes a crash dont use!
	return nullptr;
}

//seed structure

int main(int argv, char *argc[])
{
	///*
	std::vector<int> compositionC = { 1,4 };
	std::vector<std::string> elementsE = { "C","H" };
	//structure sss(8, 8, 8, compositionC, elementsE);
	//structure ss2(8, 8, 8, compositionC, elementsE);
	//RID(sss, ss2);
	//structure ss3(8, 8, 8, compositionC, elementsE);
	//RID(ss3, ss2);
	//std::cout << "starting tests" << std::endl;
	//scoreAtoms(sss);
	//basinHoppingEnergy(sss, 1000, "no", 20, 20, 20);
	//std::cout << "starting tests pass?" << std::endl;
	//return 0;
	//
	//startModifiedBasinHopping(5, 5, 5,10,20, compositionC, elementsE, basinHoppingEnergy, optimize, radialCriteria, coordinationNumber, 100000, 1, 100, 30, "newStart", "testingOut", 20, 20, 20);
	continueModifiedBasinHopping("newStart_51.txt", basinHoppingEnergy, optimize, radialCriteria, coordinationNumber, 10, 20, 100000, 1, 100, 30, "newStart_cont", "testingOut", 20, 20, 20);

	return 0;
	
	//*/
	/*
	necessary arguments:
	old or new project
	if an old project use -f filename
	if a new project dont

	basis set:
	-b GTH-Whatever

	psuedopotentials
	-p GTH-Whatever

	composition:
	-c B20

	cell size:
	-xyz x y z
	
	maximum cycles:
	-m 8000

	maximum time (hours):
	-t 24

	max coordination steps:
	-s steps


	*/
	std::string previousFilename;
	std::string basis;
	std::string potential;
	std::string compositionString;
	std::vector<int> composition;
	std::vector<std::string> elements;
	std::string taskName;
	std::string outputFileName;
	int radialCriteriaPercent;

	int cycles;
	int time;
	int coordinationSteps;
	int localMinimaTrappingSteps;

	bool continuePastRun = false;
	double x, y, z;
	double a = 0; 
	double b = 0;
	double c = 0;

	bool fflag = false;
	bool bflag = false;
	bool pflag = false;
	bool cflag = false;
	int xyzflag = 0;
	bool mflag = false;
	bool tflag = false;
	bool sflag = false;
	//ModifiedBasinHopping.exe -c CH4 -xyz 5 5 5 -abc 20 20 20 -t 1 -s 100 -l 30 -n isItWorking testOutOne -r 50
	for (int i = 1; i < argv; i++)
	{
		std::cout << "current argument: " << argc[i]  << " next argument: " << argc[i+1] << std::endl;
		if (i + 1 > argv)
		{
			std::cout << "flag out of bounds" << std::endl;
		}
		else if(argc[i+1][0] == '-') {
			std::cout << "missing argument for: " << argc[i] << std::endl;
		}
		else {
		if ((std::string)argc[i] == "-f") {
			previousFilename = argc[i + 1];
			if (DEBUG) std::cout << "previous filename: " << previousFilename << std::endl;
			continuePastRun = true;
			i++;
		}
		else if ((std::string)argc[i] == "-b")
		{
			basis = argc[i + 1];
			if (DEBUG) std::cout << "basis: " << basis << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-p")
		{
			potential = argc[i + 1];
			if (DEBUG) std::cout << "potential: " << potential << std::endl;
			i++;
		}
		else if (std::string(argc[i]) == "-c")
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
		else if ((std::string)argc[i] == "-xyz")
		{
			x = std::atof(argc[i + 1]);
			y = std::atof(argc[i + 2]);
			z = std::atof(argc[i + 3]);
			if (DEBUG) std::cout << "xyz: " << x << " " << y << " " << z << std::endl;
			i += 3;
		}
		else if ((std::string)argc[i] == "-abc")
		{
			a = std::atof(argc[i + 1]);
			b = std::atof(argc[i + 2]);
			c = std::atof(argc[i + 3]);
			if (DEBUG) std::cout << "abc: " << a << " " << b << " " << c << std::endl;
			i += 3;
		}
		else if ((std::string)argc[i] == "-m")
		{
			cycles = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "cycles?: " << cycles << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-t")
		{
			time = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "time allowed: " << time << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-s")
		{
			coordinationSteps = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "coordination steps: " << coordinationSteps << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-l")
		{
			localMinimaTrappingSteps = std::atoi(argc[i + 1]);
			if (DEBUG) std::cout << "local minima trapping steps: " << localMinimaTrappingSteps << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-n")
		{
			taskName = argc[i + 1];
			if (DEBUG) std::cout << "task name: " << taskName << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-o")
		{
			outputFileName = argc[i + 1];
			if (DEBUG) std::cout << "output filename: " << outputFileName << std::endl;
			i++;
		}
		else if ((std::string)argc[i] == "-r")
		{
			radialCriteriaPercent = std::stoi(argc[i + 1]);
			if (DEBUG) std::cout << "radial criteria percent: " << radialCriteriaPercent << std::endl;
			i++;
		}
		else {
			std::cout << argc[i] << " not identified as a flag" << std::endl;
		}
		}



		
	}

	if (DEBUG) std::cout << "opening if else " << std::endl;
	if (continuePastRun)
	{
		if (DEBUG) std::cout << "continuing modified basin hopping: " << std::endl;
		//continueModifiedBasinHopping(previousFilename, basinHoppingEnergy, optimize, radialCriteria, coordinationNumber, coordinationSteps, time, localMinimaTrappingSteps, radialCriteriaPercent,outputFileName, taskName, a, b, c);
	}
	else {
		if (DEBUG) std::cout << "start modified basin hopping: " << std::endl;
		//startModifiedBasinHopping(x, y, z, composition, elements, basinHoppingEnergy, optimize, radialCriteria, coordinationNumber, coordinationSteps, time, localMinimaTrappingSteps, radialCriteriaPercent,outputFileName, taskName, a, b, c);
	}


	//basic energy retrieval test
	/*
	srand(time(NULL));

	std::vector<int> composition = { 1,4 };
	std::vector<std::string> elements = { "C","H" };
	structure s(8, 8, 8, composition, elements);

	int num = rand() % 100;

	std::cout << energyCP2K(s, num, "methaneTest", 8, 8, 8, "");
	*/
	
	
}

