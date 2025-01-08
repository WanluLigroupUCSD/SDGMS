// ModifiedBasinHopping.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Written by Jordan Burkhardt for the Wan-Lu Li lab at UCSD

#include <iostream>
#include "chemistry.hpp"
#include "DistanceHeuristics.hpp"
#include "calculators.hpp"
#include <functional>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdlib.h>
#include <stack>
#include <queue>
#include <iomanip>

#ifndef DBL_MAX
#define DBL_MAX 1.79769e+308;
#endif

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
/*
void seedSearch(std::vector<int> composition, std::vector<std::string> elements, std::vector<int> fix, double rlimit, int charge, int state, char calculator, std::string method, std::string basis, int optimizationIterations, int SCFsteps, double time, int noptimizations, int radialCriteriaPercent, std::string outputFileName, std::string taskName,  std::vector<int> numRings, std::vector<int> nValues ,std::string converge = "1.0e-6", int timeOutMinutes = -1, int timeOutMinutesSP = -1, int bondRequirement = 0, int criteriaIterations = 10000, bool unrestricted = false, std::string seedFile = "empty", std::string seedEnergyFile = "empty",bool keepFiles = true, std::string customBasis = "empty")
{
	auto start = std::chrono::high_resolution_clock::now();
	int energyCall = 0;
	std::vector<structure*> structures = {};
	//std::vector<std::set<int>> coordinationNumbers;

	int totalStructuresSearched = 0;

	//not used in continue
	structure* s = nullptr;
	double seedEnergy = 0;

	int calculations = 0;

	std::vector<std::pair<int, bool>> symmetries = {};//used in seed mode i & r
	std::vector<structure*> seeds = {};//used in seed mode o
	
	int seedNo = 0;

	if (seedFile == "empty")
	{
		//make the seeds
		std::cout << "Creating Seeds, time(m): " << (double)time * 60 << std::endl;
		for (int i = 0; i < numRings.size(); i++) {
			std::vector<structure*> tempSeeds = proceduralSeeds(composition, elements, numRings[i], nValues, radialCriteriaPercent);
			for (std::vector<structure*>::iterator it = tempSeeds.begin(); it != tempSeeds.end(); it++) seeds.push_back(*it);
		}
		writeToXyz(seeds, outputFileName + "_unsorted_seeds" + ".xyz");
	}
	else {
		seeds = XyzToStructures(seedFile, elements);
	}

	if(seedEnergyFile == "empty") {
		//calculate seed energies
		std::cout << "Calculating energy of " << seeds.size() << " seeds, time(m) : " << (double)time * 60 << std::endl;
		for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
		{
			double timeLeftMinutes = (double)time * 60;
			double currentEnergy = energyCalculation(*it, basis, method, charge, state, calculator, taskName, ++calculations, keepFiles, converge, timeLeftMinutes, timeOutMinutesSP, unrestricted,customBasis);
			(*it)->energy = currentEnergy;
		}
		//sort seeds by energy
		std::cout << "Sorting seeds, time(m): " << (double)time * 60 << std::endl;
		std::sort(seeds.begin(), seeds.end(), energyCompare);
		std::ofstream outputfile(outputFileName + "_seed_energies" + ".txt");
		// Write to the file
		for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
		{
			std::string line = "energy: " + std::to_string((*it)->energy) + "\n";
			outputfile << line;
		}

		// Close the file
		outputfile.close();
		writeToXyz(seeds, outputFileName + "_sorted_seeds" + ".xyz");

	}
	else {
		//get seed energies from THE SEED FILE
		std::string line;
		std::ifstream file(seedEnergyFile);
		for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
		{
			std::getline(file, line);
			std::string estring = "";
			for (int i = 8; i < line.size(); i++) {
				estring += line[i];
			}
			if (line.size() > 46) (*it)->energy = DBL_MAX;//too large for stof
			//else (*it)->energy = std::stof(estring);//doesnt compine on linux
			if(line.size() <= 46) (*it)->energy = std::stof(estring);

		}
		
	}

	//optimize the top n seeds
	std::vector<structure*> optimizedSeeds = {};
	for (int n = 0; n < noptimizations; n++)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		double minutes = (double)duration.count();
		double timeLeftMinutes = (double)time * 60 - minutes;
		std::vector<double> movement = {};
		int natoms = 0;
		for (int i = 0; i < seeds[n]->elements.size(); i++)
		{
			for (int j = 0; j < seeds[n]->set[i].size(); j++)
			{
				natoms++;
				for(int m = 0; m < 3;m++) movement.push_back(0);
			}
		}
		structure* optimizedSeed = optimizationCalculation(calculator, seeds[n], ++calculations, taskName, optimizationIterations, SCFsteps, rlimit, rlimit, rlimit, movement, structures, basis, method, charge, state,keepFiles, converge, timeLeftMinutes, timeOutMinutes, unrestricted,customBasis);
		if(optimizedSeed != nullptr) optimizedSeeds.push_back(optimizedSeed);
	}
	if (noptimizations > 0)
	{
		std::sort(optimizedSeeds.begin(), optimizedSeeds.end(), energyCompare);
		writeToXyz(optimizedSeeds, outputFileName + "_optimized_seeds" + ".xyz");


		std::ofstream outputfile(outputFileName + "_optimized_seed_energies" + ".txt");
		// Write to the file
		for (std::vector<structure*>::iterator it = optimizedSeeds.begin(); it != optimizedSeeds.end(); it++)
		{
			std::string line = "energy: " + std::to_string((*it)->energy) + "\n";
			outputfile << line;
		}

		// Close the file
		outputfile.close();
	}

	for (std::vector<structure*>::iterator it = optimizedSeeds.begin(); it != optimizedSeeds.end(); it++)
	{
		delete (*it);
	}
	for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
	{
		delete (*it);
	}

	
}
*/
std::vector<structure*> createSeeds( std::vector<std::string> elements, std::vector<int> composition, double cirteriaPercent, std::string xyzout, std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> nRingsNValuesAHoles,int keepPerCombo = -1,int keepPerRingCombination = -1)
{
	std::vector<structure*> structures = {};

	//not used in continue
	structure* s = nullptr;
	double seedEnergy = 0;

	int calculations = 0;

	std::vector<std::pair<int, bool>> symmetries = {};//used in seed mode i & r

	int seedNo = 0;
	std::vector<structure*> seeds = proceduralSeeds(composition, elements, nRingsNValuesAHoles, cirteriaPercent, keepPerCombo,-1,0.01);
	//std::vector<structure*> tempSeeds = proceduralSeeds(composition, elements, numRings[i], nValues, cirteriaPercent);
	std::cout << "seeds size: " << seeds.size() << " writing to " << xyzout << std::endl;
	writeToXyz(seeds, xyzout);
	return seeds;

}

std::pair<std::vector<structure*>,int> sortSeeds(std::vector<structure*> seeds, std::vector<std::string> elements, std::vector<int> composition, std::string energyFile,std::string cancelFile, std::string calculator, double minutesElapsed, double timeOutMinutesEnergy,double time, std::string xyzfilename, double cirteriaPercent, std::string xyzout)
{
	auto start = std::chrono::high_resolution_clock::now();
	int calculations = 0;
	bool stopAlgo = false;
	for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end() && !stopAlgo; it++)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		double minutes = (double)duration.count() + minutesElapsed;
		double timeLeftMinutes = (double)time * 60 - minutes;
		double currentEnergy = energyDFT(calculator, *it, energyFile,cancelFile, timeOutMinutesEnergy, timeLeftMinutes, xyzfilename, ++calculations);
		(*it)->energy = currentEnergy;

		if (timeLeftMinutes < 5.0) stopAlgo = true;//as the calculators stop 5 minutes early, this is causing hundreds of failed calculations to occur at the end of the time limit, this code line will prevent that error.
		if (minutes > double(time) * 60 - 2) stopAlgo = true;
	}
	//sort seeds by energy
	std::cout << "Sorting seeds, time(m): " << (double)time * 60 << std::endl;
	std::sort(seeds.begin(), seeds.end(), energyCompare);
	writeToXyz(seeds, "energy_sorted_" + xyzout);


	//optimize the top n seeds
	return std::pair<std::vector<structure*>, int>(seeds, calculations);
}
std::pair<std::vector<structure*>, int> optimizeSeeds(std::pair<std::vector<structure*>, int> seedsc, std::vector<std::string> elements, std::vector<int> composition, double time, std::string optimizationFile,std::string cancelFile, std::string calculator, double minutesElapsed, double timeOutMinutesOptimization, std::string xyzfilename, double cirteriaPercent, std::string xyzout, int noptimizations)
{
	
	auto start = std::chrono::high_resolution_clock::now();
	bool stopAlgo = false;
	std::vector<structure*> optimizedSeeds = {};
	for (int n = 0; n < noptimizations && n < seedsc.first.size() && !stopAlgo; n++)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		double minutes = (double)duration.count() + minutesElapsed;
		double timeLeftMinutes = (double)time * 60 - minutes;
		std::vector<double> movement = {};
		int natoms = 0;
		for (int i = 0; i < seedsc.first[n]->elements.size(); i++)
		{
			for (int j = 0; j < seedsc.first[n]->set[i].size(); j++)
			{
				natoms++;
				for (int m = 0; m < 3; m++) movement.push_back(0);
			}
		}
		structure* optimizedSeed = optimization(calculator, seedsc.first[n], optimizationFile,cancelFile, timeOutMinutesOptimization, timeLeftMinutes, xyzfilename, ++seedsc.second);
		//structure* optimizedSeed = optimizationCalculation(calculator, seeds[n], ++calculations, taskName, optimizationIterations, SCFsteps, rlimit, rlimit, rlimit, movement, structures, basis, method, charge, state, keepFiles, converge, timeLeftMinutes, timeOutMinutes, unrestricted, customBasis);
		if (optimizedSeed != nullptr)
		{
			optimizedSeed->assignment = seedsc.first[n]->assignment;
			optimizedSeeds.push_back(optimizedSeed);
		}
		else {
			std::cout << "Optimization FAILED" << std::endl;
		}
		
		if (timeLeftMinutes < 5.0) stopAlgo = true;//as the calculators stop 5 minutes early, this is causing hundreds of failed calculations to occur at the end of the time limit, this code line will prevent that error.
		if (minutes > double(time) * 60 - 2) stopAlgo = true;
	}
	if (noptimizations > 0)
	{
		std::sort(optimizedSeeds.begin(), optimizedSeeds.end(), energyCompare);
		writeToXyz(optimizedSeeds, "optimized_" + xyzout);
	}

	/*for (std::vector<structure*>::iterator it = optimizedSeeds.begin(); it != optimizedSeeds.end(); it++)
	{
		delete (*it);
	}
	for (std::vector<structure*>::iterator it = seedsc.first.begin(); it != seedsc.first.end(); it++)
	{
		delete (*it);
	}*/
	return std::pair<std::vector<structure*>, int>(optimizedSeeds, seedsc.second);
}

void gauranteedEscape(std::vector<structure*> seeds,std::vector<std::string> elements,std::vector<int> composition, double time,std::string energyFile,std::string optimizationFile, std::string cancelFile, std::string calculator,double minutesElapsed, double timeOutMinutesEnergy,double timeOutMinutesOptimization,std::string xyzfilename,int previousCalculations, int maximumUphillSteps, double criteriaPercent, int numBranches,double stepSizeAngstroms,std::string xyzout, std::string outputdetails, double angstromCutoff, std::vector<int> handmadeAssignment,bool singleAtomDirections,int nAtomDirections)
{

	/*
	As all structures have been shown to give roughly the same RID value at a given distance in angstroms, here we calculate the RID
	*/
	std::pair<double, double> RIDavgStdev = RIDtoAngstroms(*seeds[0], angstromCutoff, 1000);//this will crash if no seeds are given
	double RIDthreshold = RIDavgStdev.first + RIDavgStdev.second;
	//create all structures to store all optimized structures and prevent exploring duplicates
	std::vector<structure*> allStructures = {};
	
	/*
	
	std::vector<structure*> seeds, - all the seeds with there energies ordered by energy
	std::vector<std::string> elements - all the elements
	std::vector<int> composition, -the number of atoms per element
	std::vector<int> fix, - an asignment of each atom to a class
	double time,
	double maxMovement,
	double eThr,
	std::string energyFile,
	std::string optimizationFile,
	std::string calculator,
	double timeOutMinutesEnergy,
	double timeOutMinutesOptimization,
	std::string xyzfilename: the input to dft calcs
	int previousCalculations, 
	
	maximumUphillSteps: the number of steps in one direction before giving up on the direction
	maximumHops: the number of directions to be explored from each seed. 
	stepSizeAngstroms: the movement per step in angstroms, ie magnitude of direction vector
	*/
	/*
	seeds is an input of optimized ordered seeds. if the seeds have fixed fragments they must also have a filled assigment variable
	max movement is the maximum movement inbetween a step
	eThr is the energy threshold
	*/

	auto start = std::chrono::high_resolution_clock::now();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double minutes = (double)duration.count() + minutesElapsed;
	double timeLeftMinutes = (double)time * 60 - minutes;
	//new algorithm todo
	//work flow
	/*
	1) Get seed(s) -> random symmetry, absolutetly random or from procedure
	1.2) Order seeds by energy and pick the best

	or

	1) Start by randomly ordering fixed fragments with acceptable covalent criteria, then perform optimizations and only move coordinates of fragments that are not fixed in place
	*/
	
	std::vector<structure*> localMinima = {};

	
	bool newSeed = true;
	int seedIterator = 0;
	bool stopAlgo = false;
	int notUnique = 0;


	//efficiency tracking variables
	double GMenergy = DBL_MAX;
	int totalStepsGO = 0;
	int totalStepsSP = 0;
	int GMstepsGO = 0;
	int GMstepsSP = 0;
	double GMtime = DBL_MAX;

	std::stack<std::pair<structure*, std::vector<double>>> theStack = {};
	for (int i = 0; i < seeds.size(); i++) {
		
		//CHECK RID
		bool passedRID = true;
		for (std::vector<structure*>::iterator ast = allStructures.begin(); ast != allStructures.end(); ast++) if (RIDthreshold >= RID(*seeds[i], *(*ast))) passedRID = false;
		//CHECK OVER

		if (passedRID) {
			localMinima.push_back(seeds[i]);
			if (!singleAtomDirections)
			{
				for (int j = 0; j < numBranches; j++) {

					std::vector<double> movement = {};
					
					int dimensions = 0;
					std::vector<int>* assignmentPtr = nullptr;
					if (seeds[i]->assignment.size() == 0) assignmentPtr = &handmadeAssignment;
					else assignmentPtr = &seeds[i]->assignment;
					for (int z = 0; z < assignmentPtr->size(); z++) if (assignmentPtr->size() > dimensions) dimensions = assignmentPtr->size();
					dimensions++;//the first assignment is 0

					if (nAtomDirections < 1) {
						
						for (int id = 0; id < dimensions; id++)
						{
							movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
							movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
							movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
						}
					}
					else {
						std::vector<bool> rncc = randomNchooseC(dimensions, nAtomDirections);
							for (int id = 0; id < dimensions; id++)
							{
								if (rncc[id]) {
									movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
									movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
									movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
								}
								else {
									movement.push_back((0));
									movement.push_back((0));
									movement.push_back((0));
								}
							}
						
					}


					double t2 = 0;
					for (int im = 0; im < movement.size(); im++) t2 += movement[im] * movement[im];
					double normalizationConstant = sqrt(t2);
					for (int im = 0; im < movement.size(); im++) movement[im] = stepSizeAngstroms * movement[im] / normalizationConstant;


					theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movement));
					
				}
				
			}
			else {
				std::cout << "DB__1" << std::endl;
				int dimensions = 0;
				std::vector<int>* assignmentPtr = nullptr;
				if (seeds[i]->assignment.size() == 0) assignmentPtr = &handmadeAssignment;
				else assignmentPtr = &seeds[i]->assignment;
				for (int z = 0; z < assignmentPtr->size(); z++) if (assignmentPtr->size() > dimensions) dimensions = assignmentPtr->size();
				dimensions++;//the first assignment is 0
				dimensions *= 3;
				std::vector<structure*> directions = {};//this will be used to make sure no directions are the same
				std::vector<std::vector<double>> movements = {};//corresponds to directions
				for (int id = 0; id < dimensions; id++)
				{
					std::vector<double> movementPositive = {};
					std::vector<double> movementNegative = {};
					for (int j = 0; j < dimensions; j++)
					{
						if (id != j) {
							movementPositive.push_back(0);
							movementNegative.push_back(0);
						}
						else {
							movementPositive.push_back(stepSizeAngstroms);
							movementNegative.push_back(-stepSizeAngstroms);
						}
					}
					structure* sP = new structure(*seeds[i], movementPositive, *assignmentPtr);
					structure* sN = new structure(*seeds[i], movementNegative, *assignmentPtr);
					bool passedRIDp = true;
					bool passedRIDn = true;
					for (int d = 0; d < directions.size(); d++)
					{
						
						if (RIDthreshold >= RID(*sP, *(directions[d]))) passedRIDp = false;
						if (RIDthreshold >= RID(*sN, *(directions[d]))) passedRIDn = false;

					}
					if (passedRIDp)
					{
						//theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movementPositive));
						movements.push_back(movementPositive);
						directions.push_back(sP);
					}
					else {
						std::cout << "positive direction " << id << " rejected for violating symmetry" << std::endl;
						delete sP;
					}
					if (passedRIDn)
					{
						//theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movementNegative));
						movements.push_back(movementNegative);
						directions.push_back(sN);
					}
					else {
						std::cout << "negative direction " << id << " rejected for violating symmetry" << std::endl;
						delete sN;
					}
				}
				for (int d = 0; d < directions.size(); d++) delete directions[d];
				std::cout << movements.size() << " branches generated from all single atom directions, adding at most " << numBranches << std::endl;
				if (numBranches == 0 || numBranches >= movements.size())
				{
					//ths is a default value, add all of the directions, complete search
					for (int m = 0; m < movements.size(); m++)
					{
						//select randomly from the movements to get the number of desired branches
						theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movements[m]));
					}
				}
				else {
					std::vector<bool> choice = randomNchooseC(movements.size(), numBranches);
					for (int ic = 0; ic < choice.size(); ic++)
					{
						if(choice[ic]) theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movements[ic]));
					}
				}
				
			}
			allStructures.push_back(seeds[i]);
		}
		else {
			std::cout << "rejected optimized structure " << i << " for not being unique" << std::endl;
			notUnique++;
		}
	}

	std::cout << "Total direction structure pairs added to the stack: " << theStack.size() << ",seeds removed for not being unique: " << notUnique << std::endl;

	bool pforp = false;
	while (!theStack.empty() && !stopAlgo)
	{
		std::cout << "GE: Getting local minima & direction vector from the stack" << std::endl;
		std::pair<structure*, std::vector<double>> csd = theStack.top();
		theStack.pop();
		structure* cStruct = csd.first;
		std::vector<double> movement = csd.second;
		if (pforp) {
			std::cout << " top of stack structure with " << cStruct->energy << std::endl;
			cStruct->print();
		}
		//localMinima.push_back(cStruct);//it is already optimized
		double lastEnergy = cStruct->energy;
		double newEnergy = DBL_MAX;

		//make a random directio

		int uphillSteps = 0;
		structure* tempS = cStruct;
		std::vector<structure*> tempStructures = {};
		std::cout << "GE: moving in the energy direction" << std::endl;
		while (newEnergy >= lastEnergy && !stopAlgo && maximumUphillSteps > uphillSteps)
		{
			std::cout << "GE: step:" << uphillSteps << std::endl;
			lastEnergy = tempS->energy;
			std::vector<int>* assignmentPtr = nullptr;
			if (cStruct->assignment.size() == 0) assignmentPtr = &handmadeAssignment;
			else assignmentPtr = &cStruct->assignment;
			//the following line is a bandaid
			//assignmentPtr = &handmadeAssignment;
			if (pforp)
			{
				std::cout << "assignment pointer: " << std::endl;
				for (int iii = 0; iii < assignmentPtr->size(); iii++)
				{
					std::cout << (*assignmentPtr)[iii] << " ";
				}
				std::cout << std::endl;
			}
			tempS = new structure(*tempS, movement, *assignmentPtr);
			if (pforp) {
				std::cout << "moved structure " << std::endl;
				tempS->print();
			}
			//std::cout << " new structure temps:" << std::endl;
			//tempS->print();
			//std::cout << " vs. old structure cstruct: " << std::endl;
			//cStruct->print();
			//std::cout << "Did cstruct even pass the criteria tests?, covalent: " << covalentCriteria(cStruct, criteriaPercent) << " , bonding: " << bondingRequirement(cStruct, criteriaPercent, 1) << std::endl;
			tempStructures.push_back(tempS);
			if (covalentCriteria(tempS, criteriaPercent) && bondingRequirement(tempS, criteriaPercent, 1))
			{
				//make sure that this direction is valid before killing time on it
				newEnergy = energyDFT(calculator, tempS, energyFile,cancelFile, timeOutMinutesEnergy, timeLeftMinutes, xyzfilename, ++previousCalculations);
				totalStepsSP++;
				tempS->energy = newEnergy;
				std::cout << "new energy: " << newEnergy << ", previous: " << lastEnergy << std::endl;
			}
			else {
				newEnergy = DBL_MAX;
				std::cout << "rejected step due to criteria, covalent: " << covalentCriteria(tempS, criteriaPercent) << ", bonding: " << bondingRequirement(tempS, criteriaPercent, 1) << std::endl;
			}
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
			minutes = (double)duration.count() + minutesElapsed;
			timeLeftMinutes = (double)time * 60 - minutes;//time variable is in hours, multiply by 60 to get minutes
			if (timeLeftMinutes < 5.0) stopAlgo = true;//as the calculators stop 5 minutes early, this is causing hundreds of failed calculations to occur at the end of the time limit, this code line will prevent that error.
			if (minutes > double(time) * 60 - 2) stopAlgo = true;
			uphillSteps++;
		}
		//std::cout << minutes << "vs. " << double(time) * 60 - 2 << std::endl;
		
		//std::cout << !stopAlgo <<" " << (newEnergy >= lastEnergy) << " " <<( maximumUphillSteps >= uphillSteps) << std::endl;
		//std::cout << newEnergy << " vs. " << lastEnergy << std::endl;

		if ((maximumUphillSteps > uphillSteps || newEnergy < lastEnergy) && !stopAlgo) {
			std::cout << "GE: Downhill movement detected with " << lastEnergy << " to " << newEnergy <<", optimizing" << std::endl;
			cStruct = optimization(calculator, tempS, optimizationFile,cancelFile, timeOutMinutesOptimization, timeLeftMinutes, xyzfilename, ++previousCalculations);
			totalStepsGO++;
			if (cStruct != nullptr) {
				cStruct->time = minutes;
				//localMinima.push_back(cStruct);
				std::cout << "new local minima optimized with energy: " << cStruct->energy << " , calculation: " << previousCalculations << std::endl;
				if (pforp)
				{
					cStruct->print();
				}
				//CHECK RID
				bool passedRID = true;
				for (std::vector<structure*>::iterator ast = allStructures.begin(); ast != allStructures.end(); ast++) if (RIDthreshold >= RID(*cStruct, *(*ast))) passedRID = false;
				if (passedRID) {
					localMinima.push_back(cStruct);
					if (cStruct->energy < GMenergy)
					{
						GMenergy = cStruct->energy;
						GMstepsGO = totalStepsGO;
						GMstepsSP = totalStepsSP;
						GMtime = cStruct->time;
					}
					if (!singleAtomDirections)
					{
						for (int j = 0; j < numBranches; j++) {
							std::vector<double> movement = {};
							int dimensions = 0;
							std::vector<int>* assignmentPtr = nullptr;
							if (cStruct->assignment.size() == 0) assignmentPtr = &handmadeAssignment;
							else assignmentPtr = &cStruct->assignment;
							for (int z = 0; z < assignmentPtr->size(); z++) if (assignmentPtr->size() > dimensions) dimensions = assignmentPtr->size();
							dimensions++;//the first assignment is 0
							for (int id = 0; id < dimensions; id++)
							{
								movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
								movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
								movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
							}


							double t2 = 0;
							for (int im = 0; im < movement.size(); im++) t2 += movement[im] * movement[im];
							double normalizationConstant = sqrt(t2);
							for (int im = 0; im < movement.size(); im++) movement[im] = stepSizeAngstroms * movement[im] / normalizationConstant;


							theStack.push(std::pair<structure*, std::vector<double>>(cStruct, movement));

						}

					}
					else {

						int dimensions = 0;
						std::vector<int>* assignmentPtr = nullptr;
						if (cStruct->assignment.size() == 0) assignmentPtr = &handmadeAssignment;
						else assignmentPtr = &cStruct->assignment;
						for (int z = 0; z < assignmentPtr->size(); z++) if (assignmentPtr->size() > dimensions) dimensions = assignmentPtr->size();
						dimensions++;//the first assignment is 0
						dimensions *= 3;
						std::vector<structure*> directions = {};//this will be used to make sure no directions are the same
						std::vector<std::vector<double>> movements = {};//corresponds to directions

						for (int id = 0; id < dimensions; id++)
						{
							std::vector<double> movementPositive = {};
							std::vector<double> movementNegative = {};
							for (int j = 0; j < dimensions; j++)
							{
								if (id != j) {
									movementPositive.push_back(0);
									movementNegative.push_back(0);
								}
								else {
									movementPositive.push_back(stepSizeAngstroms);
									movementNegative.push_back(-stepSizeAngstroms);
								}
							}
							structure* sP = new structure(*cStruct, movementPositive, *assignmentPtr);
							structure* sN = new structure(*cStruct, movementNegative, *assignmentPtr);
							bool passedRIDp = true;
							bool passedRIDn = true;
							for (int d = 0; d < directions.size(); d++)
							{

								if (RIDthreshold >= RID(*sP, *(directions[d]))) passedRIDp = false;
								if (RIDthreshold >= RID(*sN, *(directions[d]))) passedRIDn = false;

							}
							if (passedRIDp)
							{
								//theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movementPositive));
								movements.push_back(movementPositive);
								directions.push_back(sP);
							}
							else {
								std::cout << "positive direction " << id << " rejected for violating symmetry" << std::endl;
								delete sP;
							}
							if (passedRIDn)
							{
								//theStack.push(std::pair<structure*, std::vector<double>>(seeds[i], movementNegative));
								movements.push_back(movementNegative);
								directions.push_back(sN);
							}
							else {
								std::cout << "negative direction " << id << " rejected for violating symmetry" << std::endl;
								delete sN;
							}
						}
						for (int d = 0; d < directions.size(); d++) delete directions[d];
						std::cout << movements.size() << " branches generated from all single atom directions, adding at most " << numBranches << std::endl;
						if (numBranches == 0 || numBranches >= movements.size())
						{
							//ths is a default value, add all of the directions, complete search
							for (int m = 0; m < movements.size(); m++)
							{
								//select randomly from the movements to get the number of desired branches
								theStack.push(std::pair<structure*, std::vector<double>>(cStruct, movements[m]));
							}
						}
						else {
							std::vector<bool> choice = randomNchooseC(movements.size(), numBranches);
							for (int ic = 0; ic < choice.size(); ic++)
							{
								if (choice[ic]) theStack.push(std::pair<structure*, std::vector<double>>(cStruct, movements[ic]));
							}
						}

					}
					allStructures.push_back(cStruct);
				}
				else {
					std::cout << "optimized structure is similar to a structure already in the stack and will not be explored" << std::endl;
				}
			}



		}else { std::cout << "GE: no downhill movement, going to the top of the stack" << std::endl; }
		for (int ts = 0; ts < tempStructures.size(); ts++) delete tempStructures[ts];


		//check if we need to stop early
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutes = (double)duration.count() + minutesElapsed;
		timeLeftMinutes = (double)time * 60 - minutes;
		if (timeLeftMinutes < 5.0) stopAlgo = true;//same reason as in identical line above. see identical lines above for comments. line below may be deprecated.
		if (minutes > double(time) * 60 - 2) stopAlgo = true;
	}


	/*
	TODO implement DFS vs BFS stacks for directions
	std::queue<std::pair<structure*, std::vector<double>>> dQueue;
	*/
	std::cout << "GE: printing results" << std::endl;

	//Print the output
	std::sort(localMinima.begin(), localMinima.end(), energyCompare);
	writeToXyz(localMinima, xyzout);
	std::ofstream outFile(outputdetails);
	if (!outFile.is_open()) {
		std::cout << "Error: Unable to open file " << outputdetails<< " for writing." << std::endl;
		return;
	}
	else {
		if (theStack.empty()) {
			outFile << "Search completed\n";
		}
		else {
			outFile << "Unexplored structures: " << std::endl;
			while (!theStack.empty())
			{
				outFile << theStack.top().first->xyz() << "\n";
				theStack.pop();
			}
		}

	}
	outFile.close();
	std::cout << "GM energy, GM time, No. SP calculations at GM, No. GO calculations at GM: " << std::endl;
	std::cout << GMenergy << std::endl;
	std::cout << GMtime << std::endl;
	std::cout << GMstepsSP << std::endl;
	std::cout << GMstepsGO << std::endl;
	

}

void BasinHopping(std::vector<structure*> seeds, std::vector<std::string> elements, std::vector<int> composition, double time, std::string energyFile, std::string optimizationFile, std::string cancelFile, std::string calculator, double minutesElapsed, double timeOutMinutesEnergy, double timeOutMinutesOptimization, std::string xyzfilename, int previousCalculations, int maximumUphillSteps, double criteriaPercent, double stepSizeAngstroms, std::string xyzout, double angstromCutoff, std::vector<int> handmadeAssignment, bool filters)
{
	//using single atom directions to indicate RID filters are on
	/*
	As all structures have been shown to give roughly the same RID value at a given distance in angstroms, here we calculate the RID
	*/
	std::pair<double, double> RIDavgStdev = RIDtoAngstroms(*seeds[0], angstromCutoff, 1000);//this will crash if no seeds are given
	double RIDthreshold = RIDavgStdev.first + RIDavgStdev.second;
	//create all structures to store all optimized structures and prevent exploring duplicates
	std::vector<structure*> allStructures = {};

	/*

	std::vector<structure*> seeds, - all the seeds with there energies ordered by energy
	std::vector<std::string> elements - all the elements
	std::vector<int> composition, -the number of atoms per element
	std::vector<int> fix, - an asignment of each atom to a class
	double time,
	double maxMovement,
	double eThr,
	std::string energyFile,
	std::string optimizationFile,
	std::string calculator,
	double timeOutMinutesEnergy,
	double timeOutMinutesOptimization,
	std::string xyzfilename: the input to dft calcs
	int previousCalculations,

	maximumUphillSteps: the number of steps in one direction before giving up on the direction
	maximumHops: the number of directions to be explored from each seed.
	stepSizeAngstroms: the movement per step in angstroms, ie magnitude of direction vector
	*/
	/*
	seeds is an input of optimized ordered seeds. if the seeds have fixed fragments they must also have a filled assigment variable
	max movement is the maximum movement inbetween a step
	eThr is the energy threshold
	*/

	auto start = std::chrono::high_resolution_clock::now();

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
	double minutes = (double)duration.count() + minutesElapsed;
	double timeLeftMinutes = (double)time * 60 - minutes;
	//new algorithm todo
	//work flow
	/*
	1) Get seed(s) -> random symmetry, absolutetly random or from procedure
	1.2) Order seeds by energy and pick the best

	or

	1) Start by randomly ordering fixed fragments with acceptable covalent criteria, then perform optimizations and only move coordinates of fragments that are not fixed in place
	*/

	std::vector<structure*> localMinima = {};


	bool newSeed = true;
	int seedIterator = 0;
	bool stopAlgo = false;
	int notUnique = 0;

	//using seeed 1 only
	int totalStepsSP = 0;
	int totalStepsGO = 0;
	int GMstepsSP = 0;
	int GMstepsGO = 0;
	double GMtime = DBL_MAX;
	double GMenergy = DBL_MAX;
	structure* cStruct = seeds[0];
	std::vector<double> eChange = {};
	while (!stopAlgo)
	{
		std::cout << "GE: Getting local minima & direction vector from the stack" << std::endl;
		std::vector<double> movement = {};
		int dimensions = 0;
		std::vector<int>* assignmentPtr = nullptr;
		if (cStruct->assignment.size() == 0) assignmentPtr = &handmadeAssignment;
		else assignmentPtr = &cStruct->assignment;
		for (int z = 0; z < assignmentPtr->size(); z++) if (assignmentPtr->size() > dimensions) dimensions = assignmentPtr->size();
		dimensions++;//the first assignment is 0
		for (int id = 0; id < dimensions; id++)
		{
			movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
			movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
			movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
		}


		double t2 = 0;
		for (int im = 0; im < movement.size(); im++) t2 += movement[im] * movement[im];
		double normalizationConstant = sqrt(t2);
		for (int im = 0; im < movement.size(); im++) movement[im] = stepSizeAngstroms * movement[im] / normalizationConstant;


		double lastEnergy = cStruct->energy;
		double newEnergy = DBL_MAX;
		structure* tempS = new structure(*cStruct, movement, *assignmentPtr);

		
		structure* oS = nullptr;
		if (filters)
		{
			bool acceptStructure = true;
			bool covalent = true;
			bool bonding = true;
			bool uniqueness = true;
			if (!covalentCriteria(tempS, criteriaPercent)) {
				acceptStructure = false;
				covalent = false;
			}
			if (!bondingRequirement(tempS, criteriaPercent, 1)) {
				acceptStructure = false;
				bonding = false;
			}
			bool passedRID = true;
			for (std::vector<structure*>::iterator ast = allStructures.begin(); ast != allStructures.end(); ast++) if (RIDthreshold >= RID(*tempS, *(*ast))) passedRID = false;
			if (!passedRID) {
				acceptStructure = false;
				uniqueness = false;
			}
			if (acceptStructure)
			{
				oS = optimization(calculator, tempS, optimizationFile,cancelFile, timeOutMinutesOptimization, timeLeftMinutes, xyzfilename, ++previousCalculations);
				totalStepsGO++;
			}
			else {
				oS = nullptr;
				std::cout << "structure filtered out due to c:" << covalent << "b:" << bonding << "u:" << uniqueness << std::endl;
			}
			allStructures.push_back(tempS);
		}
		else {

			allStructures.push_back(tempS);
			oS = optimization(calculator, tempS, optimizationFile,cancelFile, timeOutMinutesOptimization, timeLeftMinutes, xyzfilename, ++previousCalculations);
			totalStepsGO++;
		}
		if (oS != nullptr) {
			newEnergy = oS->energy;
			
			
			
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
			minutes = (double)duration.count() + minutesElapsed;

			oS->time = minutes;

			if (newEnergy < GMenergy)
			{
				GMenergy = newEnergy;
				GMstepsGO = totalStepsGO;
				GMtime = oS->time;
			}

			localMinima.push_back(oS);
			allStructures.push_back(oS);

			/*
		Figure out the kbteff
		*/
			eChange.push_back(lastEnergy - newEnergy);
			std::sort(eChange.begin(), eChange.end());
			double ecthr = 0;
			if (eChange.size() % 2 == 0)
			{
				ecthr = (eChange[eChange.size() / 2] + eChange[eChange.size() / 2 - 1]) / 2;
			}
			else {
				ecthr = eChange[eChange.size() / 2 - 1];
			}
			//ln(0.5) = ecthr/kbteff
			double kbteff = ecthr / std::log(0.5);
			std::cout << "set kbteff = " << kbteff << std::endl;


			double p = exp((lastEnergy - newEnergy) / kbteff);

			if (rand() / RAND_MAX < p || newEnergy < lastEnergy)
			{
				std::cout << "structure accepted. moving on";
				cStruct = oS;
			}
			else {
				std::cout << "structure rejected";
			}

			//bool passedRID = true;
			//for (std::vector<structure*>::iterator ast = allStructures.begin(); ast != allStructures.end(); ast++) if (RIDthreshold >= RID(*cStruct, *(*ast))) passedRID = false;

		}
		
			
	
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutes = (double)duration.count() + minutesElapsed;
		timeLeftMinutes = (double)time * 60 - minutes;
		if (minutes > double(time) * 60 - 2) stopAlgo = true;

		
	}


	
	std::sort(localMinima.begin(), localMinima.end(), energyCompare);
	writeToXyz(localMinima, xyzout);
	std::cout << "GM energy, GM time, No. SP calculations at GM, No. GO calculations at GM: " << std::endl;
	std::cout << GMenergy << std::endl;
	std::cout << GMtime << std::endl;
	std::cout << GMstepsSP << std::endl;
	std::cout << GMstepsGO << std::endl;
}


void sharedFeatureSearch()
{
	//workflow
	/*
	1) Generate n seeds, get the energies and one optimization.
	*/

	/*
	2) search for similarities among subclusters
	2.1) report the most similar subclusters
	*/

	/*
	3) generate and search holding subclusters to be constant
	*/
}


int main(int argv, char* argc[])
{

	/*
	
	
	*/
	
	//
	srand(time(0));
	/*double stepSizeAngstroms = 0.8;
	std::vector<structure*> seedsssss = createSeeds({ "B" }, { 6 }, 30, "dm.xyz", { std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>({ 2 }, { 5 }, { 0 }) }, true, -1);
	
	std::vector<double> movement = {};
	for (int id = 0; id < 10*3; id++)
	{
		movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
		movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
		movement.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)));
	}


	double t2 = 0;
	for (int im = 0; im < movement.size(); im++) t2 += movement[im] * movement[im];
	double normalizationConstant = sqrt(t2);
	for (int im = 0; im < movement.size(); im++) movement[im] = stepSizeAngstroms * movement[im] / normalizationConstant;


	structure* tempS = new structure(*seedsssss[0], movement, { 0,1,2,3,4,5,6,7,8,9 });
	std::pair<double, double> RIDavgStdev = RIDtoAngstroms(*seedsssss[0], 0.2, 1000);
	double RIDthreshold = RIDavgStdev.first + RIDavgStdev.second;
	double RIDvalue = RID(*seedsssss[0],*tempS);
	std::cout << RIDthreshold << " < " << RIDvalue << std::endl;
	return 0;*/
	/*
	std::ofstream ofile("rid.csv");
	for (float d = 0.05; d < 10; d += 0.05) {
		std::string cline = std::to_string(d);
		std::cout << d << std::endl;
		for (int i = 4; i < 11; i++)
		{
			
			std::string fname = "b" + std::to_string(i) + ".xyz";
			std::vector<structure*> gm = structuresFromXYZ(fname);
		
				std::pair<double, double> RIDavgStdev = RIDtoAngstroms(*gm[0], d, 1000);//this will crash if no seeds are given
				double RIDthreshold = RIDavgStdev.first + RIDavgStdev.second;
				double RIDvalue = RID(*gm[0], *gm[1]);
				std::cout << "rid value " << RIDvalue << std::endl;
				std::cout << "rid threshold: " << RIDthreshold << std::endl;
		
			//std::vector<structure*> ogm = structuresFromXYZ("bout.xyz");
				cline += "," + std::to_string(RIDthreshold);
		
		
		}
		ofile << cline << std::endl;
	}
	ofile.close();
	return 0;
	*/
	/*rand(time(0));
	std::vector<structure*> gm = structuresFromXYZ("ruout.xyz");
	std::cout << covalentCriteria(gm[0], 30);
	return 0;*/
	/*
	std::vector<structure*> gm = structuresFromXYZ("b6.xyz");
	//std::vector<structure*> ogm = structuresFromXYZ("bout.xyz");

	std::pair<double, double> RIDavgStdev = RIDtoAngstroms(*gm[0], 0.2, 1000);//this will crash if no seeds are given
	double RIDthreshold = RIDavgStdev.first + RIDavgStdev.second;
	double RIDvalue = RID(*gm[0], *gm[1]);
	std::cout << "rid value " << RIDvalue << std::endl;
	std::cout << "rid threshold: " << RIDthreshold << std::endl;
	return 0;
	
	srand(time(0));
	*/
	//optimization("Gaussian", structuresFromXYZ("energy_sorted_la2b4out.xyz")[0], "la2b4O", 1000, 10000, "res.xyz", 42)->writeToXyz("hansae.xyz");
	//return 0;
	//std::vector<structure*> ses = structuresFromXYZ("agb9seeds.xyz");

	/*
	int width = 8;
	for (std::vector<structure*>::iterator it = ses.begin(); it != ses.end(); it++) {
		std::pair<double, double> val = RIDtoAngstroms(*ses[0], 2.05, 100);
		std::cout << "  RID:" << std::setw(width) << val.first << "  stdev:" << std::setw(width) << val.second << std::endl;
	}*/
	/*int width = 8;
	for (double a = 0; a < 100; a += 0.05)
	{
		std::pair<double, double> val = RIDtoAngstroms(*ses[0], a, 100);
		std::cout << "dist:" << std::setw(width) << a << "  RID:" << std::setw(width) << val.first << "  stdev:" << std::setw(width) << val.second << std::endl;
	}*/
	//return 0;
	/*
	std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> nrnva = {
		//std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>({2}, {5},{0}),
	std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>({3,2}, {8,1},{0,1,7}),
	std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>({3}, {5,1},{0,1,3})
	};
	std::vector<structure*> seeds1 = createSeeds({ "Ag","B" }, { 1,9 },10, "agb9seeds.xyz", nrnva, 20, -1);
	return 0;
	*/
	std::vector<std::string> elements = {};
	std::vector<int> composition = {};
	std::string energyFile = "";
	std::string optimizationFile = "";
	std::string cancelFile = "script.py";
	std::string xyzinp = "input.xyz";
	std::string xyzout = "out.xyz";
	double criteriaPercent = 30;
	std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> nRingsNValuesAHoles;
	int keep = 0;
	std::string calculator = "P";
	int noptimizations = -1;
	int nAtomDirections = 0;

	double timeOutMinutesEnergy;
	double timeOutMinutesOptimization;
	double time = 47.5;

	double angstromLimit = 0;

	double stepSize = 0.4;
	int maxUphillSteps = 20;
	int nBranches = 6;

	std::string energySeedsFile = "";
	std::string optimizedSeedsFile = "";
	std::vector<int> assignment;

	std::ifstream inputFile(argc[1]);

	bool printInp = true;
	bool singleAtomDirections = false;
	//std::ifstream inputFile("inp.txt");
	bool analysis = false;
	bool bhflag = false;
	double kbteff = 1;
	int rseed = 0;
	
	std::string surface = "";//file for the structure of the fixed surface
	//surface is assumed to be relative to the clusters center at 0,0,0
	

	std::string read;
	while (std::getline(inputFile, read))
	{
		if(read[read.size()-1] == '\r') read.erase(read.size() - 1);//remove unwanted \r

		//std::cout << "read: " << read << std::endl;
		if (read.substr(0, std::string("composition").size()) == "composition")
		{
			std::getline(inputFile, read);
			if (read[read.size() - 1] == '\r') read.erase(read.size() - 1);//remove unwanted \r
			do {
				std::string num = "";
				bool swap = false;
				std::string element = "";
				for (int c = 0; c < read.size(); c++)
				{
					if (read[c] == ' ') swap = true;
					else if (swap) element += read[c];
					else num += read[c];
				}
				elements.push_back(element);
				composition.push_back(std::stoi(num));

				std::getline(inputFile, read);
				if (read[read.size() - 1] == '\r') read.erase(read.size() - 1);//remove unwanted \r
			} while (read.substr(0, std::string("end").size() ) != "end");
			if (printInp) {
				std::cout << "elements:";
				for (std::vector<std::string>::iterator it = elements.begin(); it != elements.end(); it++) std::cout << *it << ",";
				std::cout << std::endl;
				std::cout << "composition:";
				for (std::vector<int>::iterator it = composition.begin(); it != composition.end(); it++) std::cout << *it << ",";
				std::cout << std::endl;
			}

		}
		else if (read.substr(0, std::string("surface").size()) == "surface")
		{
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) surface += read[c];
			}
			if (printInp) std::cout << "surface structure file:" << surface << " containing coordinates relative to the clusters starting center at 0,0,0 " << std::endl;
		}
		else if (read.substr(0, std::string("criteria").size()) == "criteria")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			criteriaPercent = std::stod(tempString);
			if (printInp) std::cout << "criteria percent: " << criteriaPercent << std::endl;
		}
		else if (read.substr(0, std::string("xyzinp").size()) == "xyzinp")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				//std::cout << "read size: " << read.size() << " c: " << c << std::endl;
				//std::cout << "readc: " << read[c] << std::endl;
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			//std::cout << "here" << std::endl;
			//std::cout << tempString << std::endl;
			xyzinp = tempString;
			//std::cout << tempString << " success" <<  std::endl;
			if (printInp) std::cout << "xyzinp: " << xyzinp << std::endl;
		}
		else if (read.substr(0, std::string("xyzout").size()) == "xyzout")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			xyzout = tempString;
			if (printInp) std::cout << "xyzout: " << xyzout << std::endl;
		}
		else if (read.substr(0, std::string("rva").size()) == "rva")
		{
			std::cout << "RVA" << std::endl;
			std::getline(inputFile, read);
			if (read[read.size() - 1] == '\r') read.erase(read.size() - 1);//remove unwanted \r
			do {
				std::string cnum = "";
				std::vector<int> r = {};
				std::vector<int> v = {};
				std::vector<int> a = {};

				int swap = 0;
				for (int c = 0; c < read.size(); c++)
				{
					std::cout << swap << " " << read[c] << std::endl;
					if (read[c] == ' ' || read[c] == ',')
					{
						switch (swap)
						{
						case 0:
							r.push_back(std::stoi(cnum));
							cnum = "";
							break;
						case 1:
							v.push_back(std::stoi(cnum));
							cnum = "";
							break;
						case 2:
							a.push_back(std::stoi(cnum));
							cnum = "";
							break;
						default:
							break;
						}
					}
					if (read[c] == ',') swap++;

					if (read[c] != ',' && read[c] != ' ' && read[c] != '\n' && read[c] != '\r') cnum += read[c];
				}
				a.push_back(std::stoi(cnum));
				std::cout << r[0] << " " << v[0] << " " << a[0] << std::endl;
				nRingsNValuesAHoles.push_back(std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>(r, v, a));
				std::getline(inputFile, read);
				if (read[read.size() - 1] == '\r') read.erase(read.size() - 1);//remove unwanted \r

			} while (read.substr(0, std::string("end").size()) != "end");
			if (printInp) {
				std::cout << "nrvah:" << std::endl;
				for (std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>>::iterator it = nRingsNValuesAHoles.begin(); it != nRingsNValuesAHoles.end(); it++) {
					std::cout << "Rings: ";
					for (std::vector<int>::iterator itt = std::get<0>(*it).begin(); itt != std::get<0>(*it).end(); itt++) std::cout << *itt << " ";
					std::cout << "\nValues:";
					for (std::vector<int>::iterator itt = std::get<1>(*it).begin(); itt != std::get<1>(*it).end(); itt++) std::cout << *itt << " ";
					std::cout << "\nHoles:";
					for (std::vector<int>::iterator itt = std::get<2>(*it).begin(); itt != std::get<2>(*it).end(); itt++) std::cout << *itt << " ";
					std::cout << std::endl;
				}
				std::cout << std::endl;
			}
		}
		else if (read.substr(0, std::string("keep").size()) == "keep")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			keep = std::stoi(tempString);
			if (printInp) std::cout << "keep:" << keep << std::endl;
		}
		else if (read.substr(0, std::string("energyFile").size()) == "energyFile")
		{
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) energyFile += read[c];
			}
			if(printInp) std::cout << "energy file:" << energyFile << std::endl;
		}
		else if (read.substr(0, std::string("optimizationFile").size()) == "optimizationFile")
		{
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) optimizationFile += read[c];
			}
			if (printInp) std::cout << "optimization files:" << optimizationFile << std::endl;
		}
		else if (read.substr(0, std::string("pythonFile").size()) == "pythonFile")
		{
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) cancelFile += read[c];
		}
		if (printInp) std::cout << "energy file:" << energyFile << std::endl;
		}
		else if (read.substr(0, std::string("rseed").size()) == "rseed")
		{
		bool swap = false;
		std::string randseed = "";
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) randseed += read[c];
		}
		rseed = std::stoi(randseed);
		}
		else if (read.substr(0, std::string("seFile").size()) == "seFile")
		{
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) energySeedsFile += read[c];
		}
		if (printInp) std::cout << "seed energy file:" << energySeedsFile << std::endl;
		}
		else if (read.substr(0, std::string("oeFile").size()) == "oeFile")
		{
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) optimizedSeedsFile += read[c];
		}
		if (printInp) std::cout << "optimized energy file:" << optimizedSeedsFile << std::endl;
		}
		else if (read.substr(0, std::string("calculator").size()) == "calculator")
		{
			calculator = "";//currently is P
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) calculator += read[c];
			}
			if (printInp) std::cout << "calculator:" << calculator << std::endl;
		}
		else if (read.substr(0, std::string("timeOutE").size()) == "timeOutE")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			timeOutMinutesEnergy = std::stod(tempString);
			if (printInp) std::cout << "timeOutMinutesEnergy:" << timeOutMinutesEnergy << std::endl;
		}
		else if (read.substr(0, std::string("timeOutO").size()) == "timeOutO")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			timeOutMinutesOptimization = std::stod(tempString);
			if (printInp) std::cout << "timeOutMinutesOptimiztion:" << timeOutMinutesOptimization << std::endl;
		}
		else if (read.substr(0, std::string("time").size()) == "time")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			time = std::stod(tempString);
			if (printInp) std::cout << "time: " << time << std::endl;
		}
		else if (read.substr(0, std::string("unique").size()) == "unique")
		{
		std::string tempString = "";
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) tempString += read[c];
		}
		angstromLimit = std::stod(tempString);
		if (printInp) std::cout << "unique structure angstrom limit: " << angstromLimit << std::endl;
		}
		else if (read.substr(0, std::string("criteria").size()) == "criteria")
		{
		std::string tempString = "";
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) tempString += read[c];
		}
		criteriaPercent = std::stod(tempString);
		if (printInp) std::cout << "cirteria percent: " << criteriaPercent << std::endl;
		}
		else if (read.substr(0, std::string("nopt").size()) == "nopt")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			noptimizations = std::stoi(tempString);
			if (printInp) std::cout << "noptimizations: " << noptimizations << std::endl;
		}
		else if (read.substr(0, std::string("stepsize").size()) == "stepsize")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			stepSize = std::stod(tempString);
			if (printInp) std::cout << "step size: " << stepSize << std::endl;
		}
		else if (read.substr(0, std::string("uphill").size()) == "uphill")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			maxUphillSteps = std::stoi(tempString);
			if (printInp) std::cout << "maxUphillSteps: " << maxUphillSteps << std::endl;
		}
		else if (read.substr(0, std::string("branches").size()) == "branches")
		{
			std::string tempString = "";
			bool swap = false;
			for (int c = 0; c < read.size(); c++)
			{
				if (read[c] == ' ') swap = true;
				else if (swap) tempString += read[c];
			}
			nBranches = std::stoi(tempString);
			if (printInp) std::cout << "nbranches: " << nBranches << std::endl;
		}
		else if (read.substr(0, std::string("assignment").size()) == "assignment")
		{
		std::string tempString = "";
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (!swap)
			{
				if (read[c] == ' ') swap = true;
			}
			else {
			//std::cout << c << ":" << read[c] << std::endl;
			if (read[c] == ' ') {
				if (tempString.size() != 0) assignment.push_back(std::stoi(tempString));
				tempString = "";
			}
			else tempString += read[c];
			}
		}
		assignment.push_back(std::stoi(tempString));
		if (printInp) {
			std::cout << "assignment: ";
			for (std::vector<int>::iterator at = assignment.begin(); at != assignment.end(); at++) std::cout << *at << " ";
			std::cout << std::endl;
		}
		}
		else if (read.substr(0, std::string("sad").size()) == "sad")
		{
		singleAtomDirections = true;
		if (printInp) std::cout << "single atom directions specified" << std::endl;
		}
		else if (read.substr(0, std::string("nad").size()) == "nad")
		{
		std::string tempString = "";
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) tempString += read[c];
		}
		nAtomDirections = std::stoi(tempString);
		if (printInp) std::cout << "n atom directions: " << nAtomDirections << std::endl;
		}
		else if (read.substr(0, std::string("bh").size()) == "bh")
		{
		bhflag = true;
		if (printInp) std::cout << "basin hopping set" << std::endl;
		}
		else if (read.substr(0, std::string("kbteff").size()) == "kbteff")
		{
		std::string tempString = "";
		bool swap = false;
		for (int c = 0; c < read.size(); c++)
		{
			if (read[c] == ' ') swap = true;
			else if (swap) tempString += read[c];
		}
		kbteff = std::stod(tempString);
		if (printInp) std::cout << "kbteff: " << kbteff << std::endl;
		}
		else if (read.substr(0, std::string("analysis").size()) == "analysis")
		{
		analysis = true;
		if (printInp) std::cout << "analysis on previous calculations specified" << std::endl;
		}
		else {
		std::cout << "unknown line:" <<  read << std::endl;
		}


	}
	if (rseed != 0) srand(rseed);
	if (!analysis) {
		auto start = std::chrono::high_resolution_clock::now();
		std::cout << "input read" << std::endl;

		std::vector<structure*> seeds = {};
		for (int i = 0; i < nRingsNValuesAHoles.size(); i++)
		{
			std::cout << "ring combination" << std::endl;
			std::cout << std::get<0>(nRingsNValuesAHoles[i])[0] << " " << std::get<1>(nRingsNValuesAHoles[i])[0] << " " << std::get<2>(nRingsNValuesAHoles[i])[0] << std::endl;
		}
		for (int i = 0; i < elements.size(); i++)
		{
			std::cout << "element: " << elements[i] << " " << composition[i] << std::endl;
		}
		if (optimizedSeedsFile == "" && energySeedsFile == "")seeds = createSeeds(elements, composition, criteriaPercent, xyzout, nRingsNValuesAHoles, keep, -1);

		std::cout << "seeds made: " << seeds.size() << std::endl;
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		double minutesElpased = (double)duration.count();

		std::pair<std::vector<structure*>, int> sortedSeeds;
		if (optimizedSeedsFile != "");//do nothing
		else if (energySeedsFile == "")sortedSeeds = sortSeeds(seeds, elements, composition, energyFile,cancelFile, calculator, minutesElpased, timeOutMinutesEnergy, time, xyzinp, criteriaPercent, xyzout);
		else sortedSeeds = std::pair<std::vector<structure*>, int>(structuresFromXYZ(energySeedsFile), 0);
		std::cout << "seeds sorted by energy" << std::endl;

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutesElpased = (double)duration.count();

		std::pair<std::vector<structure*>, int> optimizedSeeds;
		if (optimizedSeedsFile == "")optimizedSeeds = optimizeSeeds(sortedSeeds, elements, composition, time, optimizationFile,cancelFile, calculator, minutesElpased, timeOutMinutesOptimization, xyzinp, criteriaPercent, xyzout, noptimizations);
		else optimizedSeeds = std::pair<std::vector<structure*>, int>(structuresFromXYZ(optimizedSeedsFile), 0);
		/*
		Delete sorted seeds, they will not be used
		*/
		for (std::vector<structure*>::iterator st = sortedSeeds.first.begin(); st != sortedSeeds.first.end(); st++) delete* st;

		std::cout << optimizedSeeds.first.size() << " seeds optimized" << std::endl;

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
		minutesElpased = (double)duration.count();

		//if assignment is not specififed but energy or optimized seeds are provided the software will crash, therefore a fake assignment should be made where each atom gets its own assignment.
		if ((energySeedsFile != "" || optimizedSeedsFile != "") && assignment.size() == 0)
		{
			int a = 0;
			for (int i = 0; i < composition.size(); i++) for (int j = 0; j < composition[i]; j++) assignment.push_back(a++);
		}

		if(!bhflag) gauranteedEscape(optimizedSeeds.first, elements, composition, time, energyFile, optimizationFile,cancelFile, calculator, minutesElpased, timeOutMinutesEnergy, timeOutMinutesOptimization, xyzinp, optimizedSeeds.second, maxUphillSteps, criteriaPercent, nBranches, stepSize, xyzout, "output.txt", angstromLimit, assignment, singleAtomDirections,nAtomDirections);
		else BasinHopping(optimizedSeeds.first, elements, composition, time, energyFile, optimizationFile,cancelFile, calculator, minutesElpased, timeOutMinutesEnergy, timeOutMinutesOptimization, xyzinp, optimizedSeeds.second, maxUphillSteps, criteriaPercent, stepSize, xyzout, angstromLimit, assignment, singleAtomDirections);
		/*
		Delete optimized seeds, the software is done
		*/
		for (std::vector<structure*>::iterator st = optimizedSeeds.first.begin(); st != optimizedSeeds.first.end(); st++) delete* st;

	}
	else {
		double RIDthreshold = DBL_MAX;
		int calcNo = 1;
		std::vector<structure*> structures = {};
		bool stop = false;
		while (!stop)
		{
			std::string opath = "calc" + std::to_string(calcNo) + "/" + optimizationFile;
			std::string epath = "calc" + std::to_string(calcNo) + "/" + energyFile;
			if (calculator == "ADF")
			{
				opath += ".out";
				epath += ".out";
			}
			else if (calculator == "Gaussian")
			{
				opath += ".log";
				epath += ".log";
			}
			std::ifstream fileo(opath);
			std::ifstream filee(epath);
			
			if (fileo) {
				std::cout << "File exists." << std::endl;
				fileo.close();
				structure* tempStruct = readOptFile(calculator, optimizationFile, calcNo, composition, elements);
				if (tempStruct != nullptr) {
					bool passedRID = true;
					for (std::vector<structure*>::iterator ast = structures.begin(); ast != structures.end(); ast++) if (RIDthreshold >= RID(*tempStruct, *(*ast))) passedRID = false;
					if (passedRID) {
						if (structures.size() == 0)
						{
							std::pair<double, double> RIDavgStdev = RIDtoAngstroms(*tempStruct, angstromLimit, 1000);
							RIDthreshold = RIDavgStdev.first + RIDavgStdev.second;
						}
						structures.push_back(tempStruct);
					}
					else {
						delete tempStruct;
					}
				}
				calcNo += 1;
			}
			else if(filee) {
				std::cout << "Optimization file does not exist, but energy does" << std::endl;
				calcNo += 1;
			}
			else {
				stop = true;
				std::cout << "came to the end of " << calcNo << " calculations" << std::endl;
			}

			
			
		}
		std::sort(structures.begin(), structures.end(), energyCompare);
		writeToXyz(structures, "analysisResults.xyz");
	}
	
	
}

