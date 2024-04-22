#ifndef chemH
#define chemH

#include <vector>
#include <functional>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <set>
#include <string>
#include <map>
#include <queue>
# define M_PI           3.14159265358979323846

class atom {
public:
	atom(double x, double y, double z);
	atom();
	atom(const atom& t);
	atom operator=(const atom& t);
	void polarCenter(double x, double y, double z);
	bool operator<(atom b);
	bool operator<=(atom b);
	bool operator>=(atom b);
	bool operator>(atom b);
	double x;
	double y;
	double z;
	double theta;
	double psi;
	double r;
	double* scores = nullptr;
	int types = 0;
	~atom();
};

double euclideanDistance(atom a, atom b);
double angle3atoms(atom a, atom b, atom c);
double* normalVectorPlane(atom a, atom b, atom c);

double dihedralAngle(double* normalVector, atom b, atom a);


bool radiusCompare(atom a, atom b);

class structure {
public:
	structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool sort = false);
	structure(std::string input, std::vector<std::string> elements,std::vector<int> composition);
	structure(double x, double y, double z, int n, bool sort = false);//random generation in range
	structure(double x, double y, double z, std::vector<int> composition, std::vector<std::string> elements, bool sort = false);//random generation in range
	structure(structure& original, double variationX, double variationY, double variationZ, double cellX, double cellY, double cellZ, bool sort = false);//produce similar structure;
	structure(structure& original, std::vector<double> variation);
	structure(structure& original, std::vector<double> variation, double dx, double dy, double dz);
	structure(structure& original, double percentVariationX, double percentVariationY, double percentVariationZ, bool sort = false);
	bool operator<(const structure& other);
	void radiusSort();
	void rotateX(double degrees);
	void rotateY(double degrees);
	void rotateZ(double degrees);
	void generateCoordinationNumbers(int percent);
	void print();
	void makeZmatrix();
	void printZmatrix();
	void writeToObj(std::string outputName);
	void writeToXyz(std::string outputName);
	std::vector<std::vector<atom>> set;//first vector is atom types, second is atoms;
	std::vector<std::vector<int>> coordinationNumbers;
	std::vector<int> composition;
	std::vector<std::string> elements;
	double energy;
	double x, y, z;//geometric center
	bool scored = false;
	double** zmatrix = nullptr;
	std::string* zmatrixElements = nullptr;
	int zmatrixsize = 0;//equal to length of number of atoms, only set up for use in zmatrix
	~structure();
	//symmetries from group theory that do not need to be set
	int nsym = -1;
	int hsym = -1;


	//iterator
	struct Iterator
	{
		using iterator_category = std::forward_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using value_type = atom;
		using pointer = atom*;
		using reference = atom&;

		Iterator(pointer ptr , int indexI, int indexJ,std::vector<std::vector<atom>>& set);

		reference operator*() const;
		pointer operator->();
		Iterator& operator++();
		Iterator operator++(int);
		bool operator== (const Iterator& a);
		bool operator!= (const Iterator& a);

		
		pointer m_ptr;
		int indexI;
		int indexJ;
		std::vector<std::vector<atom>> set;
	};

	Iterator begin();
	Iterator end();
};

const double boltzman = 8.6173303e-5;//electron volts/kelvin
const double electronVolt = 1.602e-19;//joule

int covalentRadii(std::string element, int bond);
int vanDerWaalsRadii(std::string element);
double atomicRadius(std::string atomName);
bool bondingRequirement(structure* s, int percent, bool numBonds = 0);
double picometersToAngstrom(int picometers);
bool radialCriteria(structure& s, int percent);
bool covalentCriteria(structure* s, int percent);
bool bondingRequirement(structure* s, int percent, bool numBonds);
double MaximizeEuclideanDistance(atom a, atom b,double x, double y, double z);
double MinimizeEuclideanDistance(atom a, atom b, double x, double y, double z);
int radialCriteriaTrapping(structure& s, int percent, double x, double y, double z,int sensitivity = 5);
std::string atomicNumber(int num);




structure* seedFromGroupTheory(std::vector<int> composition, std::vector<std::string> elements,double rmax);
structure* seedWithSymmetry(int nsym, bool hsym, std::vector<int> composition, std::vector<std::string> elements, double rmax);
structure* proceduralSeedFromGroupTheory(int seed, std::vector<int> composition, std::vector<std::string> elements, double rmax);
structure* randomStructureFromGroupTheory(std::vector<int> composition, std::vector<std::string> elements, double rmax, int n, bool h, double phimax);
std::vector<std::pair<int, bool>> possibleSymmetry(std::vector<int> composition, std::vector<std::string> elements, double rmax, double rcp, int bondnum, int maxRCIterations = 10000);
structure* seedFromPossibleSymmetries(std::vector<std::pair<int, bool>> symmetry, std::vector<int> composition, std::vector<std::string> elements, double rmax, int index = -1);

std::vector<std::pair<int, bool>> possiblePartialSymmetry(std::vector<int> composition, std::vector<std::string> elements,double nMultiplier);
structure* getSeedFromSymmetryMode(char mode, std::vector<int> composition, std::vector<std::string> elements, double rmax,int call = 0, std::vector<std::pair<int, bool>> symmetries = { std::pair<int,bool>(1,0) });


std::vector<atom> partialSymmetrySet(int composition, int n, bool h, double rmax, int holes, bool alignHoles = true);



std::vector<bool*> enumeration(int size, int totalValue);

std::vector<bool*> alignedEnumeration(int size, int totalValue);

class partialSymmetryDescriptor {
	//a partial symmetry descriptor for a single element within a molecule
public:
	partialSymmetryDescriptor(int comp, int n, bool h, int holes, bool* holeLocations);
	int composition;
	int n;
	bool h;
	int holes;
	bool* holeLocations;
	~partialSymmetryDescriptor();
	void print();
};
std::vector<atom> partialSymmetrySet(partialSymmetryDescriptor descriptor, double rmax);


std::vector<structure*> seedsWithPartialSymmetry(bool extraHoles,bool alignHoles,bool variantSymmetry, double nMultiplier, std::vector<int> composition, std::vector<std::string> elements, double rmax, double radialCriteriaPercent, int bondRequirement, int criteriaIterations);

std::vector<structure*> getSeeds(std::vector<int> composition, std::vector<std::string> elements, double radius, double radialCriteriaPercent, int bondRequirement, int criteriaIterations, bool partialSymmetry = false, bool extraHoles = false, bool alignHoles = true, bool variantSymmetry = false, double nMultiplier = 1);

void writeToXyz(std::vector<structure*> structures,std::string outputName);

std::vector<structure*> structuresFromXYZ(std::string fileName);


bool energyCompare(structure* a, structure* b);
#endif

