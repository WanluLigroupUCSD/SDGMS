#ifndef chemH
#define chemH

#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <chrono>
#include <set>

class atom {
public:
	atom(double x, double y, double z);
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
	int types = 1;
	~atom();
};

double euclideanDistance(atom a, atom b);

bool radiusCompare(atom a, atom b);

class structure {
public:
	structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool sort = false);
	structure(double x, double y, double z, int n, bool sort = false);//random generation in range
	structure(double x, double y, double z, std::vector<int> composition, std::vector<std::string> elements = {"B"}, bool sort = false);//random generation in range
	structure(structure original, double variation, bool sort = false);//produce similar structure;
	structure(structure original, double percentVariationX, double percentVariationY, double percentVariationZ, bool sort = false);
	structure(structure original, double percentVariationX, double percentVariationY, double percentVariationZ, double rangeX, double rangeY, double rangeZ, bool sort = false);
	void radiusSort();
	void rotateX(double degrees);
	void rotateY(double degrees);
	void rotateZ(double degrees);
	void generateCoordinationNumbers(int percent);
	std::vector<std::vector<atom>> set;//first vector is atom types, second is atoms;
	std::vector<std::vector<int>> coordinationNumbers;
	std::vector<int> composition;
	std::vector<std::string> elements;
	double x, y, z;//geometric center
};

double boltzman = 8.6173303e-5;//electron volts/kelvin
double electronVolt = 1.602e-19;//joule

int covalentRadii(std::string element, int bond);
int vanDerWaalsRadii(std::string element);
bool radialCriteria(structure s, int percent);


#endif
