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
#include <string>

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

bool radiusCompare(atom a, atom b);

class structure {
public:
	structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool sort = false);
	structure(std::string input, std::vector<std::string> elements,std::vector<int> composition);
	structure(double x, double y, double z, int n, bool sort = false);//random generation in range
	structure(double x, double y, double z, std::vector<int> composition, std::vector<std::string> elements = {"B"}, bool sort = false);//random generation in range
	structure(structure& original, double variationX, double variationY, double variationZ, double cellX, double cellY, double cellZ, bool sort = false);//produce similar structure;
	structure(structure& original, std::vector<double> variation);
	structure(structure& original, double percentVariationX, double percentVariationY, double percentVariationZ, bool sort = false);
	void radiusSort();
	void rotateX(double degrees);
	void rotateY(double degrees);
	void rotateZ(double degrees);
	void generateCoordinationNumbers(int percent);
	void print();
	std::vector<std::vector<atom>> set;//first vector is atom types, second is atoms;
	std::vector<std::vector<int>> coordinationNumbers;
	std::vector<int> composition;
	std::vector<std::string> elements;
	double energy;
	double x, y, z;//geometric center
	bool scored = false;



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
double picometersToAngstrom(int picometers);
bool radialCriteria(structure& s, int percent);
std::string atomicNumber(int num);


#endif
