#ifndef dHeur
#define dHeur

#ifndef DBL_MAX
#define DBL_MAX 1.79769e+308;
#endif
#include "chemistry.hpp"
//no rotation algorithms
double distanceFast(structure A, structure B);
double distanceSlow(structure A, structure B);
double distanceAbsolute(structure A, structure B);

//rotational algorithms
double distanceRotate(structure A, structure B, double angle);
void scoreAtoms(structure &A);
double RID(structure& A, structure& B);//rotationally indendent distance
double RIDtoAngstroms(structure& A, double angstrom);
std::pair<double, double> RIDtoAngstroms(structure& A, double angstrom, int trials);

//determine cutoff value for structural similarity heuristic
double threshold(std::function<double(structure&, structure&)> func, std::vector<int> composition,std::vector<std::string> elements, double rangeX, double rangeY, double rangeZ, double similar, double moveX, double moveY, double moveZ, double percent = 0);
double branchAndBound(std::vector<std::vector<double>> input);
double exactDistance(std::vector<std::vector<atom>> setA, std::vector<std::vector<atom>> setB);
double threeAtomDistance(structure& a, structure& b);
#endif