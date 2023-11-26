#ifndef dHeur
#define dHeur



#include "chemistry.hpp"
//no rotation algorithms
double distanceFast(structure A, structure B);
double distanceSlow(structure A, structure B);
double distanceAbsolute(structure A, structure B);

//rotational algorithms
double distanceRotate(structure A, structure B, double angle);
void scoreAtoms(structure &A);
double RID(structure& A, structure& B);//rotationally indendent distance

//determine cutoff value for structural similarity heuristic
double threshold(std::function<double(structure&, structure&)> func, std::vector<int> composition, double rangeX, double rangeY, double rangeZ, double similar, int percent= 0 );
#endif