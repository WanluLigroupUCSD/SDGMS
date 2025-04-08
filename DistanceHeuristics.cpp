#include "DistanceHeuristics.hpp"


double distanceFast(structure A, structure B)
{
    double DistanceA = 0.0;
    double DistanceB = 0.0;
    int sets = A.set.size();
    for (int a = 0; a < sets; a++)
    {
        //for each type of atom
        int j = 0;
        int n = A.set[a].size();
        int m = B.set[a].size();
        for (int i = 0; i < n; i++)
        {
            double min = DBL_MAX;
            bool cont = true;
            //std::cout << "j: " << j << " m: " << m << " cont: " << cont << " i: " << i << std::endl;
            while (cont)
            {
                if(j < m)
                {
                    double dist = euclideanDistance(A.set[a][i], B.set[a][j]);
                    //std::cout << "current distance A to B: " << dist << std::endl;
                    if (dist < min)
                    {
                        min = dist;
                        j++;
                    }
                    else {
                        cont = false;
                        /*
                        As A.set[a][i] < A.set[a][i+1] and B.set[a][j] is closest to A.set[a][i]
                        for A.set[a][i+ 1] any B.set[a][j-n] is further than B.set[a][j] (by radius)
                        but maybe B.set[a][j]  is closer than B.set[a][j+1] (by radius) so we will check it by decrementing j
                        */
                        //std::cout << "decreasing j from: " << j << std::endl;
                        j--;
                    
                    }
                }
                else {
                    //set A stretches further than set B, so compare all atoms from now on to the last of set B
                    min = euclideanDistance(A.set[a][i], B.set[a][m-1]);
                    cont = false;
                }
            }
            DistanceA += min;
            //std::cout << "distanceA: " << DistanceA << " min: " << min << std::endl;
        }
        
        //calculate for structure B
        j = 0;
        for (int i = 0; i < m; i++)
        {
            double min = DBL_MAX;
            bool cont = true;
            while (cont)
            {
                if (j < n)
                {
                    double dist = euclideanDistance(B.set[a][i], A.set[a][j]);
                    //std::cout << "current distance B to A: " << dist << std::endl;
                    if (dist < min)
                    {
                        min = dist;
                        j++;
                    }
                    else {
                        cont = false;
                        j--;
                    }
                }
                else {
                    //set B stretches further than set A, so compare all atoms from now on to the last of set A
                    min = euclideanDistance(A.set[a][i], A.set[a][n - 1]);
                    cont = false;
                }
            }
            DistanceB += min;
        }
    }
    //std::cout << DistanceA << " " << DistanceB << std::endl;
    return DistanceA + DistanceB;
    
    
}

double distanceSlow(structure A, structure B)
{
    double DistanceA = 0.0;
    double DistanceB = 0.0;
    int sets = A.set.size();
    for (int a = 0; a < sets; a++)
    {
        //for each type of atom
        int n = A.set[a].size();
        int m = B.set[a].size();
        for (int i = 0; i < n; i++)
        {
            double min = DBL_MAX;
            for (int j = 0; j < m; j++)
            {
                double dist = euclideanDistance(A.set[a][i], B.set[a][j]);
                if (dist < min)min = dist;
            }
            DistanceA += min;
        }

        //calculate for structure B
        for (int i = 0; i < m; i++)
        {
            double min = DBL_MAX;
            for (int j = 0; j < n; j++)
            {
                double dist = euclideanDistance(B.set[a][i], A.set[a][j]);
                if (dist < min) min = dist;
            }
            DistanceB += min;
        }
    }
    //std::cout << DistanceA << " " << DistanceB << std::endl;
    return DistanceA + DistanceB;


}

double distanceAbsoluteRecursion(structure& A, structure& B, bool* lA, bool* lB, int type, int pairs)
{
    //if (lA[0] && lB[0]) std::cout << "what" << std::endl;
    pairs--;
    if (pairs == -1) return 0;
    //duplicate the lists
    int listSize = A.set[type].size();
    double min = DBL_MAX;
    for (int i = 0; i < listSize; i++)
    {
        //std::cout << "i: " << i << " pairs: " << pairs << std::endl;
        
        
        if (lA[i])
        {
            
            for (int j = 0; j < listSize; j++)
            {
                //std::cout << "j: " << j << " pairs: " << pairs << std::endl;

                //atom B[j]
                if (lB[j])
                {
                    //atom lB will be calculated so set it to false now

                    //lB[j] = false;
                    //std::cout << "calculate " << i << " " << j << " pairs: " << pairs << std::endl;
                    bool* lA2 = new bool[listSize];
                    bool* lB2 = new bool[listSize];

                    for (int z = 0; z < listSize; z++)
                    {
                        lA2[z] = lA[z];
                        lB2[z] = lB[z];
                    }
                    lA2[i] = false;
                    lB2[j] = false;
                    //std::cout << "determine " << i << " " << j << " pairs: " << pairs << std::endl;
                    double dist = distanceAbsoluteRecursion(A, B, lA2, lB2, type,pairs);
                    
                    dist += euclideanDistance(A.set[type][i], B.set[type][j]);
                    //std::cout << "dist: " << dist << " pairs: " << pairs << std::endl;
                    if (dist < min) min = dist;

                    delete[] lA2;
                    delete[] lB2;
                }
            }
        }
    }
    //std::cout << pairs << " " << min << std::endl;
    return min;
}
double distanceAbsolute(structure A, structure B)
{
    int types = A.set.size();
    double total = 0;
    for (int i = 0; i < types; i++)
    {
        int listSize = A.set[i].size();
        bool* lA = new bool[listSize];
        bool* lB = new bool[listSize];
        for (int j = 0; j < listSize; j++)
        {
            
            //initialize all atoms to be available
            lA[j] = true;
            lB[j] = true;

            
        }
        total += distanceAbsoluteRecursion(A, B, lA, lB, i, listSize);


        delete[] lA;
        delete[] lB;
    }
    return total;
}


double distanceRotate(structure A, structure B, double angle)
{
    /*
    TODO
    */
    return 0;
}
void scoreAtoms(structure &A)
{
    int types = A.set.size();
    for (std::vector<std::vector<atom>>::iterator it = A.set.begin(); it != A.set.end(); it++)
    {
        for (std::vector<atom>::iterator a = it->begin(); a != it->end(); a++)
        {
            //for each atom a
            //initialize scores
            a->scores = new double[types];
            a->types = types;
            for (int t = 0; t < types; t++)
            {
                double score = 0;
                for (std::vector<atom>::iterator b = A.set[t].begin(); b != A.set[t].end(); b++)
                {
                    double x = (a->x - b->x);
                    double y = (a->y - b->y);
                    double z = (a->z - b->z);
                    score += std::sqrt(x * x + y * y + z * z);
                    //score += euclideanDistance(*a, *b);
                }
                a->scores[t] = score;
            }
        }
    }
    A.scored = true;
}
double RID(structure& A, structure& B)
{
    /*
    TODO
    */
    int types = A.set.size();
    //structure A
    if(!A.scored) scoreAtoms(A);
    if(!B.scored) scoreAtoms(B);
    double totalA = 0;
    for (int t = 0; t < types; t++)
    {
        for (std::vector<atom>::iterator a = A.set[t].begin(); a != A.set[t].end(); a++)
        {
            //for each atom a in structure A
            //which atom in b is it closest too?
            double min = DBL_MAX;
            for (std::vector<atom>::iterator b = B.set[t].begin(); b != B.set[t].end(); b++)
            {
                double stotal = 0;
                //find distance
                for (int i = 0; i < types; i++)
                {
                    double dist = a->scores[i] - b->scores[i];
                    if (dist > 0) stotal += dist;
                    else stotal -= dist;
                }
                if (stotal < min) min = stotal;
            }
            totalA += min;

        }
    }
    //repeat for B
    double totalB = 0;
    for (int t = 0; t < types; t++)
    {
        for (std::vector<atom>::iterator b = B.set[t].begin(); b != B.set[t].end(); b++)
        {
            //for each atom b in structure B
            //which atom in a is it closest too?
            double min = DBL_MAX;
            for (std::vector<atom>::iterator a = A.set[t].begin(); a != A.set[t].end(); a++)
            {
                double stotal = 0;
                //find distance
                for (int i = 0; i < types; i++)
                {
                    double dist = b->scores[i] - a->scores[i];
                    if (dist > 0) stotal += dist;
                    else stotal -= dist;
                }
                if (stotal < min) min = stotal;
            }
            totalB += min;

        }
    }

    
    return totalA + totalB;
}

double RIDtoAngstroms(structure& A, double angstrom)
{
    int natoms = 0;
    for (int s = 0; s < A.set.size(); s++)  natoms += A.set[s].size();
    std::vector<double> movement = {};
    for (int s = 0; s < A.set.size(); s++) for (int a = 0; a < A.set[s].size(); a++) {
        movement.push_back(1 - 2*(double)rand() / (double)RAND_MAX);
        movement.push_back(1 - 2*(double)rand() / (double)RAND_MAX);
        movement.push_back(1 - 2*(double)rand() / (double)RAND_MAX);
    }

    double magnitude = 0;
    for (int m = 0; m < movement.size(); m++) magnitude += movement[m] * movement[m];
    for (int m = 0; m < movement.size(); m++) movement[m] = angstrom* movement[m] / sqrt(magnitude);

    //magnitude = 0;
    //for (int m = 0; m < movement.size(); m++) magnitude += movement[m] * movement[m];
    //std::cout << "new magnitiude: " << magnitude << std::endl;

    //for (std::vector<double>::iterator it = movement.begin(); it != movement.end(); it++) normalization += *it * *it;
    //normalization = sqrt(normalization);

    structure B(A, movement);
    return RID(A, B);
}

std::pair<double,double> RIDtoAngstroms(structure& A, double angstrom, int trials)
{
    //return the average and stdev RID at an angstrom value for trials
    std::vector<double> values = {};
    for (int t = 0; t < trials; t++) {
        std::vector<double> movement = {};
        for (int s = 0; s < A.set.size(); s++) for (int a = 0; a < A.set[s].size(); a++) {
            movement.push_back(1-2*(double)rand() / (double)RAND_MAX);
            movement.push_back(1-2*(double)rand() / (double)RAND_MAX);
            movement.push_back(1-2*(double)rand() / (double)RAND_MAX);
        }

        double magnitude = 0;
        for (int m = 0; m < movement.size(); m++) magnitude += movement[m] * movement[m];
        for (int m = 0; m < movement.size(); m++) movement[m] = angstrom * movement[m] / sqrt(magnitude);

        
        structure B(A, movement);
        values.push_back( RID(A, B));
    }
    
    double average = 0;
    for (int i = 0; i < values.size(); i++) {
        average += values[i];
    }
    average = average / (double)values.size();
    double stdev = 0;
    for (int i = 0; i < values.size(); i++) {
        stdev += (values[i] - average) * (values[i] - average);
    }
    stdev = sqrt(stdev / (double)values.size());
    return std::pair<double, double>(average, stdev);
}

double threshold(std::function<double(structure&, structure&)> heuristic, std::vector<int> composition,std::vector<std::string> elements, double rangeX, double rangeY, double rangeZ, double similar,double moveX, double moveY, double moveZ, double percent)
{
    //we want to differentiate similar structures from structures that are likely to be generated from movements.

    //composition is a list of integers. each integer represents a type of atom and its value represents the number of atoms
    //similar is the variation allowed on each point where it can be considered a similar structure (percentage of the range)
    //different is the variation on each point where it would no longer be considered a similar structure (percentage of the range)
    //percent is the percentage allowed of different structures to be taken
    //rangeX,y,z are the bounds of the structure

    bool bugging = false;

    std::vector<double> heuristicSimilar;
    std::vector<double> heuristicDifferent;

    int trials = 10000;

    for (int i = 0; i < trials; i++)
    {
        structure A = structure(rangeX, rangeY, rangeZ, composition,elements);
        std::vector<double> movementDifferent;
        std::vector<double> movementSimilar;
        int tParameters = 0;
        for (int i = 0; i < A.elements.size(); i++)
        {
            for (int j = 0; j < A.set[i].size(); j++)
            {
                //should this be randomized between -1 to 1?
                movementSimilar.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * moveX * percent / 100);
                movementDifferent.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1)* moveX);
                movementSimilar.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * moveY * percent / 100);
                movementDifferent.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * moveY);
                movementSimilar.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * moveZ * percent / 100);
                movementDifferent.push_back((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * moveZ);

            }
        }
        if ((i+1) % 1000 == 0) std::cout << i / 100 << "% of trials generated" << std::endl;
        
        structure S = structure(A, movementSimilar);
        structure D = structure(A, movementDifferent);
        heuristicSimilar.push_back(heuristic(A, S));
        heuristicDifferent.push_back(heuristic(A, D));
    }

    std::sort(heuristicSimilar.begin(), heuristicSimilar.end());
    std::sort(heuristicDifferent.begin(), heuristicDifferent.end());

    /*print distribution
    std::cout << "similar:" << std::endl;
    for (int i = 0; i < 100; i += 10)
    {
        std::cout << i << "%: " << heuristicSimilar[((trials / 100) * i) ] << std::endl;
    }
    std::cout << 100 << "%: " << heuristicSimilar[trials - 1] << std::endl;

    std::cout << "different:" << std::endl;
    for (int i = 0; i < 100; i += 10)
    {
        std::cout << i << "%: " << heuristicDifferent[((trials / 100) * i )] << std::endl;
    }
    std::cout << 100 << "%: " << heuristicDifferent[trials - 1] << std::endl;
    */

    for (int i = trials - 1; i >= 0; i--)
    {
        if (heuristicDifferent[(percent/100)*trials] > heuristicSimilar[i])
        {
            std::cout << "percent of similar lost: " << (trials - i - 1) / (trials - 1) << std::endl;
            return heuristicSimilar[i];
        }
    }
    
    
        
    return 0;

    
    //this algorithm will create a structure, generate a structure considered similar, and a structure considered different, and then calculate the heuristic differences.
    //with a distribution of these heuristics, a value will be taken
}

// Compute rotation matrix to align vector a to vector b
std::vector<std::vector<double>> rotationMatrixToAlign(atom& a, atom& b) {
    std::vector<double> v = a.cross(b).vec3();
    double c = a.dot(b);
    double k = 1 / (1 + c);

    std::vector<std::vector<double>> R = {
        {v[0] * v[0] * k + c, v[0] * v[1] * k - v[2], v[0] * v[2] * k + v[1]},
        {v[1] * v[0] * k + v[2], v[1] * v[1] * k + c, v[1] * v[2] * k - v[0]},
        {v[2] * v[0] * k - v[1], v[2] * v[1] * k + v[0], v[2] * v[2] * k + c}
    };

    return R;
}

// Apply rotation matrix to a vector
atom applyRotation(std::vector<double> point, const std::vector<std::vector<double>>& R) {
    return atom(
        R[0][0] * point[0] + R[0][1] * point[1] + R[0][2] * point[2],
        R[1][0] * point[0] + R[1][1] * point[1] + R[1][2] * point[2],
        R[2][0] * point[0] + R[2][1] * point[1] + R[2][2] * point[2]);
}




double branchAndBound(std::vector<std::vector<double>> input)
{
    //std::cout << "BAB RMS 1" << std::endl;
    if (input.size() == 1) return input[0][0];
    //std::cout << "starting new BNB" << std::endl;
    std::vector<std::vector<bool>> pattern = {};

    //compute the top row, lower and upper bounds
    std::vector<double> upperBounds = {};
    std::vector<double> lowerBounds = {};
    //std::cout << "BAB RMS 2" << std::endl;
    for (int i = 0; i < input.size(); i++)
    {
        //std::cout << "BAB RMS 3" << std::endl;
        double ub = input[0][i];
        double lb = input[0][i];
        for (int j = 1; j < input.size(); j++)
        {
            double max = 0;
            double min = DBL_MAX;
            for (int k = 0; k < input.size(); k++)
            {
                //std::cout << "BAB RMS 4" << std::endl;
                if (k != i)
                {
                    if (input[j][k] > max) max = input[j][k];
                    if (input[j][k] < min) min = input[j][k];
                }
                //std::cout << "BAB RMS 5" << std::endl;
            }
            ub += max;
            lb += min;
            //std::cout << "BAB RMS 6" << std::endl;
        }
        upperBounds.push_back(ub);
        lowerBounds.push_back(lb);
        //std::cout << "BAB RMS 7" << std::endl;
        //std::cout << "pushing upper and lower bounds: " << ub << "," << lb << std::endl;
    }
    double lowestUpperBound = DBL_MAX;
    //std::cout << "BAB RMS 8" << std::endl;
    for (int i = 0; i < upperBounds.size(); i++)
    {
        if (upperBounds[i] < lowestUpperBound) lowestUpperBound = upperBounds[i];
    }
    std::vector<bool> keep = {};
    for (int i = 0; i < lowerBounds.size(); i++)
    {
        if (lowerBounds[i] > lowestUpperBound)
        {
            keep.push_back(false);
        }
        else keep.push_back(true);
    }

    //std::cout << "BAB RMS 9" << std::endl;
    std::vector<double> smallestSolutions = {};
    for (int k = 0; k < keep.size(); k++)
    {
        //std::cout << "BAB RMS 10" << std::endl;
        if (keep[k]) {
            std::vector<std::vector<double>> smallerInput;

            for (int i = 1; i < input.size(); i++)
            {
                //skipping the first input
                std::vector<double> newVec = {};
                for (int j = 0; j < input[i].size(); j++)
                {
                    if (k != j) newVec.push_back(input[i][j]);
                }
                smallerInput.push_back(newVec);
            }
            //std::cout << "BAB RMS SMALLER SOLUTION" << std::endl;
            double smallerSolution = branchAndBound(smallerInput);
            //std::cout << "smaller solution: " << smallerSolution << " with " << input[0][k] << std::endl;
            smallestSolutions.push_back(smallerSolution + input[0][k]);
        }
        /*else {
            std::cout << "Smaller BNB Excluded by bounds" << std::endl;
        }*/
        //std::cout << "BAB RMS 11" << std::endl;
       
    }
    //std::cout << "BAB RMS 12" << std::endl;
    double minDistance = DBL_MAX;
    for (int k = 0; k < smallestSolutions.size(); k++)
    {
        if (minDistance > smallestSolutions[k]) minDistance = smallestSolutions[k];
    }
    //std::cout << "returning: " << minDistance << std::endl;
    //std::cout << "BAB RMS 13________________________________________________________-" << std::endl;
    return minDistance;
}
    

double exactDistance(std::vector<std::vector<atom>> setA, std::vector<std::vector<atom>> setB)
{
    //finds the exact distance (for the current rotation) between setA and setB
    //create a matrix computing the distance between each atom
    int natoms = 0;
    for (std::vector<std::vector<atom>>::iterator a = setA.begin(); a != setA.end(); a++) natoms += a->size();
    std::vector<std::vector<double>> distanceMatrix = {};
    for (int i = 0; i < natoms; i++) {
        distanceMatrix.push_back({});
        for (int j = 0; j < natoms; j++) {
            distanceMatrix[i].push_back(0);
        }
    }

    int it = 0;
    for (int e = 0; e < setA.size(); e++) for (int a = 0; a < setA[e].size(); a++)
    {
        int jt = 0;
        for (int eb = 0; eb < setB.size(); eb++) for (int b = 0; b < setB[eb].size(); b++)
        {
            distanceMatrix[it][jt] = euclideanDistance(setA[e][a], setB[eb][b]);
            jt++;
        }
        it++;
    }
    return branchAndBound(distanceMatrix);

}

double threeAtomDistance(structure& a, structure& b)
{
    double minDistance = DBL_MAX;
    //pick all combinations of atom a, b, c to compare
    for (int ae = 0; ae < a.set.size(); ae++) for (int aa = 0; aa < a.set[ae].size(); aa++)
    {
        std::cout << "loop RMS 1" << std::endl;
        //std::cout << "1 " << ae << " "<< aa << std::endl;
        //for each first atom in set a
        for (int aeb = 0; aeb < a.set.size(); aeb++) for (int ab = 0; ab < a.set[aeb].size(); ab++) if(aeb != ae || ab != aa)
        {
            std::cout << "loop RMS 2" << std::endl;
            //std::cout << "2" << std::endl;
            //for each second atom in a
            for (int aec = 0; aec < a.set.size(); aec++) for (int ac = 0; ac < a.set[aec].size(); ac++) if ((aec != ae || ac != aa) && (aeb != aec || ac != ab))
            {
                std::cout << "loop RMS 3" << std::endl;
                //std::cout << "3" << std::endl;

                //For each combination of 3 atoms in set a
                for (int ba = 0; ba < b.set[ae].size(); ba++)
                {
                    std::cout << "loop RMS 4" << std::endl;
                    //std::cout << "4" << std::endl;
                    //For each atom in b that can match atom a from a
                    for (int bb = 0; bb < b.set[aeb].size(); bb++) if(aeb != ae || ba != bb)
                    {
                        std::cout << "loop RMS 5" << std::endl;
                        //std::cout << "5" << std::endl;
                        //For each second atom in b
                        for (int bc = 0; bc < b.set[aec].size(); bc++) if ((aec != ae || bc != ba) && (aeb != aec || bc != bb) )
                        {
                            std::cout << "loop RMS 6" << std::endl;
                            //std::cout << "6" << std::endl;

                            //For each set of 3 atoms in a and b, rotate b to be in the same plane as a and compare distance
                            // Example points in set A and B
                            structure rotatedB = structure(b.set, b.elements);

                            // Select three points from each set
                            atom A1 = a.set[ae][aa];
                            atom A2 = a.set[aeb][ab];
                            atom A3 = a.set[aec][ac];
                            atom B1 = a.set[ae][ba];
                            atom B2 = a.set[aeb][bb];
                            atom B3 = a.set[aec][bc];
                            // Create normal vectors of the planes
                            atom normalA = (A2 - A1).cross(A3 - A1).normalize();
                            atom normalB = (B2 - B1).cross(B3 - B1).normalize();
                            std::cout << "Normal A: " << normalA.x << " " << normalA.y << " " << normalA.z << std::endl;
                            std::cout << "Normal B: " << normalB.x << " " << normalB.y << " " << normalB.z << std::endl;

                            // Compute rotation matrix
                            std::vector<std::vector<double>> R = rotationMatrixToAlign(normalB, normalA);
                            std::cout << "RMS?1" << std::endl;
                            // Rotate points in set B
                            std::vector<std::vector<atom>> rotatedSet = {};
                            for (std::vector<std::vector<atom>>::iterator e = b.set.begin(); e != b.set.end(); e++){
                                std::vector<atom> rset = {};
                                for (std::vector<atom>::iterator at = e->begin(); at != e->end(); at++) {
                                    rset.push_back(applyRotation(at->vec3(), R));
                                    
                                }
                                rotatedSet.push_back(rset);
                            }
                            std::cout << "RMS?2" << std::endl;
                            // Print rotated points in set B
                            
                            //now find the distance between the two sets.
                            double dist = exactDistance(a.set, rotatedSet);
                            std::cout << "RMS?3" << std::endl;
                            structure(rotatedSet, a.elements).print();
                            std::cout << "current distance from exact distance: " << dist << std::endl;
                            if (dist < minDistance)  minDistance = dist;
                            std::cout << "RMS?4" << std::endl;

                        }
                        std::cout << "outer loop RMS 5" << std::endl;
                    }
                    std::cout << "outer loop RMS 4" << std::endl;
                }
                std::cout << "outer loop RMS 3" << std::endl;

            }
            std::cout << "outer loop RMS 2" << std::endl;
        }
        std::cout << "outer loop RMS 1" << std::endl;
    }
    std::cout << "outer loop RMS 0" << std::endl;
    return minDistance;
}
