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

double threshold(std::function<double(structure&, structure&)> heuristic, std::vector<int> composition, double rangeX, double rangeY, double rangeZ, double similar, int percent)
{
    //composition is a list of integers. each integer represents a type of atom and its value represents the number of atoms
    //similar is the variation allowed on each point where it can be considered a similar structure (percentage of the range)
    //different is the variation on each point where it would no longer be considered a similar structure (percentage of the range)
    //percent is the percentage allowed of different structures to be taken
    //rangeX,y,z are the bounds of the structure

    std::vector<double> heuristicSimilar;
    std::vector<double> heuristicDifferent;

    int trials = 10000;

    for (int i = 0; i < trials; i++)
    {
        if ((i+1) % 1000 == 0) std::cout << i / 100 << "% of trials generated" << std::endl;
        structure A = structure(rangeX, rangeY, rangeZ, composition);
        structure S = structure(A, similar/100 * rangeX, similar/100 * rangeY, similar/100 * rangeZ);
        structure D = structure(A, similar/100 * rangeX, similar/100 * rangeY, similar/100 * rangeZ,rangeX,rangeY,rangeZ);

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


