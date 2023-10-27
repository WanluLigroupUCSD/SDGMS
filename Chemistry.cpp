#include "chemistry.h"
atom::atom(double x, double y, double z) :x(x), y(y), z(z)
{
    r = std::sqrt(x * x + y * y + z * z);
    psi = std::acos((x / std::sqrt(x * x + y * y)));
    theta = z / r;
}
atom::~atom()
{
    delete[] scores;
}

void atom::polarCenter(double x, double y, double z)
{
    double newX = this->x - x;
    double newY = this->y - y;
    double newZ = this->z - z;
    r = std::sqrt(newX * newX + newY * newY + newZ * newZ);
    psi = std::acos((newX / std::sqrt(newX * newX + newY * newY)));
    theta = newZ / r;
}

bool atom::operator<(atom b) { return this->r < b.r; };
bool atom::operator<=(atom b) { return this->r <= b.r; };
bool atom::operator>=(atom b) { return this->r >= b.r; };
bool atom::operator>(atom b) { return this->r > b.r; };


bool radiusCompare(atom a, atom b) { return a.r < b.r; };

double euclideanDistance(atom a, atom b) {
	double x = (a.x - b.x);
	double y = (a.y - b.y);
	double z = (a.z - b.z);
	double d = std::sqrt(x * x + y * y + z * z);
	return d;
}



structure::structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool sort) :set(set), elements(elements)
{
    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int n = 0;
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            n++;
            totalX += itt->x;
            totalY += itt->y;
            totalZ += itt->z;
        }
    }
    x = totalX / (double)n;
    y = totalY / (double)n;
    z = totalZ / (double)n;

    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        composition.push_back(it->size());
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            itt->polarCenter(x, y, z);
        }
    }
    if (sort)this->radiusSort();
}
structure::structure(double x, double y, double z, int n, bool sort)
{
    std::vector<atom> atoms;
    for (int i = 0; i < n; i++)
    {

        double newX = x * (rand() % RAND_MAX) / RAND_MAX;
        double newY = y * (rand() % RAND_MAX) / RAND_MAX;
        double newZ = z * (rand() % RAND_MAX) / RAND_MAX;
        atoms.push_back(atom(newX, newY, newZ));
        //std::cout << "new atom at " << newX << "," << newY << "," << newZ << std::endl;
    }
    set.push_back(atoms);


    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int m = 0;
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            m++;
            totalX += itt->x;
            totalY += itt->y;
            totalZ += itt->z;
        }
    }
    this->x = totalX / (double)m;
    this->y = totalY / (double)m;
    this->z = totalZ / (double)m;

    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        composition.push_back(it->size());
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            itt->polarCenter(this->x, this->y, this->z);
        }
    }
    if (sort)this->radiusSort();
    /*
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            std::cout << "sorted: "<< itt->r << ":" << itt->x << "," << itt->y << "," << itt->z << std::endl;
        }
    }
    */
}
structure::structure(double x, double y, double z, std::vector<int> composition, std::vector<std::string> elements, bool sort) :elements(elements)
{

    for (std::vector<int>::iterator it = composition.begin(); it != composition.end(); it++)
    {
        std::vector<atom> atoms;
        for (int i = 0; i < *it; i++)
        {

            double newX = x * (rand() % RAND_MAX) / RAND_MAX;
            double newY = y * (rand() % RAND_MAX) / RAND_MAX;
            double newZ = z * (rand() % RAND_MAX) / RAND_MAX;
            atoms.push_back(atom(newX, newY, newZ));
            //std::cout << "new atom at " << newX << "," << newY << "," << newZ << std::endl;
        }
        set.push_back(atoms);
    }



    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int m = 0;
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            m++;
            totalX += itt->x;
            totalY += itt->y;
            totalZ += itt->z;
        }
    }
    this->x = totalX / (double)m;
    this->y = totalY / (double)m;
    this->z = totalZ / (double)m;

    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        composition.push_back(it->size());
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            itt->polarCenter(this->x, this->y, this->z);
        }
    }
    if (sort)this->radiusSort();
}
structure::structure(structure original, double variation, bool sort) :elements(original.elements)
{
    for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
    {
        std::vector<atom> atoms;
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double varX = 2 * variation * rand() / RAND_MAX - 1;
            double varY = 2 * variation * rand() / RAND_MAX - 1;
            double varZ = 2 * variation * rand() / RAND_MAX - 1;
            atom at = atom(itt->x + varX, itt->y + varY, itt->z + varZ);
            //std::cout << "new atom at " << at.x << "," << at.y << "," << at.z << std::endl;
            atoms.push_back(at);
        }
        set.push_back(atoms);
    }

    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int n = 0;
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            n++;
            totalX += itt->x;
            totalY += itt->y;
            totalZ += itt->z;
        }
    }
    x = totalX / (double)n;
    y = totalY / (double)n;
    z = totalZ / (double)n;

    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        composition.push_back(it->size());
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            itt->polarCenter(x, y, z);
        }
    }
    if (sort)this->radiusSort();
}
structure::structure(structure original, double variationX, double variationY, double variationZ, bool sort) :elements(original.elements)
{
    //makes variants of at most variation x,y,z
    for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
    {
        std::vector<atom> atoms;
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double randX = ((2 * rand()) / RAND_MAX) - 1;
            double randY = ((2 * rand()) / RAND_MAX) - 1;
            double randZ = ((2 * rand()) / RAND_MAX) - 1;
            atom at = atom(itt->x + variationX * randX, itt->y + variationY * randY, itt->z + variationZ * randZ);
            atoms.push_back(at);
        }
        set.push_back(atoms);
    }

    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int n = 0;
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            n++;
            totalX += itt->x;
            totalY += itt->y;
            totalZ += itt->z;
        }
    }
    x = totalX / (double)n;
    y = totalY / (double)n;
    z = totalZ / (double)n;

    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        composition.push_back(it->size());
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            itt->polarCenter(x, y, z);
        }
    }
    if (sort)this->radiusSort();
}
structure::structure(structure original, double variationX, double variationY, double variationZ, double rangeX, double rangeY, double rangeZ, bool sort) :elements(original.elements)
{
    //makes variants of structure original where the minimum variation is varX,y& z and maximum is range
    for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
    {
        std::vector<atom> atoms;
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double randX = ((2 * rand()) / RAND_MAX) - 1;
            double randY = ((2 * rand()) / RAND_MAX) - 1;
            double randZ = ((2 * rand()) / RAND_MAX) - 1;
            atom at = atom(itt->x + variationX + (rangeX - variationX) * randX, itt->y + variationY + (rangeY - variationY) * randY, itt->z + variationZ + (rangeZ - variationZ) * randZ);
            atoms.push_back(at);
        }
        set.push_back(atoms);
    }

    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int n = 0;
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            n++;
            totalX += itt->x;
            totalY += itt->y;
            totalZ += itt->z;
        }
    }
    x = totalX / (double)n;
    y = totalY / (double)n;
    z = totalZ / (double)n;

    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        composition.push_back(it->size());
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            itt->polarCenter(x, y, z);
        }
    }
    if (sort)this->radiusSort();
}

void structure::radiusSort()
{
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        std::sort(it->begin(), it->end(), radiusCompare);
    }
}

void structure::rotateX(double x)
{
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double oldX = itt->x;
            double oldY = itt->y;
            double oldZ = itt->z;

            double newX = oldX;
            double newY = oldY * std::cos(x) - oldZ * std::sin(x);
            double newZ = oldY * std::sin(x) + oldZ * std::cos(x);

            itt->x = newX;
            itt->y = newY;
            itt->z = newZ;
        }
    }
}

void structure::rotateY(double y)
{
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double oldX = itt->x;
            double oldY = itt->y;
            double oldZ = itt->z;

            double newX = oldX * std::cos(y) + oldZ * std::sin(y);
            double newY = oldY;
            double newZ = -oldX * std::sin(y) + oldZ * std::cos(y);

            itt->x = newX;
            itt->y = newY;
            itt->z = newZ;
        }
    }
}

void structure::rotateZ(double z)
{
    for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
    {
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double oldX = itt->x;
            double oldY = itt->y;
            double oldZ = itt->z;

            double newX = oldX * std::cos(z) - oldY * std::sin(z);
            double newY = oldX * std::sin(z) + oldY * std::cos(z);
            double newZ = oldZ;

            itt->x = newX;
            itt->y = newY;
            itt->z = newZ;
        }
    }
}

void structure::generateCoordinationNumbers(int percent)
{
	for (int i = 0; i < this->set.size(); i++)
	{
		std::vector<int> celement;
		for (int a = 0; a < this->set[i].size(); a++)
		{
			//for each atom b in s
			for (int j = 0; j < this->set.size(); j++)
			{
				int bonds = 0;
				for (int b = 0; b < this->set[j].size(); b++)
				{
					bool accept = false;
					if (i == j && b == a) {
					}
					else {
						//do not perform on the same atom
						double dist = euclideanDistance(this->set[i][a], this->set[j][b]);
						
						double bl1 = covalentRadii(this->elements[i], 1);
						double bl1r = bl1 * ((double)percent / (double)100);
						if (dist > bl1 - bl1r && dist < bl1 + bl1r) accept = true;
						else {
							double bl2 = covalentRadii(this->elements[i], 2);
							if (bl2 != 0)
							{
								double bl2r = bl2 * ((double)percent / (double)100);
								if (dist > bl2 - bl2r && dist < bl2 + bl2r) accept = true;
								else {
									double bl3 = covalentRadii(this->elements[i], 3);
									if (bl3 != 0)
									{
										double bl3r = bl3 * ((double)percent / (double)100);
										if (dist > bl3 - bl3r && dist < bl3 + bl3r) accept = true;
									}
								}
							}




						}
					}
					if (accept) bonds++;



				}
				celement.push_back(bonds);
			}

		}
		coordinationNumbers.push_back(celement);
	}
}


//chemistry
int covalentRadii(std::string element, int bond)
{
	if (element == "H") {
		if (bond == 1) return 32;
	}
	else if (element == "He") {
		if (bond == 1) return 46;
	}
	else if (element == "Li") {
		if (bond == 1) return 133;
		else if (bond == 2) return 124;
	}
	else if (element == "Be") {
		if (bond == 1) return 102;
		else if (bond == 2) return 90;
		else if (bond == 3) return 85;
	}
	else if (element == "B") {
		if (bond == 1) return 85;
		else if (bond == 2) return 78;
		else if (bond == 3) return 73;
	}
	else if (element == "C") {
		if (bond == 1) return 75;
		else if (bond == 2) return 67;
		else if (bond == 3) return 60;
	}
	else if (element == "N") {
		if (bond == 1) return 71;
		else if (bond == 2) return 60;
		else if (bond == 3) return 54;
	}
	else if (element == "O") {
		if (bond == 1) return 63;
		else if (bond == 2) return 57;
		else if (bond == 3) return 53;
	}
	else if (element == "F") {
		if (bond == 1) return 64;
		else if (bond == 2) return 59;
		else if (bond == 3) return 53;
	}
	else if (element == "Ne") {
		if (bond == 1) return 67;
		else if (bond == 2) return 96;
	}
	else if (element == "Na") {
		if (bond == 1) return 155;
		else if (bond == 2) return 160;
	}
	else if (element == "Mg") {
		if (bond == 1) return 139;
		else if (bond == 2) return 132;
		else if (bond == 3) return 127;
	}
	else if (element == "Al") {
		if (bond == 1) return 126;
		else if (bond == 2) return 113;
		else if (bond == 3) return 111;
	}
	else if (element == "Si") {
		if (bond == 1) return 116;
		else if (bond == 2) return 107;
		else if (bond == 3) return 102;
	}
	else if (element == "P") {
		if (bond == 1) return 111;
		else if (bond == 2) return 102;
		else if (bond == 3) return 94;
	}
	else if (element == "S") {
		if (bond == 1) return 103;
		else if (bond == 2) return 94;
		else if (bond == 3) return 95;
	}
	else if (element == "Cl") {
		if (bond == 1) return 99;
		else if (bond == 2) return 95;
		else if (bond == 3) return 93;
	}
	else if (element == "Ar") {
		if (bond == 1) return 96;
		else if (bond == 2) return 107;
		else if (bond == 3) return 96;
	}
	else if (element == "K") {
		if (bond == 1) return 196;
		else if (bond == 2) return 193;
	}
	else if (element == "Ca") {
		if (bond == 1) return 171;
		else if (bond == 2) return 147;
		else if (bond == 3) return 133;
	}
	else if (element == "Sc") {
		if (bond == 1) return 148;
		else if (bond == 2) return 116;
		else if (bond == 3) return 114;
	}
	else if (element == "Ti") {
		if (bond == 1) return 136;
		else if (bond == 2) return 117;
		else if (bond == 3) return 108;
	}
	else if (element == "V") {
		if (bond == 1) return 134;
		else if (bond == 2) return 112;
		else if (bond == 3) return 106;
	}
	else if (element == "Cr") {
		if (bond == 1) return 122;
		else if (bond == 2) return 111;
		else if (bond == 3) return 103;
	}
	else if (element == "Mn") {
		if (bond == 1) return 119;
		else if (bond == 2) return 105;
		else if (bond == 3) return 103;
	}
	else if (element == "Fe") {
		if (bond == 1) return 116;
		else if (bond == 2) return 109;
		else if (bond == 3) return 102;
	}
	else if (element == "Co") {
		if (bond == 1) return 111;
		else if (bond == 2) return 103;
		else if (bond == 3) return 96;
	}
	else if (element == "Ni") {
		if (bond == 1) return 110;
		else if (bond == 2) return 101;
		else if (bond == 3) return 101;
	}
	else if (element == "Cu") {
		if (bond == 1) return 112;
		else if (bond == 2) return 115;
		else if (bond == 3) return 120;
	}
	else if (element == "Zn") {
		if (bond == 1) return 118;
		else if (bond == 2) return 120;
	}
	else if (element == "Ga") {
		if (bond == 1) return 124;
		else if (bond == 2) return 117;
		else if (bond == 3) return 121;
	}
	else if (element == "Ge") {
		if (bond == 1) return 121;
		else if (bond == 2) return 111;
		else if (bond == 3) return 114;
	}
	else if (element == "As") {
		if (bond == 1) return 121;
		else if (bond == 2) return 114;
		else if (bond == 3) return 106;
	}
	else if (element == "Se") {
		if (bond == 1) return 116;
		else if (bond == 2) return 107;
		else if (bond == 3) return 107;
	}
	else if (element == "Br") {
		if (bond == 1) return 114;
		else if (bond == 2) return 109;
		else if (bond == 3) return 110;
	}
	else if (element == "Kr") {
		if (bond == 1) return 117;
		else if (bond == 2) return 121;
		else if (bond == 3) return 108;
	}
	else if (element == "Rb") {
		if (bond == 1) return 210;
		else if (bond == 2) return 202;
	}
	else if (element == "Sr") {
		if (bond == 1) return 185;
		else if (bond == 2) return 157;
		else if (bond == 3) return 139;
	}
	else if (element == "Y") {
		if (bond == 1) return 163;
		else if (bond == 2) return 130;
		else if (bond == 3) return 124;
	}
	else if (element == "Zr") {
		if (bond == 1) return 154;
		else if (bond == 2) return 127;
		else if (bond == 3) return 121;
	}
	else if (element == "Nb") {
		if (bond == 1) return 147;
		else if (bond == 2) return 125;
		else if (bond == 3) return 116;
	}
	else if (element == "Mo") {
		if (bond == 1) return 138;
		else if (bond == 2) return 121;
		else if (bond == 3) return 113;
	}
	else if (element == "Tc") {
		if (bond == 1) return 128;
		else if (bond == 2) return 120;
		else if (bond == 3) return 110;
	}
	else if (element == "Ru") {
		if (bond == 1) return 125;
		else if (bond == 2) return 114;
		else if (bond == 3) return 103;
	}
	else if (element == "Rh") {
		if (bond == 1) return 125;
		else if (bond == 2) return 110;
		else if (bond == 3) return 106;
	}
	else if (element == "Pd") {
		if (bond == 1) return 120;
		else if (bond == 2) return 117;
		else if (bond == 3) return 112;
	}
	else if (element == "Ag") {
		if (bond == 1) return 128;
		else if (bond == 2) return 139;
		else if (bond == 3) return 137;
	}
	else if (element == "Cd") {
		if (bond == 1) return 136;
		else if (bond == 2) return 144;
	}
	else if (element == "In") {
		if (bond == 1) return 142;
		else if (bond == 2) return 136;
		else if (bond == 3) return 146;
	}
	else if (element == "Sn") {
		if (bond == 1) return 140;
		else if (bond == 2) return 130;
		else if (bond == 3) return 132;
	}
	else if (element == "Sb") {
		if (bond == 1) return 140;
		else if (bond == 2) return 133;
		else if (bond == 3) return 127;
	}
	else if (element == "Te") {
		if (bond == 1) return 136;
		else if (bond == 2) return 128;
		else if (bond == 3) return 121;
	}
	else if (element == "I") {
		if (bond == 1) return 133;
		else if (bond == 2) return 129;
		else if (bond == 3) return 125;
	}
	else if (element == "Xe") {
		if (bond == 1) return 131;
		else if (bond == 2) return 135;
		else if (bond == 3) return 122;
	}
	else if (element == "Cs") {
		if (bond == 1) return 232;
		else if (bond == 2) return 209;
	}
	else if (element == "Ba") {
		if (bond == 1) return 196;
		else if (bond == 2) return 161;
		else if (bond == 3) return 149;
	}
	else if (element == "La") {
		if (bond == 1) return 180;
		else if (bond == 2) return 139;
		else if (bond == 3) return 139;
	}
	else if (element == "Ce") {
		if (bond == 1) return 163;
		else if (bond == 2) return 137;
		else if (bond == 3) return 131;
	}
	else if (element == "Pr") {
		if (bond == 1) return 176;
		else if (bond == 2) return 138;
		else if (bond == 3) return 128;
	}
	else if (element == "Nd") {
		if (bond == 1) return 174;
		else if (bond == 2) return 137;
	}
	else if (element == "Pm") {
		if (bond == 1) return 173;
		else if (bond == 2) return 135;
	}
	else if (element == "Sm") {
		if (bond == 1) return 172;
		else if (bond == 2) return 134;
	}
	else if (element == "Eu") {
		if (bond == 1) return 168;
		else if (bond == 2) return 134;
	}
	else if (element == "Gd") {
		if (bond == 1) return 169;
		else if (bond == 2) return 135;
		else if (bond == 3) return 132;
	}
	else if (element == "Tb") {
		if (bond == 1) return 168;
		else if (bond == 2) return 135;
	}
	else if (element == "Dy") {
		if (bond == 1) return 167;
		else if (bond == 2) return 133;
	}
	else if (element == "Ho") {
		if (bond == 1) return 166;
		else if (bond == 2) return 133;
	}
	else if (element == "Er") {
		if (bond == 1) return 165;
		else if (bond == 2) return 133;
	}
	else if (element == "Tm") {
		if (bond == 1) return 164;
		else if (bond == 2) return 131;
	}
	else if (element == "Yd") {
		if (bond == 1) return 170;
		else if (bond == 2) return 129;
	}
	else if (element == "Lu") {
		if (bond == 1) return 162;
		else if (bond == 2) return 131;
		else if (bond == 3) return 131;
	}
	else if (element == "Hf") {
		if (bond == 1) return 152;
		else if (bond == 2) return 128;
		else if (bond == 3) return 122;
	}
	else if (element == "Ta") {
		if (bond == 1) return 146;
		else if (bond == 2) return 126;
		else if (bond == 3) return 119;
	}
	else if (element == "W") {
		if (bond == 1) return 137;
		else if (bond == 2) return 120;
		else if (bond == 3) return 115;
	}
	else if (element == "Re") {
		if (bond == 1) return 131;
		else if (bond == 2) return 119;
		else if (bond == 3) return 110;
	}
	else if (element == "Os") {
		if (bond == 1) return 129;
		else if (bond == 2) return 116;
		else if (bond == 3) return 109;
	}
	else if (element == "Ir") {
		if (bond == 1) return 122;
		else if (bond == 2) return 115;
		else if (bond == 3) return 107;
	}
	else if (element == "Pt") {
		if (bond == 1) return 123;
		else if (bond == 2) return 112;
		else if (bond == 3) return 110;
	}
	else if (element == "Au") {
		if (bond == 1) return 124;
		else if (bond == 2) return 121;
		else if (bond == 3) return 123;
	}
	else if (element == "Hg") {
		if (bond == 1) return 133;
		else if (bond == 2) return 142;
	}
	else if (element == "Tl") {
		if (bond == 1) return 144;
		else if (bond == 2) return 142;
		else if (bond == 3) return 150;
	}
	else if (element == "Pb") {
		if (bond == 1) return 144;
		else if (bond == 2) return 135;
		else if (bond == 3) return 137;
	}
	else if (element == "Bi") {
		if (bond == 1) return 151;
		else if (bond == 2) return 141;
		else if (bond == 3) return 135;
	}
	else if (element == "Po") {
		if (bond == 1) return 145;
		else if (bond == 2) return 135;
		else if (bond == 3) return 129;
	}
	else if (element == "At") {
		if (bond == 1) return 147;
		else if (bond == 2) return 138;
		else if (bond == 3) return 138;
	}
	else if (element == "Rn") {
		if (bond == 1) return 142;
		else if (bond == 2) return 145;
		else if (bond == 3) return 133;
	}
	else if (element == "Fr") {
		if (bond == 1) return 223;
		else if (bond == 2) return 218;
	}
	else if (element == "Ra") {
		if (bond == 1) return 201;
		else if (bond == 2) return 173;
		else if (bond == 3) return 159;
	}
	else if (element == "Ac") {
		if (bond == 1) return 186;
		else if (bond == 2) return 153;
		else if (bond == 3) return 140;

	}
	else if (element == "Th") {
		if (bond == 1) return 175;
		else if (bond == 2) return 143;
		else if (bond == 3) return 136;

	}
	else if (element == "Pa") {
		if (bond == 1) return 169;
		else if (bond == 2) return 138;
		else if (bond == 3) return 129;

	}
	else if (element == "U") {
		if (bond == 1) return 170;
		else if (bond == 2) return 134;
		else if (bond == 3) return 118;
	}
	else if (element == "Np") {
		if (bond == 1) return 171;
		else if (bond == 2) return 136;
		else if (bond == 3) return 116;
	}
	else if (element == "Pu") {
		if (bond == 1) return 172;
		else if (bond == 2) return 135;
	}
	else if (element == "Am") {
		if (bond == 1) return 166;
		else if (bond == 2) return 135;
	}
	else if (element == "Cm") {
		if (bond == 1) return 166;
		else if (bond == 2) return 136;
	}
	else if (element == "Bk") {
		if (bond == 1) return 168;
		else if (bond == 2) return 139;
	}
	else if (element == "Cf") {
		if (bond == 1) return 168;
		else if (bond == 2) return 140;
	}
	else if (element == "Es") {
		if (bond == 1) return 165;
		else if (bond == 2) return 140;
	}
	else if (element == "Fm") {
		if (bond == 1) return 167;
	}
	else if (element == "Md") {
		if (bond == 1) return 173;
		else if (bond == 2) return 139;
	}
	else if (element == "No") {
		if (bond == 1) return 176;
		else if (bond == 2) return 159;
	}
	else if (element == "Lr") {
		if (bond == 1) return 161;
		else if (bond == 2) return 141;
	}
	else if (element == "Rf") {
		if (bond == 1) return 157;
		else if (bond == 2) return 140;
		else if (bond == 3) return 131;

	}
	else if (element == "Db") {
		if (bond == 1) return 149;
		else if (bond == 2) return 136;
		else if (bond == 3) return 126;

	}
	else if (element == "Sg") {
		if (bond == 1) return 143;
		else if (bond == 2) return 128;
		else if (bond == 3) return 121;

	}
	else if (element == "Bh") {
		if (bond == 1) return 141;
		else if (bond == 2) return 128;
		else if (bond == 3) return 119;

	}
	else if (element == "Hs") {
		if (bond == 1) return 134;
		else if (bond == 2) return 125;
		else if (bond == 3) return 118;

	}
	else if (element == "Mt") {
		if (bond == 1) return 129;
		else if (bond == 2) return 125;
		else if (bond == 3) return 113;
	}
	else if (element == "Ds") {
		if (bond == 1) return 128;
		else if (bond == 2) return 116;
		else if (bond == 3) return 112;
	}
	else if (element == "Rg") {
		if (bond == 1) return 121;
		else if (bond == 2) return 116;
		else if (bond == 3) return 118;
	}
	else if (element == "112") {
		if (bond == 1) return 122;
		else if (bond == 2) return 137;
		else if (bond == 3) return 130;
	}
	else if (element == "113") {
		if (bond == 1) return 136;
	}
	else if (element == "114") {
		if (bond == 1) return 143;
	}
	else if (element == "115") {
		if (bond == 1) return 162;
	}
	else if (element == "116") {
		if (bond == 1) return 175;
	}
	else if (element == "117") {
		if (bond == 1) return 165;
	}
	else if (element == "118") {
		if (bond == 1) return 157;
	}
	return 0;
}


int vanDerWaalsRadii(std::string element)
{
	//source is in angstroms	, converted to pm by conversion of 1/100		  return 0.0;
	// Batsanov, S.S., Van der Waals Radii Evaluated from
	//Structural Parameters of Metals, Zh.Fiz.Khim., 2000,
	//	vol. 74, no. 7, pp. 1273–1276.
	//& Rowland, 1996 - crystallographic non metals. 
	if (element == "H") return 110;
	else if (element == "He") return 00;
	else if (element == "Li") return 224;
	else if (element == "Be") return 186;
	else if (element == "B")  return 174;
	else if (element == "C")  return 177;
	else if (element == "N")  return 164;
	else if (element == "O")  return 158;
	else if (element == "F")  return 146;
	else if (element == "Ne") return 00;
	else if (element == "Na") return 257;
	else if (element == "Mg") return 227;
	else if (element == "Al") return 221;
	else if (element == "Si") return 206;
	else if (element == "P")  return 00;
	else if (element == "S")  return 181;
	else if (element == "Cl") return 176;
	else if (element == "Ar") return 00;
	else if (element == "K")  return 300;
	else if (element == "Ca") return 261;
	else if (element == "Sc") return 228;
	else if (element == "Ti") return 214;
	else if (element == "V")  return 203;
	else if (element == "Cr") return 197;
	else if (element == "Mn") return 196;
	else if (element == "Fe") return 196;
	else if (element == "Co") return 195;
	else if (element == "Ni") return 194;
	else if (element == "Cu") return 200;
	else if (element == "Zn") return 202;
	else if (element == "Ga") return 208;
	else if (element == "Ge") return 213;
	else if (element == "As") return 216;
	else if (element == "Se") return 00;
	else if (element == "Br") return 187;
	else if (element == "Kr") return 00;
	else if (element == "Rb") return 312;
	else if (element == "Sr") return 278;
	else if (element == "Y")  return 245;
	else if (element == "Zr") return 225;
	else if (element == "Nb") return 213;
	else if (element == "Mo") return 206;
	else if (element == "Tc") return 204;
	else if (element == "Ru") return 202;
	else if (element == "Rh") return 202;
	else if (element == "Pd") return 205;
	else if (element == "Ag") return 213;
	else if (element == "Cd") return 117;
	else if (element == "In") return 224;
	else if (element == "Sn") return 229;
	else if (element == "Sb") return 233;
	else if (element == "Te") return 00;
	else if (element == "I")  return 203;
	else if (element == "Xe") return 00;
	else if (element == "Cs") return 331;
	else if (element == "Ba") return 285;
	else if (element == "La") return 251;
	else if (element == "Ce") return 00;
	else if (element == "Pr") return 00;
	else if (element == "Nd") return 00;
	else if (element == "Pm") return 00;
	else if (element == "Sm") return 00;
	else if (element == "Eu") return 00;
	else if (element == "Gd") return 00;
	else if (element == "Tb") return 00;
	else if (element == "Dy") return 00;
	else if (element == "Ho") return 00;
	else if (element == "Er") return 00;
	else if (element == "Tm") return 00;
	else if (element == "Yd") return 00;
	else if (element == "Lu") return 00;
	else if (element == "Hf") return 224;
	else if (element == "Ta") return 00;
	else if (element == "W")  return 207;
	else if (element == "Re") return 205;
	else if (element == "Os") return 203;
	else if (element == "Ir") return 203;
	else if (element == "Pt") return 206;
	else if (element == "Au") return 213;
	else if (element == "Hg") return 217;
	else if (element == "Tl") return 225;
	else if (element == "Pb") return 00;
	else if (element == "Bi") return 242;
	else if (element == "Po") return 00;
	else if (element == "At") return 00;
	else if (element == "Rn") return 00;
	else if (element == "Fr") return 00;
	else if (element == "Ra") return 00;
	else if (element == "Ac") return 00;
	else if (element == "Th") return 243;
	else if (element == "Pa") return 00;
	else if (element == "U")  return 217;
	else if (element == "Np") return 00;
	else if (element == "Pu") return 00;
	else if (element == "Am") return 00;
	else if (element == "Cm") return 00;
	else if (element == "Bk") return 00;
	else if (element == "Cf") return 00;
	else if (element == "Es") return 00;
	else if (element == "Fm") return 00;
	else if (element == "Md") return 00;
	else if (element == "No") return 00;
	else if (element == "Lr") return 00;
	else if (element == "Rf") return 00;
	else if (element == "Db") return 00;
	else if (element == "Sg") return 00;
	else if (element == "Bh") return 00;
	else if (element == "Hs") return 00;
	else if (element == "Mt") return 00;
	else if (element == "Ds") return 00;
	else if (element == "Rg") return 00;
	else if (element == "112")return 00;
	else if (element == "113")return 00;
	else if (element == "114")return 00;
	else if (element == "115")return 00;
	else if (element == "116")return 00;
	else if (element == "117")return 00;
	else if (element == "118")return 00;
	return 0;
}

bool radialCriteria(structure s, int percent)
{
	/*there will be a radial critera based on covalentand van der waals radii.
	//single, double and triple bonding must be allowed
	//another atom cannot get within the van der waals radii, unless its bonded

	//check if an atom is within the van der waals radii, if it is bonded. if not get rid of it

	//the percent refers to the accepted deviation from the covalent bonding radius
	*/

	//for each atom a in s
	for (int i = 0; i < s.set.size(); i++)
	{
		for (int a = 0; a < s.set[i].size(); a++)
		{
			//for each atom b in s
			for (int j = 0; j < s.set.size(); j++)
			{
				for (int b = 0; b < s.set[j].size(); b++)
				{
					if (i == j && b == a) {
					}
					else {
						//do not perform on the same atom
						double dist = euclideanDistance(s.set[i][a], s.set[j][b]);
						if (dist < vanDerWaalsRadii(s.elements[i]))
						{
							bool accept = false;
							double bl1 = covalentRadii(s.elements[i], 1);
							double bl1r = bl1 * ((double)percent / (double)100);
							if (dist > bl1 - bl1r && dist < bl1 + bl1r) accept = true;
							else {
								double bl2 = covalentRadii(s.elements[i], 2);
								if (bl2 != 0)
								{
									double bl2r = bl2 * ((double)percent / (double)100);
									if (dist > bl2 - bl2r && dist < bl2 + bl2r) accept = true;
									else {
										double bl3 = covalentRadii(s.elements[i], 3);
										if (bl3 != 0)
										{
											double bl3r = bl3 * ((double)percent / (double)100);
											if (dist > bl3 - bl3r && dist < bl3 + bl3r) accept = true;
										}
									}
								}




							}
							if (!accept) return false;
						}


					}

				}
			}

		}
	}
	//all passed, return true
	return true;


}
