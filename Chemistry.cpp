#include "chemistry.hpp"


atom::atom()
{
//i guess we have to have a default constructor
	x, y, z = 0;
	r, psi, theta = 0;
	types = 0;
	scores = nullptr;
}
atom::atom(double x, double y, double z) :x(x), y(y), z(z)
{
    r = std::sqrt(x * x + y * y + z * z);
    psi = std::acos((x / std::sqrt(x * x + y * y)));
    theta = z / r;
}
atom::atom(double x, double y, double z,int group) :x(x), y(y), z(z),group(group)
{
	r = std::sqrt(x * x + y * y + z * z);
	psi = std::acos((x / std::sqrt(x * x + y * y)));
	theta = z / r;
}
atom::atom(const atom& t):x(t.x),y(t.y),z(t.z),theta(t.theta),psi(t.psi),r(t.r),types(t.types)
{
	if (t.scores != nullptr)
	{
		this->scores = new double[t.types];
		for (int i = 0; i < t.types; i++)
		{
			this->scores[i] = t.scores[i];
		}
	}
	else {
		this->scores = nullptr;
	}
	
}
atom atom::operator=(const atom& t)
{
	this->x = t.x;
	this->y = t.y;
	this->z = t.z;
	this->theta = t.theta;
	this->psi = t.psi;
	this->r = t.r;
	this->types = t.types;
	if (t.scores != nullptr)
	{
		this->scores = new double[t.types];
		for (int i = 0; i < t.types; i++)
		{
			this->scores[i] = t.scores[i];
		}
	}
	else {
		this->scores = nullptr;
	}
	return *this;
}
atom::~atom()
{
	//std::cout << "atom deleter called" << std::endl;
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

atom atom::operator-(const atom& other) {
	return atom( this->x - other.x, this->y - other.y, this->z - other.z );
}

atom atom::cross(const atom& other){
	return atom(
		this->y * other.z - this->z * other.y,
		this->z * other.x - this->x * other.z,
		this->x * other.y - this->y * other.x
	);
}

std::vector<double> atom::vec3() { return{this->x,this->y,this->z}; }

double atom::dot(const atom& other) {
	return this->x * other.x + this->y * other.y + this->z * other.z;
}

atom atom::normalize(){
	double n = std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	return atom( this->x / n, this->y / n, this->z / n );
}


bool radiusCompare(atom a, atom b) { return a.r < b.r; };

double euclideanDistance(atom a, atom b) {
	double x = (a.x - b.x);
	double y = (a.y - b.y);
	double z = (a.z - b.z);
	double d = std::sqrt(x * x + y * y + z * z);
	return d;
}
double angle3atoms(atom a, atom b, atom c)
{
	double bx = b.x - a.x;
	double by = b.y - a.y;
	double bz = b.z - a.z;
	double cx = c.x - a.x;
	double cy = c.y - a.y;
	double cz = c.z - a.z;

	double dotproduct = bx * cx + by * cy + bz * cz;
	double magnitudebvector = std::sqrt(bx * bx + by * by + bz * bz);
	double magnitudecvector = std::sqrt(cx * cx + cy * cy + cz * cz);

	return std::acos(dotproduct / (magnitudebvector * magnitudecvector));
}

double* normalVectorPlane(atom a, atom b, atom c)
{
	//returns the angle between the three atoms and the normal vector of their plane
	
	//new coordinates for b and c
	double bx = b.x - a.x;
	double by = b.y - a.y;
	double bz = b.z - a.z;
	double cx = c.x - a.x;
	double cy = c.y - a.y;
	double cz = c.z - a.z;

	double dotproduct = bx * cx + by * cy + bz * cz;
	double magnitudebvector = std::sqrt(bx * bx + by * by + bz * bz);
	double magnitudecvector = std::sqrt(cx * cx + cy * cy + cz * cz);

	//double angle = std::acos(dotproduct / (magnitudebvector * magnitudecvector));
	double crossproducti = by*cz-bz*cy;//determinant of by bz, cy cz
	double crossproductj = bz * cx - bx * cz;//negative deterimant of bx bz, cx cz
	double crossproductk = bx * cy - cx * by;//determinant of bx by, cx cy

	double* normalVector = new double[3];
	normalVector[0] = crossproducti;
	normalVector[1] = crossproductj;
	normalVector[2] = crossproductk;

	return normalVector;


}

double dihedralAngle(double* normalVector, atom b, atom a)
{
	double bx = b.x - a.x;
	double by = b.y - a.y;
	double bz = b.z - a.z;

	double dotproduct = bx * normalVector[0] + by * normalVector[1] + bz * normalVector[2];
	double magnitudebvector = std::sqrt(bx * bx + by * by + bz * bz);
	double magnitudenvector = std::sqrt(normalVector[0] * normalVector[0] + normalVector[1] * normalVector[1] + normalVector[2] * normalVector[2]);

	return std::acos(dotproduct / (magnitudebvector * magnitudenvector));
}



structure::structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool center) :set(set), elements(elements)
{
	for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
	{
		composition.push_back(it->size());
	}
	//std::cout << "set sizes::" << std::endl;
	//for (int i = 0; i < set.size(); i++) std::cout << set[i].size() << " ";
	//std::cout << std::endl;
    double totalX = 0.0;
    double totalY = 0.0;
    double totalZ = 0.0;
    int n = 0;
	//Today I learned that these iterators create copies of everything and then delete them...
	//so if you want to use these iterators you need copy constructors for the objects inside that will not just copy pointers but instead recreate the objects
	if (center) {
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
			//composition.push_back(it->size());
			for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
			{
				itt->polarCenter(x, y, z);
			}
		}
	}
	//this->print();
	//this->zmatrix = nullptr;
	//this->zmatrixElements = nullptr;
}
structure::structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, std::vector<std::vector<int>> assignment, bool center )
{
	//the set input will have empty elements corresponding to the fragments, which must not be included
	for (std::vector<std::vector<atom>>::iterator at = set.begin(); at != set.end(); at++) if (at->size() != 0) {
		this->set.push_back(*at);
		this->composition.push_back(at->size());
	}
	for (std::vector<std::vector<int>>::iterator it = assignment.begin(); it != assignment.end(); it++) if (it->size() != 0) for(std::vector<int>::iterator itt = it->begin(); itt!= it->end(); itt++) this->assignment.push_back(*itt);

	for (std::vector<std::string>::iterator et = elements.begin(); et != elements.end(); et++) if (et->size() <= 2) this->elements.push_back(*et);//removes .xyz files

	//std::cout << "creating structure from fragments with assignment size: " << assignment.size() << std::endl;
	int a = 0;
	//std::cout << "sets: ";
	/*for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++)
	{
		std::cout << it->size() << ",";
	}
	std::cout << std::endl;
	for (int s = 0; s < this->set.size(); s++)
	{
		std::cout << "element: " << this->elements[s] << ", " << " assignments: ";
		for (int ss = 0; ss < this->set[s].size(); ss++)
		{
			//std::cout << ss << ":ss,s:" << s << ",a:" << a << std::endl;
			std::cout << this->assignment[a++] << " ";
		}
		std::cout << std::endl;

	}*/
	double totalX = 0.0;
	double totalY = 0.0;
	double totalZ = 0.0;
	int n = 0;



	if (center) {
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
			//composition.push_back(it->size());
			for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
			{
				itt->polarCenter(x, y, z);
			}
		}
	}
	//this->print();
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

            double newX = x * (2.0*(rand() % RAND_MAX) / RAND_MAX - 1);
            double newY = y * (2.0*(rand() % RAND_MAX) / RAND_MAX - 1);
            double newZ = z * (2.0*(rand() % RAND_MAX) / RAND_MAX - 1);
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
structure::structure(std::string input, std::vector<std::string> elements, std::vector<int> composition):elements(elements),composition(composition)
{
	bool bugging = 0;
	int lastSpace = 0;
	int space = 0;
	int element = -1;
	if(bugging) std::cout << " debugging: " << input << std::endl;
	if (input[space] == 'e')
	{
		//it says energy: first, ignore
		while (input[++space] != ' ') if (bugging) std::cout << "debugging? " << input[space] << std::endl;
	}
	lastSpace = space;
	if (bugging) std::cout << "input:" << input[space] << std::endl;
	while (input[++space] != ' ');
	this->energy = std::stod(input.substr(lastSpace, space - lastSpace));
	for (int e = 0; e < composition.size(); e++)
	{
		//for all types of elements
		std::vector<atom> element;
		for (int a = 0; a < composition[e]; a++)
		{
			//for each atom of type e
			double cx, cy, cz;
			lastSpace = space;
			while (input[++space] != ' '){}
			//passed the letter
			lastSpace = space;
			while (input[++space] != ' '){}
			cx = std::stof(input.substr(lastSpace, space - lastSpace));
			lastSpace = space;
			while (input[++space] != ' '){}
			cy = std::stof(input.substr(lastSpace, space - lastSpace));
			lastSpace = space;
			if (bugging) std::cout << "cy: " << cy << std::endl;
			while (input[++space] != ' ' && input[space] != '\0') {}
			if (bugging) std::cout << "another problem?" << std::endl;
			cz = std::stof(input.substr(lastSpace, space - lastSpace));
			element.push_back(atom(cx,cy,cz));
			if (bugging) std::cout << "cz: " << cz << std::endl;
		}
		this->set.push_back(element);
	}
}
structure::structure(structure& original, double variationX, double variationY, double variationZ, double cellX, double cellY, double cellZ, bool sort) :elements(original.elements)
{


    for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
    {
        std::vector<atom> atoms;
        for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
        {
            double varX = variationX * (2.0* rand() / RAND_MAX - 1);
			while((itt->x + varX) < -cellX || itt->x + varX > cellX)
			{
				varX = variationX * (2.0 * rand() / RAND_MAX - 1);
			}
            double varY = variationY * (2.0* rand() / RAND_MAX - 1);
			while ((itt->y + varY) < -cellY || itt->y + varY > cellY)
			{
				varY = variationY * (2.0 * rand() / RAND_MAX - 1);
			}
            double varZ = variationZ * (2.0* rand() / RAND_MAX - 1);
			while ((itt->z + varZ) < -cellZ || itt->z + varZ > cellZ)
			{
				varZ = variationZ * (2.0 * rand() / RAND_MAX - 1);
			}

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
structure::structure(structure& original, double variationX, double variationY, double variationZ, bool sort) :elements(original.elements)
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
structure::structure(structure& original, std::vector<double> variation) :elements(original.elements)
{
	//makes variants of variation
	int varIterator = 0;
	for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
	{
		std::vector<atom> atoms;
		for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
		{
			atom at = atom(itt->x + variation[varIterator], itt->y + variation[varIterator+1], itt->z + variation[varIterator+2]);
			atoms.push_back(at);
			varIterator += 3;
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
}
structure::structure(structure& original, std::vector<double> variation,std::vector<int> fix, bool bound, float xb, float yb, float zb) :elements(original.elements),assignment(fix)
{
	
	//std::cout << "MOVEMENT: ";
	//for (std::vector<double>::iterator it = variation.begin(); it != variation.end(); it++ ) std::cout << *it << " ";
	std::cout << std::endl;
	//this one differs as variation is size = fix.size()*3 instead of atoms.size()*3
	int atomIterator = 0;
	for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
	{
		std::vector<atom> atoms;
		for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
		{
			//std::cout << "atom iterator: " << atomIterator << std::endl;
			//std::cout << "variations for atom " << atomIterator << ": " << variation[fix[atomIterator] * 3] << " " << variation[fix[atomIterator] * 3 + 1] << " " <<  variation[fix[atomIterator] * 3 + 2] << std::endl;
			//std::cout << "indexes for atom " << atomIterator << ": " << fix[atomIterator] * 3 << " " << fix[atomIterator] * 3 + 1 << " " << fix[atomIterator] * 3 + 2 << std::endl;
			double ax, ay, az;
			if (fix[atomIterator] < 0)
			{
				//this is used for the support, so it never moves.
				ax = itt->x;
				ay = itt->y;
				az = itt->z;
				//std::cout << "NOT MOVING SUPPORT" << std::endl;
			}
			else {
				//std::cout << "MOVING NON SUPPORT index " << fix[atomIterator] << " from atom iterator " << atomIterator << std::endl;
				//std::cout << "variation: " << variation[fix[atomIterator] * 3] << " " << variation[fix[atomIterator] * 3 + 1] << " " << variation[fix[atomIterator] * 3 + 2] << std::endl;
				ax = itt->x + variation[fix[atomIterator] * 3];
				ay = itt->y + variation[fix[atomIterator] * 3 + 1];
				az = itt->z + variation[fix[atomIterator] * 3 + 2];
				if (bound) {
					if (ax > xb / 2) ax -= xb;
					if (ax < -xb / 2) ax += xb;
					if (ay > yb / 2) ay -= yb;
					if (ay < -yb / 2) ay += yb;
					if (az > zb / 2) az -= zb;
					if (az < -zb / 2) az += zb;
				}
			}

			atom at = atom(ax,ay,az);
			atoms.push_back(at);

			atomIterator += 1;
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
}
structure::structure(structure& original, std::vector<double> variation, double dx, double dy, double dz) :elements(original.elements)
{
	//makes variants of at most variation dx,dy,dz based on variation vector with values ranging 1 to -1
	int varIterator = 0;
	for (std::vector<std::vector<atom>>::iterator it = original.set.begin(); it != original.set.end(); it++)
	{
		std::vector<atom> atoms;
		for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
		{
			atom at = atom(itt->x + dx * variation[varIterator], itt->y + dy * variation[varIterator+1], itt->z + dz * variation[varIterator+2]);
			varIterator += 3;
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
}
/*structure::structure(const structure& other) :elements(other.elements), assignment(other.assignment), set(other.set), zmatrix(nullptr), zmatrixElements(nullptr)
{
	std::cout << "copy constructor called" << std::endl;
}*/
/*structure::structure(structure& original, double variationX, double variationY, double variationZ, double rangeX, double rangeY, double rangeZ, bool sort) :elements(original.elements)
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
}*/

using iterator_category = std::forward_iterator_tag;
using difference_type = std::ptrdiff_t;
using value_type = atom;
using pointer = atom*;
using reference = atom&;




structure::Iterator::Iterator(pointer ptr,int indexI, int indexJ,std::vector<std::vector<atom>>& set) : m_ptr(ptr),indexI(indexI),indexJ(indexJ),set(set) {}
reference structure::Iterator::operator*() const { return *m_ptr; }
pointer structure::Iterator::operator->() { return m_ptr; }
structure::Iterator& structure::Iterator::operator++() { 
	//std::cout << "iteration; indexI: " << indexI << " indexI.size(): " << set[indexI].size() << " , IndexJ : " << indexJ << std::endl;
	if (indexJ + 1 >= set[indexI].size() && indexI + 1 >= set.size())
	{
		//end pointer
		indexI = -1;
		indexJ = -1;
		m_ptr = nullptr;
	}
	else if (indexJ + 1 >= set[indexI].size())
	{
		indexI = indexI + 1;
		indexJ = 0;
		m_ptr = &set[indexI][indexJ];
	}
	else {
		indexJ++;
		m_ptr = &set[indexI][indexJ];
	}
	return *this; 
}
structure::Iterator structure::Iterator::operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
bool structure::Iterator::operator== (const structure::Iterator& a) { //return a.m_ptr == m_ptr;
	return (a.indexI == indexI && a.indexJ == indexJ);
};
bool structure::Iterator::operator!= (const structure::Iterator& a) { 
	//std::cout << "!=: " << a.indexI << " " << a.indexJ << " vs " << indexI << " " << indexJ << " " << " pointsers: " << a.m_ptr << " " << m_ptr  << std::endl;
	/*if (a.indexI == indexI && a.indexJ == indexJ)
	{
		std::cout << "a vs b: X:" << (a.m_ptr)->x << " vs " << (m_ptr)->x << ",Y: " << (a.m_ptr)->y << " vs " << (m_ptr)->y << ",z: " << (a.m_ptr)->z << " vs " << (m_ptr)->z << std::endl;
	}
	*/
	//return a.m_ptr != m_ptr; 
	return !(a.indexI == indexI && a.indexJ == indexJ);
};


//first and last values
structure::Iterator structure::begin() { return structure::Iterator(&set[0][0],0,0,set); }
structure::Iterator structure::end() {
	//int a = set.size() - 1;
	//int b = set[a].size() - 1;
	return structure::Iterator(nullptr,-1,-1,set); 
	//technically the end is not supposed to be the last atom, but in this case we will set it to that to avoid other issues
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

void structure::print()
{
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			std::cout << this->elements[i] << " " << this->set[i][j].x << " " << this->set[i][j].y << " " << this->set[i][j].z << std::endl;
		}
	}
}


void structure::makeZmatrix()
{
	//create the z matrix
	if (zmatrix == nullptr)
	{
		zmatrixsize = 0;
		for (int i = 0; i < this->set.size(); i++) zmatrixsize += this->set[i].size();
		zmatrix = new double* [zmatrixsize];
		zmatrixElements = new std::string [zmatrixsize];
		for (int i = 0; i < zmatrixsize; i++) zmatrix[i] = new double[3];//bond length in angstroms, angle, dihedral angle

		double xshift = 0;
		double yshift = 0;
		double zshift = 0;

		double a2rx = 0;
		double a2ry = 0;
		double a2rz = 0;

		double a3rx = 0;
		double a3ry = 0;
		double a3rz = 0;

		int catom = 1;
		double* normalVector = nullptr;


		//reorder by euclidean distance
		std::vector<std::tuple<double, atom, std::string>> ordered;

		for (int i = 0; i < this->elements.size(); i++)
		{
			for (int j = 0; j < this->set[i].size(); j++)
			{
				ordered.push_back(std::tuple<double,atom,std::string>(euclideanDistance(this->set[0][0], this->set[i][j]), this->set[i][j], this->elements[i]));
			}
		}
		std::sort(ordered.begin(), ordered.end(), [](std::tuple<double, atom, std::string>& left, std::tuple<double, atom, std::string>& right) {
			return std::get<0>(left) < std::get<0>(right);
		});


		for (int i = 0; i < ordered.size(); i++)
		{
				switch (catom)
				{
				case 1:

					zmatrix[0][0] = 0;
					zmatrix[0][1] = 0;
					zmatrix[0][2] = 0;

					//no data for this atom, this is atom 0
					break;
				case 2:

					zmatrix[catom - 1][0] = std::get<0>(ordered[catom-1]);
					zmatrix[catom - 1][1] = 0;
					zmatrix[catom - 1][2] = 0;

					break;
				case 3:

					normalVector = normalVectorPlane(std::get<1>(ordered[0]), std::get<1>(ordered[1]), std::get<1>(ordered[2]));

					zmatrix[catom - 1][0] = std::get<0>(ordered[catom - 1]);
					zmatrix[catom - 1][1] = (180.0/M_PI)*angle3atoms(std::get<1>(ordered[0]), std::get<1>(ordered[1]), std::get<1>(ordered[catom-1]));
					zmatrix[catom - 1][2] = 0;

					break;
				default:
					//greater than 4
					zmatrix[catom - 1][0] = std::get<0>(ordered[catom - 1]);
					zmatrix[catom - 1][1] = (180.0 / M_PI)*angle3atoms(std::get<1>(ordered[0]), std::get<1>(ordered[1]), std::get<1>(ordered[catom - 1]));
					zmatrix[catom - 1][2] = (180.0 / M_PI)*dihedralAngle(normalVector, std::get<1>(ordered[catom - 1]), std::get<1>(ordered[0]));
					break;
				}
				zmatrixElements[catom - 1] = std::get<2>(ordered[catom - 1]);
				catom++;
		}

		delete[] normalVector;


	}


	//delete the z matrix
}

structure::~structure()
{
	//std::cout << "DELETE STRUCTURE" << std::endl;
	//delete the z matrix
	if (zmatrix != nullptr) {
		for (int i = 0; i < zmatrixsize; i++) delete[] zmatrix[i];
		delete[] zmatrix;
	}
	delete[] zmatrixElements;
}

void structure::printZmatrix()
{
	int zmatrixiterator = 0;
	std::cout << "Zmatrix: Ra - Rn , angle A-B-n ,  dihedral angle A-B-C-n" << std::endl;
	if (zmatrixsize == 0 || zmatrix == nullptr) this->makeZmatrix();
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			std::cout << this->elements[i] << ": " << zmatrix[zmatrixiterator][0] << " " << zmatrix[zmatrixiterator][1] << " " << zmatrix[zmatrixiterator][2] << std::endl;
			zmatrixiterator++;
		}
	}
}

void structure::writeToObj(std::string outputName)
{
	
	int natoms = 0;
	//find center of mass
	double tx = 0;
	double ty = 0;
	double tz = 0;
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			tx += this->set[i][j].x;
			ty += this->set[i][j].y;
			tz += this->set[i][j].z;
			natoms++;
		}
	}
	tx = tx / natoms;
	ty = ty / natoms;
	tz = tz / natoms;

	/*
	
	*/
	std::ofstream outFile(outputName);
	if (!outFile.is_open()) {
		std::cout << "Error: Unable to open file " << outputName << " for writing." << std::endl;
		return;
	}

	const double PI = 3.14159265359;

	outFile << "# Structure OBJ File" << std::endl;

	int vertexOffset = 1;

	int numSegments = 32;


	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			double radius = atomicRadius(this->elements[i]);

			double x_offset = 1 * (this->set[i][j].x - tx);
			double y_offset = 1 * (this->set[i][j].y - ty);
			double z_offset = 1 * (this->set[i][j].z - tz);
			std::cout << x_offset << " " << y_offset << " " << z_offset << std::endl;
			//mtl file
			outFile << "mtllib " << this->elements[i] << ".mtl" << std::endl;
			//vertices
			for (int j = 0; j <= numSegments; ++j) {
				for (int i = 0; i <= numSegments; ++i) {
					double theta = (double)j / numSegments * PI;
					double phi = (double)i / numSegments * 2 * PI;
					double vx = x_offset + radius * sin(theta) * cos(phi);
					double vy = y_offset + radius * sin(theta) * sin(phi);
					double vz = z_offset + radius * cos(theta);
					outFile << "v " << vx << " " << vy << " " << vz << std::endl;
				}
			}

			//faces
			for (int j = 0; j < numSegments; ++j) {
				for (int i = 0; i < numSegments; ++i) {
					int first = vertexOffset + (j * (numSegments + 1)) + i;
					int second = first + numSegments + 1;
					outFile << "f " << first << " " << second << " " << (second + 1) % ((numSegments + 1) * (numSegments + 1)) + vertexOffset << std::endl;
					outFile << "f " << first << " " << (second + 1) % ((numSegments + 1) * (numSegments + 1)) + vertexOffset << " " << (first + 1) % ((numSegments + 1) * (numSegments + 1)) + vertexOffset << std::endl;
				}
			}

			vertexOffset += (numSegments + 1) * (numSegments + 1);
		}
	}

	outFile.close();
	

}

void structure::writeToXyz(std::string outputName)
{

	int natoms = 0;
	//find center of mass
	double tx = 0;
	double ty = 0;
	double tz = 0;
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			tx += this->set[i][j].x;
			ty += this->set[i][j].y;
			tz += this->set[i][j].z;
			natoms++;
		}
	}
	tx = tx / natoms;
	ty = ty / natoms;
	tz = tz / natoms;

	/*

	*/
	std::ofstream outFile(outputName);
	if (!outFile.is_open()) {
		std::cout << "Error: Unable to open file " << outputName << " for writing." << std::endl;
		return;
	}

	const double PI = 3.14159265359;

	//outFile << "# Structure OBJ File" << std::endl;

	outFile << std::to_string(natoms) << std::endl << "Energy: " << this->energy << " , time: " << this->time;
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			outFile << "\n" << this->elements[i] << " " << std::to_string(this->set[i][j].x - tx) << " " << std::to_string(this->set[i][j].y - ty) << " " << std::to_string(this->set[i][j].z - tz);
		}
	}

	outFile.close();


}
std::string structure::xyz()
{
	int natoms = 0;
	//find center of mass
	double tx = 0;
	double ty = 0;
	double tz = 0;
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			tx += this->set[i][j].x;
			ty += this->set[i][j].y;
			tz += this->set[i][j].z;
			natoms++;
		}
	}
	tx = tx / natoms;
	ty = ty / natoms;
	tz = tz / natoms;

	/*

	*/
	

	const double PI = 3.14159265359;

	//outFile << "# Structure OBJ File" << std::endl;
	std::string retVal = std::to_string(natoms) + "\nEnergy: " + std::to_string(this->energy) + " , time: " + std::to_string(this->time);
	for (int i = 0; i < this->set.size(); i++)
	{
		for (int j = 0; j < this->set[i].size(); j++)
		{
			retVal += "\n" + this->elements[i] + " " + std::to_string(this->set[i][j].x - tx) + " " + std::to_string(this->set[i][j].y - ty) + " " + std::to_string(this->set[i][j].z - tz);
		}
	}

	return retVal;
}
std::vector<structure*>  XyzToStructures(std::string fileName,std::vector<std::string> elements)
{
	std::ifstream file(fileName);
	std::string line;
	std::vector<structure*> retValue;
	int numAtoms = 0;
	while (std::getline(file, line))
	{
		if (numAtoms == 0) numAtoms = std::stoi(line);//the first line is a number
		if (std::stoi(line) == numAtoms)
		{
		std::vector<int> composition;

		std::getline(file, line);//blank line
		std::vector<std::vector<atom>> set;
		for (int i = 0; i < elements.size(); i++) set.push_back({});
		for (int i = 0; i < numAtoms; i++)
		{
			std::getline(file, line);
			std::string element;
			std::string x;
			std::string y;
			std::string z;
			int nspace = 0;
			for (int c = 0; c < line.size(); c++)
			{
				if (line[c] == ' ')
				{
					nspace++;
				}
				else {
					switch (nspace)
					{
					case 0:
						element += line[c];
						break;
					case 1:
						x += line[c];
						break;
					case 2:
						y += line[c];
						break;
					case 3:
						z += line[c];
						break;
					default:
						break;
					}
				}
			}
			int e = 0;
			while (element != elements[e]) { e++; };
			set[e].push_back(atom(std::stof(x), std::stof(y), std::stof(z)));
		}

		structure* newStructure = new structure(set, elements);
		retValue.push_back(newStructure);
		}
		else {
			int lineskip = std::stoi(line);
			for(int i = 0; i <= lineskip; i++) std::getline(file, line);
			
		}
	}

	file.close();
	return retValue;
}

structure* merge(structure* a, structure* b, double shiftZb)
{
	//this function is intended to merge a support (structure b) to prior structure and can shift the support on the z axis if necessary.
	std::vector<std::string> elements = {};
	std::vector<std::vector<atom>> set = {};
	std::vector<std::vector<int>> assignment = {};
	for (int e = 0; e < b->elements.size(); e++)
	{
		//check if this element is already in elements
		bool contained = false;
		int elementIndex = -1;
		for (int ee = 0; ee <  elements.size();  ee++)
		{
			if (elements[ee] == b->elements[e])
			{
				contained = true;
				elementIndex = ee;
			}
		}
		if (!contained) {
			elementIndex = elements.size();//size is about to increase by 1 making current size the index.
			elements.push_back(b->elements[e]);
			set.push_back({});
			assignment.push_back({});
		}

		for (std::vector<atom>::iterator at = b->set[e].begin(); at != b->set[e].end(); at++)
		{
			//std::cout << "adding element type " << e << " in merge from b" << std::endl;
			set[elementIndex].push_back(atom(at->x, at->y, at->z + shiftZb));
			assignment[elementIndex].push_back(-1);//-1 for support group
		}

	}

	int aindex = 0;
	for (int e = 0; e < a->elements.size(); e++)
	{
		//check if this element is already in elements
		bool contained = false;
		int elementIndex = -1;
		for (int ee = 0; ee < elements.size(); ee++)
		{
			if (elements[ee] == a->elements[e])
			{
				contained = true;
				elementIndex = ee;
			}
		}
		if (!contained) {
			elementIndex = elements.size();//size is about to increase by 1 making current size the index.
			elements.push_back(a->elements[e]);
			set.push_back({});
			assignment.push_back({});
		}

		for (std::vector<atom>::iterator at = a->set[e].begin(); at != a->set[e].end(); at++) {
			//std::cout << "adding atom from a" << std::endl;
			set[elementIndex].push_back(atom(at->x, at->y, at->z));
			assignment[elementIndex].push_back(a->assignment[aindex++]);
		}
	}
	//std::cout << "set size:" << set.size() << std::endl;
	structure* retval = new structure(set, elements,assignment);
	//retval->assignment = assignment;
	//retval->zmatrix = nullptr;
	//retval->zmatrixElements = nullptr;
	return retval;


}

//chemistry
int covalentRadii(std::string element, int bond, bool horizontal)
{
	//std::cout << element << std::endl;
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
	if (bond == 1) return 146;
	else if (bond == 2) return 146;// 114;
	else if (bond == 3) return 146;// 103;
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
	else {
		//this is not an atom, instead its a collection of atoms in an xyz file
		int horizontalCovalentRadii = 0;
		int verticalCovalentRadii = 0;
		std::ifstream file(element); // Open the file # of hole location and element distribution combinations:
		std::string line;
		std::getline(file, line);//num atoms
		std::getline(file, line);//horizontal_vertical
		std::string hori = "";
		std::string vert = "";
		bool verttime = false;
		//std::cout << line << std::endl;
		for (int s = 0; s < line.size(); s++)
		{
			//std::cout << line[s] << std::endl;
			if (line[s] == ' ') verttime = true;
			else if (verttime) vert += line[s];
			else hori += line[s];
		}
		file.close();
		if (horizontal) return std::stoi(hori);
		else return std::stoi(vert);

		
	
		}
	return 0;
}

std::string atomicNumber(int num)
{
	switch (num)
	{
	case 1:
		return "H";
	case 2:
		return "He";
	case 3:
		return "Li";
	case 4:
		return "Be";
	case 5:
		return "B";
	case 6:
		return "C";
	case 7:
		return "N";
	case 8:
		return "O";
	case 9:
		return "F";
	case 10:
		return "Ne";
	case 11:
		return "Na";
	case 12:
		return "Mg";
	case 13:
		return "Al";
	case 14:
		return "Si";
	case 15:
		return "P";
	case 16:
		return "S";
	case 17:
		return "Cl";
	case 18:
		return "Ar";
	case 19:
		return "K";
	case 20:
		return "Ca";
	case 21:
		return "Sc";
	case 22:
		return "Ti";
	case 23:
		return "V";
	case 24:
		return "Cr";
	case 25:
		return "Mn";
	case 26:
		return "Fe";
	case 27:
		return "Co";
	case 28:
		return "Ni";
	case 29:
		return "Cu";
	case 30:
		return "Zn";
	case 31:
		return "Ga";
	case 32:
		return "Ge";
	case 33:
		return "As";
	case 34:
		return "Se";
	case 35:
		return "Br";
	case 36:
		return "Kr";
	case 37:
		return "Rb";
	case 38:
		return "Sr";
	case 39:
		return "Y";
	case 40:
		return "Zr";
	case 41:
		return "Nb";
	case 42:
		return "Mo";
	case 43:
		return "Tc";
	case 44:
		return "Ru";
	case 45:
		return "Rh";
	case 46:
		return "Pd";
	case 47:
		return "Ag";
	case 48:
		return "Cd";
	case 49:
		return "In";
	case 50:
		return "Sn";
	case 51:
		return "Sb";
	case 52:
		return "Te";
	case 53:
		return "I";
	case 54:
		return "Xe";
	case 55:
		return "Cs";
	case 56:
		return "Ba";
	case 57:
		return "La";
	case 58:
		return "Ce";
	case 59:
		return "Pr";
	case 60:
		return "Nd";
	case 61:
		return "Pm";
	case 62:
		return "Sm";
	case 63:
		return "Eu";
	case 64:
		return "Gd";
	case 65:
		return "Tb";
	case 66:
		return "Dy";
	case 67:
		return "Ho";
	case 68:
		return "Er";
	case 69:
		return "Tm";
	case 70:
		return "Yb";
	case 71:
		return "Lu";
	case 72:
		return "Hf";
	case 73:
		return "Ta";
	case 74:
		return "W";
	case 75:
		return "Re";
	case 76:
		return "Os";
	case 77:
		return "Ir";
	case 78:
		return "Pt";
	case 79:
		return "Au";
	case 80:
		return "Hg";
	case 81:
		return "Tl";
	case 82:
		return "Pb";
	case 83:
		return "Bi";
	case 84:
		return "Po";
	case 85:
		return "At";
	case 86:
		return "Rn";
	case 87:
		return "Fr";
	case 88:
		return "Ra";
	case 89:
		return "Ac";
	case 90:
		return "Th";
	case 91:
		return "Pa";
	case 92:
		return "U";
	case 93:
		return "Np";
	case 94:
		return "Pu";
	case 95:
		return "Am";
	case 96:
		return "Cm";
	case 97:
		return "Bk";
	case 98:
		return "Cf";
	case 99:
		return "Es";
	case 100:
		return "Fm";
	case 101:
		return "Md";
	case 102:
		return "No";
	case 103:
		return "Lr";
	case 104:
		return "Rf";
	case 105:
		return "Db";
	case 106:
		return "Sg";
	case 107:
		return "Bh";
	case 108:
		return "Hs";
	case 109:
		return "Mt";
	case 110:
		return "Ds";
	case 111:
		return "Rg";
	case 112:
		return "Cn";
	case 113:
		return "Nh";
	case 114:
		return "Fl";
	case 115:
		return "Mc";
	case 116:
		return "Lv";
	case 117:
		return "Ts";
	case 118:
		return "Og";
	default:
		return "Unknown";
	}
	return "todo";
}

double atomicRadius(std::string atomName)
{
	std::map<std::string, double> atomicRadii{
		{"H", 0.53},
		{"He", 0.31},
		{"Li", 1.67},
		{"Be", 1.12},
		{"B", 0.87},
		{"C", 0.67},
		{"N", 0.56},
		{"O", 0.48},
		{"F", 0.42},
		{"Ne", 0.38},
		{"Na", 1.90},
		{"Mg", 1.45},
		{"Al", 1.18},
		{"Si", 1.11},
		{"P", 0.98},
		{"S", 0.88},
		{"Cl", 0.79},
		{"Ar", 0.71},
		{"K", 2.43},
		{"Ca", 1.94},
		{"Sc", 1.84},
		{"Ti", 1.76},
		{"V", 1.71},
		{"Cr", 1.66},
		{"Mn", 1.61},
		{"Fe", 1.56},
		{"Co", 1.52},
		{"Ni", 1.49},
		{"Cu", 1.45},
		{"Zn", 1.42},
		{"Ga", 1.36},
		{"Ge", 1.25},
		{"As", 1.14},
		{"Se", 1.03},
		{"Br", 0.94},
		{"Kr", 0.88},
		{"Rb", 2.65},
		{"Sr", 2.19},
		{"Y", 2.12},
		{"Zr", 2.06},
		{"Nb", 1.98},
		{"Mo", 1.90},
		{"Tc", 1.83},
		{"Ru", 1.78},
		{"Rh", 1.73},
		{"Pd", 1.69},
		{"Ag", 1.65},
		{"Cd", 1.61},
		{"In", 1.56},
		{"Sn", 1.45},
		{"Sb", 1.33},
		{"Te", 1.23},
		{"I", 1.15},
		{"Xe", 1.08},
		{"Cs", 2.98},
		{"Ba", 2.53},
		{"La", 1.95},
		{"Ce", 1.85},
		{"Pr", 2.47},
		{"Nd", 2.06},
		{"Pm", 2.05},
		{"Sm", 2.38},
		{"Eu", 2.31},
		{"Gd", 2.33},
		{"Tb", 2.25},
		{"Dy", 2.28},
		{"Ho", 2.26},
		{"Er", 2.26},
		{"Tm", 2.22},
		{"Yb", 2.22},
		{"Lu", 2.17},
		{"Hf", 2.08},
		{"Ta", 2.00},
		{"W", 1.93},
		{"Re", 1.88},
		{"Os", 1.85},
		{"Ir", 1.80},
		{"Pt", 1.77},
		{"Au", 1.74},
		{"Hg", 1.71},
		{"Tl", 1.56},
		{"Pb", 1.54},
		{"Bi", 1.43},
		{"Po", 1.35},
		{"At", 1.27},
		{"Rn", 1.20},
		{"Fr", 2.87},
		{"Ra", 2.43},
		{"Ac", 1.95},
		{"Th", 1.80},
		{"Pa", 1.80},
		{"U", 1.75}
	};
	return atomicRadii.find(atomName)->second;
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

double picometersToAngstrom(int picometers)
{
	return (double)picometers * 0.01;
}

bool radialCriteria(structure& s, int percent)
{



	//std::cout << "RCp: " << s.set.size() << std::endl;
	//bond requirement forces at least that number of covalent bonds per atom. recommended value is 1
	/*there will be a radial critera based on covalentand van der waals radii.
	//single, double and triple bonding must be allowed
	//another atom cannot get within the van der waals radii, unless its bonded

	//check if an atom is within the van der waals radii, if it is bonded. if not get rid of it

	//the percent refers to the accepted deviation from the covalent bonding radius
	*/
	bool bugging = false;//debug output
	std::string failureCause = "";
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
						
						double diameterVDW = picometersToAngstrom(vanDerWaalsRadii(s.elements[i])) + picometersToAngstrom(vanDerWaalsRadii(s.elements[j]));
						double diameterR1 = picometersToAngstrom(covalentRadii(s.elements[i], 1)) + picometersToAngstrom(covalentRadii(s.elements[j], 1));
						double diameterR2 = picometersToAngstrom(covalentRadii(s.elements[i], 2)) + picometersToAngstrom(covalentRadii(s.elements[j], 2));
						if (picometersToAngstrom(covalentRadii(s.elements[i], 2)) == 0 || picometersToAngstrom(covalentRadii(s.elements[j], 2)) == 0) diameterR2 = 0;
						double diameterR3 = picometersToAngstrom(covalentRadii(s.elements[i], 3)) + picometersToAngstrom(covalentRadii(s.elements[j], 3));
						if (picometersToAngstrom(covalentRadii(s.elements[i], 3)) == 0 || picometersToAngstrom(covalentRadii(s.elements[j], 3)) == 0) diameterR2 = 0;
						double boundR1 = diameterR1 * ((double)percent / (double)100);
						double boundR2 = diameterR2 * ((double)percent / (double)100);
						double boundR3 = diameterR3 * ((double)percent / (double)100);
						//do not perform on the same atom
						double dist = euclideanDistance(s.set[i][a], s.set[j][b]);

						//std::cout << "atom " << a << " to " << b << " distance " << dist << ", diamterR1 - boundR1: " << diameterR1 - boundR1 << ", diameterR1 + boundR1: " << diameterR1 + boundR1 << std::endl;

						if (dist < diameterVDW)
						{
							bool accept = false;
							if (dist > diameterR1 - boundR1 && dist < diameterR1 + boundR1) accept = true;
							else {
								if (diameterR2 != 0)
								{
									if (dist > diameterR2 - boundR2 && dist < diameterR2 + boundR2) accept = true;
									else {
										if (diameterR3 != 0)
										{
											if (dist > diameterR3 - boundR3 && dist < diameterR3 + boundR3) accept = true;
										}
									}
								}
								else {
									if (bugging) failureCause += "bl2 analysis skipped for missing data\n";
									
								}




							}
							if (!accept) {
								if(bugging) std::cout << "failed for no covalent bond between " << s.elements[i] << a << " and " << s.elements[j] << b << std::endl;
								if (bugging) s.print();
								return false;
							}
						}
						//else if (bugging) std::cout << "outside of van der waals for " << s.elements[i] << a << " and " << s.elements[j] << b << " being longer than the van der waals radii," << picometersToAngstrom(vanDerWaalsRadii(s.elements[i])) << ", of " << s.elements[i] << a << " at " << dist << std::endl;

						
					}
					
				}
				
			}
			

		}
		
	}
	//all passed, return true
	return true;


}

bool covalentCriteria(structure* s, int percent)
{
	/*
	This only checks that atoms are not within the covalent radius, so that they do not become too close during steps and create poor structures
	
	*/
	bool bugging = false;//debug output
	std::string failureCause = "";
	//for each atom a in s
	for (int i = 0; i < s->set.size(); i++)
	{
		for (int a = 0; a < s->set[i].size(); a++)
		{
			//for each atom b in s
			int covalentBonds = 0;
			for (int j = 0; j < s->set.size(); j++)
			{
				for (int b = 0; b < s->set[j].size(); b++)
				{
					if (i == j && b == a) {
					}
					else {

						double diameterR1 = picometersToAngstrom(covalentRadii(s->elements[i], 1)) + picometersToAngstrom(covalentRadii(s->elements[j], 1));
						double diameterR2 = picometersToAngstrom(covalentRadii(s->elements[i], 2)) + picometersToAngstrom(covalentRadii(s->elements[j], 2));
						if (picometersToAngstrom(covalentRadii(s->elements[i], 2)) == 0 || picometersToAngstrom(covalentRadii(s->elements[j], 2)) == 0) diameterR2 = 0;
						double diameterR3 = picometersToAngstrom(covalentRadii(s->elements[i], 3)) + picometersToAngstrom(covalentRadii(s->elements[j], 3));
						if (picometersToAngstrom(covalentRadii(s->elements[i], 3)) == 0 || picometersToAngstrom(covalentRadii(s->elements[j], 3)) == 0) diameterR2 = 0;
						double boundR1 = diameterR1 * ((double)percent / (double)100);
						double boundR2 = diameterR2 * ((double)percent / (double)100);
						double boundR3 = diameterR3 * ((double)percent / (double)100);
						//do not perform on the same atom
						double dist = euclideanDistance(s->set[i][a], s->set[j][b]);

						double diameterVDW = picometersToAngstrom(vanDerWaalsRadii(s->elements[i])) + picometersToAngstrom(vanDerWaalsRadii(s->elements[j]));
							bool accept = false;
							//std::cout << "atom " << a << " to " << b << " distance " << dist << ", diamterR1 - boundR1: " << diameterR1 - boundR1 << ", diameterR1 + boundR1: " << diameterR1 + boundR1 << std::endl;

							if (dist > diameterR1 - boundR1) accept = true;
							else if(false) {
								if (diameterR2 != 0)
								{
									if (dist > diameterR2 - boundR2) accept = true;
									else {
										if (diameterR3 != 0)
										{
											if (dist > diameterR3 - boundR3) accept = true;
										}
									}
								}
								else {
									if (bugging) failureCause += "bl2 analysis skipped for missing data\n";

								}




							}
							if (!accept) {
								return false;
							}
						//else if (bugging) std::cout << "outside of van der waals for " << s.elements[i] << a << " and " << s.elements[j] << b << " being longer than the van der waals radii," << picometersToAngstrom(vanDerWaalsRadii(s.elements[i])) << ", of " << s.elements[i] << a << " at " << dist << std::endl;


					}

				}
			}

		}
	}
	//all passed, return true
	return true;


}


bool bondingRequirement(structure* s, int percent, bool numBonds)
{
	// step 1: produce a graph of the structure where nodes are atoms and edges represent covalent bonds
	//get number of atoms
	int natoms = 0;
	for (int i = 0; i < s->set.size(); i++) natoms += s->set[i].size();
	if (natoms == 0) {
		std::cout << "cryptid pointer error detected from createXnew loop in BHC algo" << std::endl;
		/*
		TODO FIX THIS POINTER ERROR. THe software works with this unfixed, so its ok for now.

		*/
		return 0; //glitch case, not understood pointer error where structure is deleted between loops in basin hopping
	}
	//else std::cout << "did not have a cryptid pointer error" << std::endl;
	if (numBonds == 0) return 1;
	bool** graph = new bool*[natoms];
	for (int i = 0; i < natoms; i++) graph[i] = new bool[natoms];
	int index_a = 0;
	for (int i = 0; i < s->set.size(); i++)
	{
		for (int a = 0; a < s->set[i].size(); a++)
		{
			double radiusA1 = picometersToAngstrom(covalentRadii(s->elements[i], 1));
			double radiusA2 = picometersToAngstrom(covalentRadii(s->elements[i], 2));
			double radiusA3 = picometersToAngstrom(covalentRadii(s->elements[i], 3));
			

			int index_b = 0;
			for (int j = 0; j < s->composition.size(); j++)
			{
				for (int b = 0; b < s->composition[j]; b++)
				{
					if (i == j && b == a) {
						graph[index_a][index_b] = 0;//same atom

					}else{
						double radiusB1 = picometersToAngstrom(covalentRadii(s->elements[j], 1));
						double radiusB2 = picometersToAngstrom(covalentRadii(s->elements[j], 2));
						double radiusB3 = picometersToAngstrom(covalentRadii(s->elements[j], 3));


						double diameterR1 = radiusA1 + radiusB1;
						double boundR1 = diameterR1 * ((double)percent / (double)100);
						double dist = euclideanDistance(s->set[i][a], s->set[j][b]);

						bool bond = false;

						if (dist < diameterR1 + boundR1 && dist > diameterR1 - boundR1)
						{
							bond = true;
						}
						else if (false) {
							double diameterR2 = radiusA2 + radiusB2;
							double boundR2 = diameterR2 * ((double)percent / (double)100);
							if (dist < diameterR2 + boundR2 && dist > diameterR2 - boundR2)
							{
								bond = true;
							}
							else if (radiusB3 != 0 && radiusA3 != 0) {
								double diameterR3 = radiusA3 + radiusB3;
								double boundR3 = diameterR3 * ((double)percent / (double)100);
								if (dist < diameterR3 + boundR3 && dist > diameterR3 - boundR3)
								{
									bond = true;
								}

							}
						}
						graph[index_a][index_b] = bond;
						//if(!bond) std::cout << "a: " << index_a << " b: " << index_b << std::endl;
					}
					index_b++;
				}
			}

			index_a++;
		}
	}
	//step 2: ensure the graph is connected
	bool* visited = new bool[natoms];
	for (int b = 0; b < natoms; b++) visited[b] = false;
	std::queue<int> atomQueue;
	atomQueue.push(0);
	while (!atomQueue.empty())
	{
		int a = atomQueue.front();
		visited[a] = true;
		for (int b = 0; b < natoms; b++)
		{
			if (graph[a][b] && !visited[b]) atomQueue.push(b);
		}


		//remove the current atom
		atomQueue.pop();
	}
	//ensuring all atoms were visited will ensure the graph is connected
	bool connected = true;
	for (int b = 0; b < natoms; b++) if (!visited[b]) {
		connected = false;
		//get information about the disconnection:
		/*std::cout << "atom: " << b << " not connected";
		int tempb = b + 1;
		for (int i = 0; i < s->composition.size(); i++)
		{
			if (tempb > s->composition[i])
			{
				tempb -= s->composition[i];
			}
			else {
				std::cout << ", element: " << s->elements[i] << std::endl;
				i = s->composition.size();
			}
		}
		*/
	}
	//step 3: number of bonds requirement
	bool bondRequirement = true;
	for (int a = 0; a < natoms; a++)
	{
		int bondNo = 0;
		for (int b = 0; b < natoms; b++)
		{
			if (graph[a][b]) bondNo++;
		}
		if (bondNo < numBonds) bondRequirement = false;
	}
	//step 4: cleanup
	for (int i = 0; i < natoms; i++) delete[] graph[i];
	delete[] graph;
	delete[] visited;
	//step 5: return
	//std::cout << connected << bondRequirement << std::endl;
	return (connected && bondRequirement);
}

double MaximizeEuclideanDistance(atom a, atom b,double x, double y, double z) {
	//returns the maximum euclidean distance possible based on possible x y z movement
	double x1 = (a.x - b.x);
	double x2 = 0;
	if (x1 > 0) {
		//we want to stretch a and b as far as possible by making b as small as possible and a as large as possible. 
		//ie 
		//x2 = (a.x + x) - (b.x - x);
		//This is the equivalent of adding 2x.
		x2 = x1 + 2 * x;
	}
	else {
		x2 = x1 - 2 * x;
	}
	double y1 = (a.y - b.y);
	double y2 = 0;
	if (y1 > 0) {
		y2 = y1 + 2 * y;
	}
	else {
		y2 = y1 - 2 * y;
	}
	double z1 = (a.z - b.z);
	double z2 = 0;
	if (z1 > 0) {
		z2 = z1 + 2 * z;
	}
	else {
		z2 = z1 - 2 * z;
	}
	double d = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2);
	return d;
}

double MinimizeEuclideanDistance(atom a, atom b, double x, double y, double z) {
	//returns the minimum euclidean distance possible based on possible x y z movement
	double x1 = (a.x - b.x);
	double x2 = 0;
	if (x1 > 0) {
		//we want to squeeze a and b as close as possible by making b as large as possible and a as small as possible,
		// but we do not want to pass 0 
		//ie 
		//x2 = (a.x - x) - (b.x + x);
		//This is the equivalent of subtracting 2x.
		x2 = x1 - 2 * x;
		if (x2 < 0) x2 = 0;//if we pass 0 go to 0
	}
	else {
		x2 = x1 + 2 * x;
		if (x2 > 0) x2 = 0;
	}
	double y1 = (a.y - b.y);
	double y2 = 0;
	if (y1 > 0) {
		y2 = y1 - 2 * y;
		if (y2 < 0) y2 = 0;
	}
	else {
		y2 = y1 + 2 * y;
		if (y2 > 0) y2 = 0;
	}
	double z1 = (a.z - b.z);
	double z2 = 0;
	if (z1 > 0) {
		z2 = z1 - 2 * z;
		if (z2 < 0) z2 = 0;
	}
	else {
		z2 = z1 + 2 * z;
		if (z2 > 0) z2 = 0;
	}
	double d = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2);
	return d;
}

int radialCriteriaTrapping(structure& s, int percent,double x, double y, double z, int sensitivity)
{
	/*
	This function evaluates if its possible to escape the radial criteria or if you will get trapped based on the range of movement allowed. 
	returns 0 for moveable, returns a new recommended radial criteria if it is not
	*/
	
	bool failure = false;
	int suggestedPercent = 0;
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
						double maxDist = MaximizeEuclideanDistance(s.set[i][a], s.set[j][b],x,y,z);
						double minDist = MinimizeEuclideanDistance(s.set[i][a], s.set[j][b], x, y, z);
						if (minDist < picometersToAngstrom(vanDerWaalsRadii(s.elements[i])))//use minimum as this is a less than
						{
							bool accept = false;
							double bl1 = picometersToAngstrom(covalentRadii(s.elements[i], 1));
							double bl1r = bl1 * ((double)percent / (double)100);
							if ((maxDist > bl1 - bl1r && maxDist < bl1 + bl1r)|| (minDist > bl1 - bl1r && minDist < bl1 + bl1r)) accept = true;//both less than and greater than requires both
							else {
								double bl2 = picometersToAngstrom(covalentRadii(s.elements[i], 2));
								if (bl2 != 0)
								{
									double bl2r = bl2 * ((double)percent / (double)100);
									if ((maxDist > bl2 - bl2r && maxDist < bl2 + bl2r) || (minDist > bl2 - bl2r && minDist < bl2 + bl2r)) accept = true;
									else {
										double bl3 = picometersToAngstrom(covalentRadii(s.elements[i], 3));
										if (bl3 != 0)
										{
											double bl3r = bl3 * ((double)percent / (double)100);
											if ((maxDist > bl3 - bl3r && maxDist < bl3 + bl3r) || (minDist > bl3 - bl3r && minDist < bl3 + bl3r)) accept = true;
										}
									}
								}




							}
							if (!accept) {
								failure = true;
							}
						}
						//else if (bugging) std::cout << "outside of van der waals for " << s.elements[i] << a << " and " << s.elements[j] << b << " being longer than the van der waals radii," << picometersToAngstrom(vanDerWaalsRadii(s.elements[i])) << ", of " << s.elements[i] << a << " at " << dist << std::endl;


					}

				}
			}

		}
	}
	if (failure)
	{
		suggestedPercent = radialCriteriaTrapping(s, percent + sensitivity, x, y, z,sensitivity) + sensitivity;
	}
	else {
		suggestedPercent = 0;
	}

	if (suggestedPercent > 0)
	{
		//std::cout << "WARNING: Structure found which will trap radial criteria with no hope. Percentage temporarily increased to " << suggestedPercent << std::endl;
	}
	//all passed, return true
	return suggestedPercent;


}
structure* seedFromGroupTheory(std::vector<int> composition, std::vector<std::string> elements, double rmax)
{
	//bool remove = true;
	bool symPrint = false;
	//find Cmax
	int Cmax = 0;
	for (int i = 0; i < composition.size(); i++) if (composition[i] > Cmax) Cmax = composition[i];
	if (symPrint)std::cout << "Cmax: " << Cmax<< std::endl;

	//in this smarter way to find just m, n and h
	int n = rand() % Cmax + 1;
	
	if (symPrint)std::cout << "n: " << n << std::endl;
	bool h = false;

	//we need to prevent a symmetry where atoms from multiple elements are forced to be in the center. this will happen when the left over elements that cannot be assigned to Cn symmetry are odd and there is mirror symmetry.
	//but i dont think this is a real problem, so i will not prevent it.
	bool Hallowed = true;
	//figure out if h symmetry is allowed
	//we must ensure only one element has an odd number of atoms to be placed along the central rotational axis
	
	//if(remove) n = 6;

	int nodd = 0;
	for (int i = 0; i < composition.size(); i++)
	{
		int CRAatoms = composition[i] % n;
		if (CRAatoms % 2 == 1) nodd++;
	}
	if (nodd > 0) Hallowed = false;
	if (Hallowed) h = rand() % 2;
	else h = 0;
	//if(remove) if (h == 1) std::cout << "h set to 1" << std::endl;
	//if(remove) std::cout << "h allowed: " << Hallowed << " with n: " << n << std::endl;

	if (symPrint)std::cout << "h: " << h << std::endl;
	int m = n * pow(2, h);
	if (symPrint)std::cout << "m: " << m << std::endl;

	/*
	if (true)
	{
		n = 2;
		h = 0;
		m = 2;
	}
	*/
	/*
	//get m <= Cmax
	int m = rand() % Cmax + 1;
	if (symPrint)std::cout << "m: " << m << std::endl;
	//n must be chosen first as a number that m is divisible by.
	std::vector<int> divisors = {};
	for (int i = 1; i <= m / 2; i++) if (m % i == 0) divisors.push_back(i);
	divisors.push_back(m);
	int n = divisors[rand() % divisors.size()];
	if (symPrint)std::cout << "n: " << n << std::endl;
	//double hv = log2(m / n);
	//if (symPrint)std::cout << "hv: " << hv << std::endl;
	//we do not need hv right now, but if it is to be included, then log2(m/n) must be a whole number, cannot be 1.5 or something. it will not be compatible with symmetry. cannot be cast to the nearest integer either.
	//if m/n = 2, h = 1. otherwise, n must be forced to be m.
	bool h = false;
	//the following conditional is just one way to solve for compatable m, h, and n. really the way m was found above is stupid.
	if (m == n * 2) h = true;
	
	else if (m >= n * 2)
	{
		h = true;
		m = n * 2;
	}
	else {
		n = m;
	}
	
	//if (hv > 0) h = rand() % 2;
	//int v = rand() % (hv - (int)h);
	//for now it doesn't seem like there is a real difference between v and n symmetry, so v will be ingored
	//v = 0;
	if (symPrint)std::cout << "h: " << h << std::endl;
	//construct the partition based on this data
	*/
	int Cangle = (360 /n);
	if (symPrint)std::cout << "Cn angle: " << Cangle << std::endl;
	//the partition is a slice of a sphere, a hemisphere if h. 
	//the position of each atom must be defined with a radius, and two angles.
	std::vector<std::vector<atom>> set;//required at the end to actually produce the structure

	std::vector<std::vector<std::tuple<double, double, double>>> partition;
	std::vector<std::vector<std::tuple<double, double, double>>> partitionMirrorPlane;
	std::vector<std::vector<std::tuple<double, double, double>>> partitionCentralRotationalAxis;
	for (int e = 0; e < composition.size(); e++)
	{
		std::vector<std::tuple<double, double, double>> celement;
		std::vector<std::tuple<double, double, double>> celementMP;
		std::vector<std::tuple<double, double, double>> celementCRA;

		int s = 0;
		if (composition[e] >= m)
		{
			if (symPrint)std::cout << "Composition of "  << elements[e] << " is  above or equal to m"  << std::endl;
			s = composition[e] % m;//these atoms must be along the central rotational axis, as all others are duplicated, even those lying on planes
			if (symPrint)std::cout << "initial s: " << s << std::endl;
			//determine the phi angle, (not the same as Cangle plane)
			for (int a = 0; a < (composition[e] - s) / m; a++) {
				if (symPrint)std::cout << "Cn atoms: " << a << " " << Cangle << std::endl;
				double r = rmax * double(rand()) / double(RAND_MAX);
				double phi = 0;
				double theta = 0;
				if (h) phi = (double)90 * double(rand()) / double(RAND_MAX);
				else phi = (double)180 * double(rand()) / double(RAND_MAX);
				if (phi == 0)
				{
					//these must be moved to s
					s += n * pow(2, h) - 1;//1 is already here
					celementCRA.push_back(std::tuple<double, double, double>(r, phi, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
					if (symPrint)std::cout << "added n to s:  " << s << " " << std::endl;
				}
				else if (phi == 90 && h == 1)
				{
					//there will be n of these, 
					//we can either add n to s or add another set on the plane
					//if (rand() % 2)
					//{
					//	celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 0, Cangle * double(rand()) / double(RAND_MAX)));
					//}
					//else {
					s += n;
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celementMP.push_back(std::tuple<double, double, double>(r, 90, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
					if (symPrint)std::cout << "added n to s:  " << s<< " " << std::endl;
					//}
				}
				else {
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celement.push_back(std::tuple<double, double, double>(r, phi, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
				}

				//place atom

			}
		}
		else s = composition[e];
		if (h)
		{
			//if s >= n we could place some on the mirror plane
			while (s >= n)
			{
				
				s -= n;
				celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 90, (double)Cangle * double(rand()) / double(RAND_MAX)));
				if (symPrint)std::cout << "from s placing n atoms on the mirror plane " << std::endl;
			}
			if (s % 2 == 1)
			{
				//place 1 atom at the center, and place (s-1)/2 along the central rotational axis
				celementCRA.push_back(std::tuple<double, double, double>(0, 0, 0));
				if (symPrint)std::cout << "s%2 == 1; placing atom at the center" << std::endl;
				for (int i = 1; i < s; i++)
				{
					
					double r = rmax * double(rand()) / double(RAND_MAX);
					if (symPrint)std::cout << "s%2 == 1; placing atom at  "<< r << " on CRA" << std::endl;
					celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
				}
			}
			else {
				//place s/2 atoms along the central rotational axis
				for (int i = 0; i < s/2; i++)
				{
					double r = rmax * double(rand()) / double(RAND_MAX);
					if (symPrint)std::cout << "s%2 == 0; placing at omat  " << r << " on CRA" << std::endl;
					celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
				}
				
			}

		}
		else
		{
			//place s atoms along the central rotational axis
			//while (s > n)
			//{
			//	s -= n;
			//	celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 0, Cangle * double(rand()) / double(RAND_MAX)));
			//}

			for (int i = 0; i < s; i++)
			{
				
				double r = (2 * double(rand()) / double(RAND_MAX) - 1) * rmax;
				if (symPrint)std::cout << "in s; placing at omat  " << r << " on CRA" << std::endl;
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}
		}
		partition.push_back(celement);
		partitionMirrorPlane.push_back(celementMP);
		partitionCentralRotationalAxis.push_back(celementCRA);
	}

	/*
	populate
	*/
	for (int e = 0; e < partition.size(); e++)
	{
		std::vector<atom> eSet;
		for (int a = 0; a < partition[e].size(); a++)
		{
			for (int i = 0; i < n; i++)
			{
				double r = std::get<0>(partition[e][a]);
				double theta = M_PI / (double)180 * (std::get<2>(partition[e][a]) + (double)i / (double)n * 360);
				double phi = M_PI / (double)180 * std::get<1>(partition[e][a]);
				double x = r * cos(theta) * sin(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(phi);
				if(symPrint)std::cout << "From Cn symmetry placed atom at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
				z = -z;
				if(h && symPrint)std::cout << "From Cn symmetry (mirror plane) placed atom at " << x << " " << y << " " << z << std::endl;
				if(h) eSet.push_back(atom(x, y, z));//mirror plane duplicate
			}
		}
		for (int a = 0; a < partitionMirrorPlane[e].size(); a++)
		{
			for (int i = 0; i < n; i++)
			{
				double r = std::get<0>(partitionMirrorPlane[e][a]);
				double theta = M_PI / (double)180 * (std::get<2>(partitionMirrorPlane[e][a]) + (double)i / (double)n * 360);
				double phi = 90;
				double x = r * cos(theta);
				double y = r * sin(theta);
				double z = 0;
				if (symPrint)std::cout << "From Cn symmetry placed atom on the mirror plane at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
			}
		}
		for (int a = 0; a < partitionCentralRotationalAxis[e].size(); a++)
		{
				double r = std::get<0>(partitionCentralRotationalAxis[e][a]);
				double x = 0;
				double y = 0;
				double z = r;
				if (symPrint)std::cout << "Placed atom on the central rotational axis at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
				if(h) if(z != 0) eSet.push_back(atom(x, y, -z));
		}
		set.push_back(eSet);
	}


	structure* retValue = new structure(set, elements);
	retValue->nsym = n;
	retValue->hsym = h;
	return retValue;
}

structure* seedWithSymmetry(int n, bool h, std::vector<int> composition, std::vector<std::string> elements, double rmax)
{
	bool possible = true;
	//bool remove = true;
	bool symPrint = false;
	//find Cmax
	int Cmax = 0;
	for (int i = 0; i < composition.size(); i++) if (composition[i] > Cmax) Cmax = composition[i];
	if (n > Cmax) return nullptr;

	//in this smarter way to find just m, n and h
	if (h) {
		int nodd = 0;
		for (int i = 0; i < composition.size(); i++)
		{
			int CRAatoms = composition[i] % n;
			if (CRAatoms % 2 == 1) nodd++;
		}
		if (nodd > 0) return nullptr;
	}
	//if nullptr was returned then the symmetry provided is unacceptable
	if (symPrint)std::cout << "h: " << h << std::endl;
	int m = n * pow(2, h);
	if (symPrint)std::cout << "m: " << m << std::endl;

	/*
	if (true)
	{
		n = 2;
		h = 0;
		m = 2;
	}
	*/
	/*
	//get m <= Cmax
	int m = rand() % Cmax + 1;
	if (symPrint)std::cout << "m: " << m << std::endl;
	//n must be chosen first as a number that m is divisible by.
	std::vector<int> divisors = {};
	for (int i = 1; i <= m / 2; i++) if (m % i == 0) divisors.push_back(i);
	divisors.push_back(m);
	int n = divisors[rand() % divisors.size()];
	if (symPrint)std::cout << "n: " << n << std::endl;
	//double hv = log2(m / n);
	//if (symPrint)std::cout << "hv: " << hv << std::endl;
	//we do not need hv right now, but if it is to be included, then log2(m/n) must be a whole number, cannot be 1.5 or something. it will not be compatible with symmetry. cannot be cast to the nearest integer either.
	//if m/n = 2, h = 1. otherwise, n must be forced to be m.
	bool h = false;
	//the following conditional is just one way to solve for compatable m, h, and n. really the way m was found above is stupid.
	if (m == n * 2) h = true;

	else if (m >= n * 2)
	{
		h = true;
		m = n * 2;
	}
	else {
		n = m;
	}

	//if (hv > 0) h = rand() % 2;
	//int v = rand() % (hv - (int)h);
	//for now it doesn't seem like there is a real difference between v and n symmetry, so v will be ingored
	//v = 0;
	if (symPrint)std::cout << "h: " << h << std::endl;
	//construct the partition based on this data
	*/
	int Cangle = (360 / n);
	if (symPrint)std::cout << "Cn angle: " << Cangle << std::endl;
	//the partition is a slice of a sphere, a hemisphere if h. 
	//the position of each atom must be defined with a radius, and two angles.
	std::vector<std::vector<atom>> set;//required at the end to actually produce the structure

	std::vector<std::vector<std::tuple<double, double, double>>> partition;
	std::vector<std::vector<std::tuple<double, double, double>>> partitionMirrorPlane;
	std::vector<std::vector<std::tuple<double, double, double>>> partitionCentralRotationalAxis;
	for (int e = 0; e < composition.size(); e++)
	{
		std::vector<std::tuple<double, double, double>> celement;
		std::vector<std::tuple<double, double, double>> celementMP;
		std::vector<std::tuple<double, double, double>> celementCRA;

		int s = 0;
		if (composition[e] >= m)
		{
			if (symPrint)std::cout << "Composition of " << elements[e] << " is  above or equal to m" << std::endl;
			s = composition[e] % m;//these atoms must be along the central rotational axis, as all others are duplicated, even those lying on planes
			if (symPrint)std::cout << "initial s: " << s << std::endl;
			//determine the phi angle, (not the same as Cangle plane)
			for (int a = 0; a < (composition[e] - s) / m; a++) {
				if (symPrint)std::cout << "Cn atoms: " << a << " " << Cangle << std::endl;
				double r = rmax * double(rand()) / double(RAND_MAX);
				double phi = 0;
				double theta = 0;
				if (h) phi = (double)90 * double(rand()) / double(RAND_MAX);
				else phi = (double)180 * double(rand()) / double(RAND_MAX);
				if (phi == 0)
				{
					//these must be moved to s
					s += n * pow(2, h) - 1;//1 is already here
					celementCRA.push_back(std::tuple<double, double, double>(r, phi, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
					if (symPrint)std::cout << "added n to s:  " << s << " " << std::endl;
				}
				else if (phi == 90 && h == 1)
				{
					//there will be n of these, 
					//we can either add n to s or add another set on the plane
					//if (rand() % 2)
					//{
					//	celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 0, Cangle * double(rand()) / double(RAND_MAX)));
					//}
					//else {
					s += n;
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celementMP.push_back(std::tuple<double, double, double>(r, 90, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
					if (symPrint)std::cout << "added n to s:  " << s << " " << std::endl;
					//}
				}
				else {
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celement.push_back(std::tuple<double, double, double>(r, phi, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
				}

				//place atom

			}
		}
		else s = composition[e];
		if (h)
		{
			//if s >= n we could place some on the mirror plane
			while (s >= n)
			{

				s -= n;
				celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 90, (double)Cangle * double(rand()) / double(RAND_MAX)));
				if (symPrint)std::cout << "from s placing n atoms on the mirror plane " << std::endl;
			}
			if (s % 2 == 1)
			{
				//place 1 atom at the center, and place (s-1)/2 along the central rotational axis
				celementCRA.push_back(std::tuple<double, double, double>(0, 0, 0));
				if (symPrint)std::cout << "s%2 == 1; placing atom at the center" << std::endl;
				for (int i = 1; i < s; i++)
				{

					double r = rmax * double(rand()) / double(RAND_MAX);
					if (symPrint)std::cout << "s%2 == 1; placing atom at  " << r << " on CRA" << std::endl;
					celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
				}
			}
			else {
				//place s/2 atoms along the central rotational axis
				for (int i = 0; i < s / 2; i++)
				{
					double r = rmax * double(rand()) / double(RAND_MAX);
					if (symPrint)std::cout << "s%2 == 0; placing at omat  " << r << " on CRA" << std::endl;
					celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
				}

			}

		}
		else
		{
			//place s atoms along the central rotational axis
			//while (s > n)
			//{
			//	s -= n;
			//	celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 0, Cangle * double(rand()) / double(RAND_MAX)));
			//}

			for (int i = 0; i < s; i++)
			{

				double r = (2 * double(rand()) / double(RAND_MAX) - 1) * rmax;
				if (symPrint)std::cout << "in s; placing at omat  " << r << " on CRA" << std::endl;
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}
		}
		partition.push_back(celement);
		partitionMirrorPlane.push_back(celementMP);
		partitionCentralRotationalAxis.push_back(celementCRA);
	}

	/*
	populate
	*/
	for (int e = 0; e < partition.size(); e++)
	{
		std::vector<atom> eSet;
		for (int a = 0; a < partition[e].size(); a++)
		{
			for (int i = 0; i < n; i++)
			{
				double r = std::get<0>(partition[e][a]);
				double theta = M_PI / (double)180 * (std::get<2>(partition[e][a]) + (double)i / (double)n * 360);
				double phi = M_PI / (double)180 * std::get<1>(partition[e][a]);
				double x = r * cos(theta) * sin(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(phi);
				if (symPrint)std::cout << "From Cn symmetry placed atom at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
				z = -z;
				if (h && symPrint)std::cout << "From Cn symmetry (mirror plane) placed atom at " << x << " " << y << " " << z << std::endl;
				if (h) eSet.push_back(atom(x, y, z));//mirror plane duplicate
			}
		}
		for (int a = 0; a < partitionMirrorPlane[e].size(); a++)
		{
			for (int i = 0; i < n; i++)
			{
				double r = std::get<0>(partitionMirrorPlane[e][a]);
				double theta = M_PI / (double)180 * (std::get<2>(partitionMirrorPlane[e][a]) + (double)i / (double)n * 360);
				double phi = 90;
				double x = r * cos(theta);
				double y = r * sin(theta);
				double z = 0;
				if (symPrint)std::cout << "From Cn symmetry placed atom on the mirror plane at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
			}
		}
		for (int a = 0; a < partitionCentralRotationalAxis[e].size(); a++)
		{
			double r = std::get<0>(partitionCentralRotationalAxis[e][a]);
			double x = 0;
			double y = 0;
			double z = r;
			if (symPrint)std::cout << "Placed atom on the central rotational axis at " << x << " " << y << " " << z << std::endl;
			eSet.push_back(atom(x, y, z));
			if (h) if (z != 0) eSet.push_back(atom(x, y, -z));
		}
		set.push_back(eSet);
	}


	structure* retValue = new structure(set, elements);
	retValue->nsym = n;
	retValue->hsym = h;
	return retValue;
}
std::vector<std::pair<int, bool>> possibleSymmetry(std::vector<int> composition, std::vector<std::string> elements, double rmax,double rcp,int bondnum, int criteriaIterations)
{
	//bool remove = true;
	bool symPrint = false;
	//find Cmax
	int Cmax = 0;
	for (int i = 0; i < composition.size(); i++) if (composition[i] > Cmax) Cmax = composition[i];
	if (symPrint)std::cout << "Cmax: " << Cmax << std::endl;

	//in this smarter way to find just m, n and h


	std::vector<std::pair<int, bool>> theoreticalSymmetry = {};
	for (int n = Cmax; n > 0; n--)
	{
		theoreticalSymmetry.push_back(std::pair<int, bool>(n, 0));
		bool h = false;

		bool Hallowed = true;
		int nodd = 0;
		for (int i = 0; i < composition.size(); i++)
		{
			int CRAatoms = composition[i] % n;
			if (CRAatoms % 2 == 1) nodd++;
		}
		if (nodd > 0) Hallowed = false;
		if (Hallowed) theoreticalSymmetry.push_back(std::pair<int, bool>(n, 1));
	}

	std::vector<std::pair<int, bool>> reasonableSymmetry = {};


	double rad = rmax;
	srand(time(0));
	for (std::vector<std::pair<int, bool>>::iterator it = theoreticalSymmetry.begin(); it != theoreticalSymmetry.end(); it++)
	{
		int nsym = it->first;
		int hsym = it->second;
		//std::cout << "n,h: " << nsym << " " << hsym << std::endl;
		structure* seed = nullptr;
		int iter = 0;
		bool found = false;
		do {
			if (seed != nullptr) delete seed;
			seed = seedWithSymmetry(nsym, hsym, composition, elements, rmax);
			iter++;

			if (radialCriteria(*seed, rcp))
			{
				if (bondingRequirement(seed, rcp, bondnum)) found = true;
				else found = false;
			}
			else found = false;


		} while (!found && iter < criteriaIterations);
		//std::cout << subit << " " << iter << std::endl;
		if (iter < criteriaIterations) reasonableSymmetry.push_back(std::pair<int,bool>(nsym, hsym));
		delete seed;
	}

	return reasonableSymmetry;
	
}
structure* proceduralSeedFromGroupTheory(int seedNo, std::vector<int> composition, std::vector<std::string> elements, double rmax)
{
	//bool remove = true;
	bool symPrint = false;
	//find Cmax
	int Cmax = 0;
	for (int i = 0; i < composition.size(); i++) if (composition[i] > Cmax) Cmax = composition[i];
	if (symPrint)std::cout << "Cmax: " << Cmax << std::endl;

	//in this smarter way to find just m, n and h
	

	std::vector<std::pair<int, bool>> possibleSymmetry = {};
	for (int n = Cmax; n > 0; n--)
	{
		possibleSymmetry.push_back(std::pair<int, bool>(n, 0));
		bool h = false;

		bool Hallowed = true;
		int nodd = 0;
		for (int i = 0; i < composition.size(); i++)
		{
			int CRAatoms = composition[i] % n;
			if (CRAatoms % 2 == 1) nodd++;
		}
		if (nodd > 0) Hallowed = false;
		if (Hallowed) possibleSymmetry.push_back(std::pair<int, bool>(n, 1));
	}
	int index = seedNo % possibleSymmetry.size();
	int n = possibleSymmetry[index].first;
	bool h = possibleSymmetry[index].second;
	
	return seedWithSymmetry(n, h, composition, elements, rmax);
}

structure* randomStructureFromGroupTheory(std::vector<int> composition, std::vector<std::string> elements, double rmax, int n = 1, bool h = 0, double phimin = 90)
{
	bool symPrint = false;
	
	int m = n * pow(2, h);

	int Cangle = (360 / n);
	if (symPrint)std::cout << "Cn angle: " << Cangle << std::endl;
	std::vector<std::vector<atom>> set;//required at the end to actually produce the structure

	std::vector<std::vector<std::tuple<double, double, double>>> partition;
	std::vector<std::vector<std::tuple<double, double, double>>> partitionMirrorPlane;
	std::vector<std::vector<std::tuple<double, double, double>>> partitionCentralRotationalAxis;
	for (int e = 0; e < composition.size(); e++)
	{
		std::vector<std::tuple<double, double, double>> celement;
		std::vector<std::tuple<double, double, double>> celementMP;
		std::vector<std::tuple<double, double, double>> celementCRA;

		int s = 0;
		if (composition[e] >= m)
		{
			if (symPrint)std::cout << "Composition of " << elements[e] << " is  above or equal to m" << std::endl;
			s = composition[e] % m;//these atoms must be along the central rotational axis, as all others are duplicated, even those lying on planes
			if (symPrint)std::cout << "initial s: " << s << std::endl;
			//determine the phi angle, (not the same as Cangle plane)
			for (int a = 0; a < (composition[e] - s) / m; a++) {
				if (symPrint)std::cout << "Cn atoms: " << a << " " << Cangle << std::endl;
				double r = rmax * double(rand()) / double(RAND_MAX);
				double phi = 0;
				double theta = 0;
				if (h) phi = (double)90 * double(rand()) / double(RAND_MAX);
				else phi = (double)180 * double(rand()) / double(RAND_MAX);
				if (phi < phimin) phi = phimin;
				if (phi == 0)
				{
					//these must be moved to s
					s += n * pow(2, h) - 1;//1 is already here
					celementCRA.push_back(std::tuple<double, double, double>(r, phi, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
					if (symPrint)std::cout << "added n to s:  " << s << " " << std::endl;
				}
				else if (phi == 90 && h == 1)
				{
					//there will be n of these, 
					//we can either add n to s or add another set on the plane
					//if (rand() % 2)
					//{
					//	celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 0, Cangle * double(rand()) / double(RAND_MAX)));
					//}
					//else {
					s += n;
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celementMP.push_back(std::tuple<double, double, double>(r, 90, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
					if (symPrint)std::cout << "added n to s:  " << s << " " << std::endl;
					//}
				}
				else {
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celement.push_back(std::tuple<double, double, double>(r, phi, theta));
					if (symPrint)std::cout << "r phi theta:  " << r << " " << phi << " " << theta << std::endl;
				}

				//place atom

			}
		}
		else s = composition[e];
		if (h)
		{
			//if s > n we could place some on the mirror plane
			while (s > n)
			{

				s -= n;
				celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 90, (double)Cangle * double(rand()) / double(RAND_MAX)));
				if (symPrint)std::cout << "from s placing n atoms on the mirror plane " << std::endl;
			}
			if (s % 2 == 1)
			{
				//place 1 atom at the center, and place (s-1)/2 along the central rotational axis
				celementCRA.push_back(std::tuple<double, double, double>(0, 0, 0));
				if (symPrint)std::cout << "s%2 == 1; placing atom at the center" << std::endl;
				for (int i = 1; i < s; i++)
				{

					double r = rmax * double(rand()) / double(RAND_MAX);
					if (symPrint)std::cout << "s%2 == 1; placing atom at  " << r << " on CRA" << std::endl;
					celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
				}
			}
			else {
				//place s/2 atoms along the central rotational axis
				for (int i = 0; i < s; i++)
				{
					double r = rmax * double(rand()) / double(RAND_MAX);
					if (symPrint)std::cout << "s%2 == 0; placing at omat  " << r << " on CRA" << std::endl;
					celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
				}

			}

		}
		else
		{
			//place s atoms along the central rotational axis
			//while (s > n)
			//{
			//	s -= n;
			//	celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 0, Cangle * double(rand()) / double(RAND_MAX)));
			//}

			for (int i = 0; i < s; i++)
			{

				double r = (2 * double(rand()) / double(RAND_MAX) - 1) * rmax;
				if (symPrint)std::cout << "in s; placing at omat  " << r << " on CRA" << std::endl;
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}
		}
		partition.push_back(celement);
		partitionMirrorPlane.push_back(celementMP);
		partitionCentralRotationalAxis.push_back(celementCRA);
	}

	/*
	populate
	*/
	std::cout << partition.size() << std::endl;
	for (int e = 0; e < partition.size(); e++)
	{
		std::vector<atom> eSet;
		for (int a = 0; a < partition[e].size(); a++)
		{
			for (int i = 0; i < n; i++)
			{
				double r = std::get<0>(partition[e][a]);
				double theta = M_PI / (double)180 * (std::get<2>(partition[e][a]) + (double)i / (double)n * 360);
				double phi = M_PI / (double)180 * std::get<1>(partition[e][a]);
				double x = r * cos(theta) * sin(phi);
				double y = r * sin(theta) * sin(phi);
				double z = r * cos(phi);
				if (symPrint)std::cout << "From Cn symmetry placed atom at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
				z = -z;
				if (h && symPrint)std::cout << "From Cn symmetry (mirror plane) placed atom at " << x << " " << y << " " << z << std::endl;
				if (h) eSet.push_back(atom(x, y, z));//mirror plane duplicate
			}
		}
		for (int a = 0; a < partitionMirrorPlane[e].size(); a++)
		{
			for (int i = 0; i < n; i++)
			{
				double r = std::get<0>(partitionMirrorPlane[e][a]);
				double theta = M_PI / (double)180 * (std::get<2>(partitionMirrorPlane[e][a]) + (double)i / (double)n * 360);
				double phi = 90;
				double x = r * cos(theta);
				double y = r * sin(theta);
				double z = 0;
				if (symPrint)std::cout << "From Cn symmetry placed atom on the mirror plane at " << x << " " << y << " " << z << std::endl;
				eSet.push_back(atom(x, y, z));
			}
		}
		for (int a = 0; a < partitionCentralRotationalAxis[e].size(); a++)
		{
			double r = std::get<0>(partitionCentralRotationalAxis[e][a]);
			double x = 0;
			double y = 0;
			double z = r;
			if (symPrint)std::cout << "Placed atom on the central rotational axis at " << x << " " << y << " " << z << std::endl;
			eSet.push_back(atom(x, y, z));
		}
		set.push_back(eSet);
	}


	structure* retValue = new structure(set, elements);
	retValue->nsym = n;
	retValue->hsym = h;
	return retValue;
}

structure* seedFromPossibleSymmetries(std::vector<std::pair<int, bool>> symmetry, std::vector<int> composition, std::vector<std::string> elements, double rmax, int index )
{
	if (index == -1) index = rand() % symmetry.size();
	else index = index % symmetry.size();
	int n = symmetry[index].first;
	int h = symmetry[index].second;
	return seedWithSymmetry(n, h, composition, elements, rmax);
}



std::vector<std::pair<int,bool>> possiblePartialSymmetry(std::vector<int> composition, std::vector<std::string> elements, double nmultiplyer)
{
	/*
	//many structures are not symmetrical but have aspects of symmetry, such as one element having symmetry but others not.
	//here are the types of partial symmetry this method will account for:

	0. same symmetry for all elements
	1. Individual symmetries for each element all on the same CRA
	2. Individual symmetry by elements, but random placement with respect to eachother
	3. Individual symmetry for 1 or more elements, and random placement for the others

	1. Combine sets of elements before symmetry
	2. 
	
	

	Within these 4 subdivisions it is also possible that an element will resemble the symmetry of a larger molecule with atoms removed.
	Therefore we allow for the n value of Cn symmetry to exceed the number of atoms.
	This is represented by nmultiplyer

	There are an infinite number of these partial symmetries so all cannot be procedurally explored and therefore we seek to randomly generate such a partial symmetry, and then return it so it can be used to generate structures. if no reasonable structure can be found a different partial symmetry can be used
	The following parameters must be used to define a partial symmetry:

	A vector of pairs for each elements n & h    \\ when n & h are 0 this is random
	an integer describing which of the 4 modes will be used
	*/
	std::vector<std::pair<int, bool>> retvalue = {};
	for (int i = 0; i < elements.size(); i++)
	{
		int Cmax = composition[i];
		int n = rand() % (int)(Cmax*nmultiplyer);
		int CRAatoms = composition[i] % n;
		bool h = true;
		if (CRAatoms % 2 == 1) h = false;
		if (h)h = rand() % 2;
		retvalue.push_back(std::pair<int,bool>(n, h));
	}
	return retvalue;
	
}


/*
Partial symmetry set can be used to make the symmetry for an individual element based on provided details
*/
std::vector<atom> partialSymmetrySet(int composition,int n, bool h, double rmax, int holes, bool alignHoles)
{
	int Cangle = (360 / n);

	std::vector<std::tuple<double, double, double>> partition;
	std::vector<std::tuple<double, double, double>> partitionMirrorPlane;
	std::vector<std::tuple<double, double, double>> partitionCentralRotationalAxis;

	//to use the same symmetry procedure as before to populate the mirror plane, partition, and CRA, we must consider a true n and fake n.
	//the fake n is n - holes, what it should appear like to the simulator.
	int nh = n - holes;//if cmax is off, then nh = n
	//now use nh when distributing numbers of atoms, and n when calcualting coordinates


	std::vector<std::tuple<double, double, double>> celement;
	std::vector<std::tuple<double, double, double>> celementMP;
	std::vector<std::tuple<double, double, double>> celementCRA;

	int s = 0;
	//m is related to the distribution of atoms so we use nh
	int m = nh * pow(2, h);


	if (composition >= m)
	{
		s = composition % m;//these atoms must be along the central rotational axis, as all others are duplicated, even those lying on planes
		//determine the phi angle, (not the same as Cangle plane)
		for (int a = 0; a < (composition - s) / m; a++) {
			double r = rmax * double(rand()) / double(RAND_MAX);
			double phi = 0;
			double theta = 0;
			if (h) phi = (double)90 * double(rand()) / double(RAND_MAX);
			else phi = (double)180 * double(rand()) / double(RAND_MAX);
			if (phi == 0)
			{
				//these must be moved to s
				s += nh * pow(2, h) - 1;//1 is already here
				celementCRA.push_back(std::tuple<double, double, double>(r, phi, theta));
			}
			else if (phi == 90 && h == 1)
			{
				s += nh;
				theta = (double)Cangle * double(rand()) / double(RAND_MAX);//cangle has information on n rather than nh, which is good as its for coordinates
				celementMP.push_back(std::tuple<double, double, double>(r, 90, theta));
			}
			else {
				theta = (double)Cangle * double(rand()) / double(RAND_MAX);
				celement.push_back(std::tuple<double, double, double>(r, phi, theta));
			}
		}
	}
	else s = composition;
	if (h)
	{
		//if s >= nh we could place some on the mirror plane
		while (s >= nh)
		{

			s -= nh;
			celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 90, (double)Cangle * double(rand()) / double(RAND_MAX)));
		}
		if (s % 2 == 1)
		{
			celementCRA.push_back(std::tuple<double, double, double>(0, 0, 0));
			for (int i = 0; i < (s-1)/2; i++)
			{

				double r = rmax * double(rand()) / double(RAND_MAX);
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}
		}
		else {
			//place s/2 atoms along the central rotational axis
			for (int i = 0; i < s / 2; i++)
			{
				double r = rmax * double(rand()) / double(RAND_MAX);
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}

		}

	}
	else
	{
		for (int i = 0; i < s; i++)
		{

			double r = (2 * double(rand()) / double(RAND_MAX) - 1) * rmax;
			celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
		}
	}
	partition = celement;
	partitionMirrorPlane = celementMP;
	partitionCentralRotationalAxis = celementCRA;

	std::set<int> partitionSet = {};//the nh of n partitions that will be filled
	std::vector<int> partitions = {};
	if (!alignHoles)
	{
		//if we don't align holes then we need to skip random partitions
		int iterator = 0;
		while (partitionSet.size() != nh)
		{
			if (rand() % n == 0) partitionSet.insert(iterator % nh);
		}
		for (std::set<int>::iterator it = partitionSet.begin(); it != partitionSet.end(); it++) partitions.push_back(*it);
	}
	else {
		for (int i = 0; i < nh; i++) partitions.push_back(i);
	}
	std::vector<atom> eSet;
	for (int a = 0; a < partition.size(); a++)
	{
		for (int ih = 0; ih < nh; ih++)
		{
			int i = partitions[ih];
			double r = std::get<0>(partition[a]);
			double theta = M_PI / (double)180 * (std::get<2>(partition[a]) + (double)i / (double)n * 360);
			double phi = M_PI / (double)180 * std::get<1>(partition[a]);
			double x = r * cos(theta) * sin(phi);
			double y = r * sin(theta) * sin(phi);
			double z = r * cos(phi);
			eSet.push_back(atom(x, y, z));
			z = -z;

			if (h) eSet.push_back(atom(x, y, z));//mirror plane duplicate
		}
	}
	for (int a = 0; a < partitionMirrorPlane.size(); a++)
	{
		for (int ih = 0; ih < nh; ih++)
		{
			int i = partitions[ih];
			double r = std::get<0>(partitionMirrorPlane[a]);
			double theta = M_PI / (double)180 * (std::get<2>(partitionMirrorPlane[a]) + (double)i / (double)n * 360);
			double phi = 90;
			double x = r * cos(theta);
			double y = r * sin(theta);
			double z = 0;
			eSet.push_back(atom(x, y, z));
		}
	}
	for (int a = 0; a < partitionCentralRotationalAxis.size(); a++)
	{
		double r = std::get<0>(partitionCentralRotationalAxis[a]);
		double x = 0;
		double y = 0;
		double z = r;
		eSet.push_back(atom(x, y, z));
		if (h) if (z != 0) eSet.push_back(atom(x, y, -z));
	}
	return eSet;
}
std::vector<atom> partialSymmetrySet(partialSymmetryDescriptor descriptor, double rmax)
{
	int composition = descriptor.composition;
	int n = descriptor.n;
	bool h = descriptor.h;
	int holes = descriptor.holes;
	bool* holeLocations = descriptor.holeLocations;

	int Cangle = (360 / n);

	std::vector<std::tuple<double, double, double>> partition;
	std::vector<std::tuple<double, double, double>> partitionMirrorPlane;
	std::vector<std::tuple<double, double, double>> partitionCentralRotationalAxis;

	//to use the same symmetry procedure as before to populate the mirror plane, partition, and CRA, we must consider a true n and fake n.
	//the fake n is n - holes, what it should appear like to the simulator.
	int nh = n - holes;//if cmax is off, then nh = n
	//now use nh when distributing numbers of atoms, and n when calcualting coordinates


	std::vector<std::tuple<double, double, double>> celement;
	std::vector<std::tuple<double, double, double>> celementMP;
	std::vector<std::tuple<double, double, double>> celementCRA;

	int s = 0;
	//m is related to the distribution of atoms so we use nh
	int m = nh * pow(2, h);

	if (m != 0) {
		if (composition >= m)
		{
			s = composition % m;//these atoms must be along the central rotational axis, as all others are duplicated, even those lying on planes
			//determine the phi angle, (not the same as Cangle plane)
			for (int a = 0; a < (composition - s) / m; a++) {
				double r = rmax * double(rand()) / double(RAND_MAX);
				double phi = 0;
				double theta = 0;
				if (h) phi = (double)90 * double(rand()) / double(RAND_MAX);
				else phi = (double)180 * double(rand()) / double(RAND_MAX);
				if (phi == 0)
				{
					//these must be moved to s
					s += nh * pow(2, h) - 1;//1 is already here
					celementCRA.push_back(std::tuple<double, double, double>(r, phi, theta));
				}
				else if (phi == 90 && h == 1)
				{
					s += nh;
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);//cangle has information on n rather than nh, which is good as its for coordinates
					celementMP.push_back(std::tuple<double, double, double>(r, 90, theta));
				}
				else {
					theta = (double)Cangle * double(rand()) / double(RAND_MAX);
					celement.push_back(std::tuple<double, double, double>(r, phi, theta));
				}
			}
		}
		else s = composition;
	}
	else s = composition;
	if (h)
	{
		//if s >= nh we could place some on the mirror plane
		if (nh != 0) {
			while (s >= nh)
			{

				s -= nh;
				celementMP.push_back(std::tuple<double, double, double>(rmax * double(rand()) / double(RAND_MAX), 90, (double)Cangle * double(rand()) / double(RAND_MAX)));
			}
		}
		if (s % 2 == 1)
		{
			celementCRA.push_back(std::tuple<double, double, double>(0, 0, 0));
			for (int i = 0; i < (s-1)/2; i++)
			{

				double r = rmax * double(rand()) / double(RAND_MAX);
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}
		}
		else {
			//place s/2 atoms along the central rotational axis
			for (int i = 0; i < s / 2; i++)
			{
				double r = rmax * double(rand()) / double(RAND_MAX);
				celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
			}

		}

	}
	else
	{
		for (int i = 0; i < s; i++)
		{

			double r = (2 * double(rand()) / double(RAND_MAX) - 1) * rmax;
			celementCRA.push_back(std::tuple<double, double, double>(r, 0, 0));
		}
	}
	partition = celement;
	partitionMirrorPlane = celementMP;
	partitionCentralRotationalAxis = celementCRA;

	
	
	std::vector<atom> eSet;
	for (int a = 0; a < partition.size(); a++)
	{
		for (int p = 0; p < n; p++)
		{
			if (!holeLocations[p])
			{
			double r = std::get<0>(partition[a]);
			double theta = M_PI / (double)180 * (std::get<2>(partition[a]) + (double)p / (double)n * 360);
			double phi = M_PI / (double)180 * std::get<1>(partition[a]);
			double x = r * cos(theta) * sin(phi);
			double y = r * sin(theta) * sin(phi);
			double z = r * cos(phi);
			eSet.push_back(atom(x, y, z));
			z = -z;

			if (h) eSet.push_back(atom(x, y, z));//mirror plane duplicate

			}
		}
	}
	for (int a = 0; a < partitionMirrorPlane.size(); a++)
	{
		for (int p = 0; p < n; p++)
		{
			if (!holeLocations[p])
			{
				double r = std::get<0>(partitionMirrorPlane[a]);
				double theta = M_PI / (double)180 * (std::get<2>(partitionMirrorPlane[a]) + (double)p / (double)n * 360);
				double phi = 90;
				double x = r * cos(theta);
				double y = r * sin(theta);
				double z = 0;
				eSet.push_back(atom(x, y, z));
			}
		}
	}
	for (int a = 0; a < partitionCentralRotationalAxis.size(); a++)
	{
		double r = std::get<0>(partitionCentralRotationalAxis[a]);
		double x = 0;
		double y = 0;
		double z = r;
		eSet.push_back(atom(x, y, z));
		if (h) if (z != 0) eSet.push_back(atom(x, y, -z));
	}
	return eSet;
}


partialSymmetryDescriptor::partialSymmetryDescriptor(int comp, int n, bool h, int holes, bool* holeLocations) :composition(comp), n(n), h(h), holes(holes), holeLocations(holeLocations)
{

}
partialSymmetryDescriptor::~partialSymmetryDescriptor()
{
	//this might be a bad idea because the hole locations are being made elsewhere now
	//delete[] holeLocations;
}

std::vector<structure*> seedsWithPartialSymmetry(bool extraHoles, bool alignHoles, bool variantSymmetry, double nMultiplier, std::vector<int> composition, std::vector<std::string> elements, double rmax, double radialCriteriaPercent, int bondRequirement = -1, int criteriaIterations = 10000, int nMin = 0)
{
	std::vector<int> nValues = {};
	int Cmax = 0;
	for (int e = 0; e < composition.size(); e++) if (composition[e] > Cmax) Cmax = composition[e];
	int maxN = int(Cmax * nMultiplier);
	for (int n = maxN; n > nMin; n--)
	{
		nValues.push_back(n);
	}
	return seedsWithPartialSymmetry(extraHoles, alignHoles, variantSymmetry, nValues, composition, elements, rmax, radialCriteriaPercent, bondRequirement, criteriaIterations);

}

std::vector<structure*> seedsWithPartialSymmetry(bool extraHoles,bool alignHoles,bool variantSymmetry, std::vector<int> nValues,  std::vector<int> composition, std::vector<std::string> elements, double rmax, double radialCriteriaPercent, int bondRequirement = -1, int criteriaIterations = 10000)
{
	std::vector<structure*> retValue;
	
	std::vector<bool*> pointers;
	/*
	Features to add in the future:
	-Combinations of elements into the same symmetry descriptor
	*/

	int Cmax = 0;
	for (int e = 0; e < composition.size(); e++) if (composition[e] > Cmax) Cmax = composition[e];
	std::vector<std::vector<partialSymmetryDescriptor>> elementDescriptorSets = {};
	for (int e = 0; e < composition.size(); e++) elementDescriptorSets.push_back({});
	for (int it = 0; it < nValues.size(); it++)
	{
		int n = nValues[it];
		bool Hallowed = true;
		int nodd = 0;
		for (int i = 0; i < composition.size(); i++)
		{
			int CRAatoms = composition[i] % n;
			if (CRAatoms % 2 == 1) nodd++;
		}
		if (nodd > 0) Hallowed = false;

		
		
		//create a set for each element
		for (int e = 0; e < composition.size(); e++)
		{
			//the minimum number of holes must be n - composition.size() but can be higher if atoms are put on the CRA
			int holes = 0;
			if (n > composition[e]) holes += n - composition[e];
			std::vector<int> holeValues = {0,holes};
			if (extraHoles) for (int i = 1; i <= composition[e]; i++) holeValues.push_back(holes + i);

			//std::vector<partialSymmetryDescriptor> descriptorSets = {};

			for (int i = 0; i < holeValues.size(); i++) {

				/*
				Todo, unfortunately we actually need to save all these parameters so that these sets can be continuously regenerated later :(
				
				*/
				std::vector<bool*> holeLocations = {};
				if (e == 0) {
					if (!alignHoles)holeLocations = enumeration(n, holeValues[i], true);
					else holeLocations = alignedEnumeration(n, holeValues[i], true);
				}
				else {
					if (!alignHoles)holeLocations = enumeration(n, holeValues[i], false);
					else holeLocations = alignedEnumeration(n, holeValues[i], false);
				}
				std::cout << "hole locations size for n: " << n << ", holes: "<< holeValues[i] << ":" << holeLocations.size() << std::endl;
				for (std::vector<bool*>::iterator pit = holeLocations.begin(); pit != holeLocations.end(); pit++) pointers.push_back(*pit);
				if (Hallowed)
				{
					if(!variantSymmetry && e == 0){
						//if all components have the same symmetry then there is no reason to rotate the holes in the first element
						elementDescriptorSets[e].push_back(partialSymmetryDescriptor(composition[e], n, 1, holeValues[i], holeLocations[0]));
					}
					else {
						for (std::vector<bool*>::iterator it = holeLocations.begin(); it != holeLocations.end(); it++)
						{
							elementDescriptorSets[e].push_back(partialSymmetryDescriptor(composition[e], n, 1, holeValues[i], *it));
						}
					}
				}
				if (!variantSymmetry && e == 0) {
					//as above if all components have the same symmetry then there is no reason to rotate the holes in the first element
					elementDescriptorSets[e].push_back(partialSymmetryDescriptor(composition[e], n, 0, holeValues[i], holeLocations[0]));
				}
				else {
					for (std::vector<bool*>::iterator it = holeLocations.begin(); it != holeLocations.end(); it++)
					{
						elementDescriptorSets[e].push_back(partialSymmetryDescriptor(composition[e], n, 0, holeValues[i], *it));
					}
				}
			}
			//elementDescriptorSets[e].push_back(descriptorSets);
		}
		/*
		COMPOUNDS
		
		todo, for now this is unnecessary.
		A future feature which will combine and divide elements into symmetries
		*/
		if (!variantSymmetry) {
			int numStructures = 1;
			for (int i = 0; i < elementDescriptorSets.size(); i++) {
				numStructures *= elementDescriptorSets[i].size();
			}
			std::cout << "Number of structure types from partial symmetry " << numStructures << std::endl;

			for (int i = 0; i < numStructures; i++)
			{
				if (i % 25 == 0) std::cout << "calculated so far: " << i << ", structures found: " << retValue.size() << std::endl;
				//i = 84241;
				//std::cout << "WPI 1" << std::endl;
				structure* testStruc = nullptr;
				int iter = 0;
				bool found = false;
				//std::cout << "0 " << std::endl;
				std::vector<int> iterators = {};
				for (int e = 0; e < composition.size(); e++) iterators.push_back(0);
				//std::cout << "WPI 2" << std::endl;
				//std::cout << "1 " << std::endl;
				for (int it = 0; it < i; it++)
				{
					//std::cout << "WPI 3 " << it << std::endl;
					int e = 0;
					iterators[e] += 1;
					//std::cout << "it " << it << std::endl;
					while (iterators[e] >= elementDescriptorSets[e].size())
					{
						//std::cout << "moving up at " << it << std::endl;
						iterators[e] = 0;
						e++;
						iterators[e] += 1;
					}

				}



				//std::cout << "2 " << std::endl;
				//std::cout << "current iterators: ";
				//for (int it = 0; it < iterators.size(); it++) std::cout << iterators[it] << " ";
				//std::cout << std::endl;
				//std::cout << "starting do loop for " << i << std::endl;
				do {
					std::vector<std::vector<atom>> structureSets;

					//std::cout << "a" << std::endl;
					for (int e = 0; e < composition.size(); e++)
					{
						//std::cout << "loading PSD for " << elements[e] << std::endl;
						//elementDescriptorSets[e][iterators[e]].print();
						/*
						the setOfE is nonsense
						. we need to check that the correct amount of atoms is being chosen
						*/
						std::vector<atom> setOfE = partialSymmetrySet(elementDescriptorSets[e][iterators[e]],rmax);
						//std::cout << "for " << elements[e] << ":" << setOfE.size();
						structureSets.push_back(setOfE);
					}
					//std::cout << "b" << std::endl;
					if (testStruc != nullptr) delete testStruc;
					/*
					What is being passed to structure is nonsense.
					*/

					testStruc = new structure(structureSets, elements);
					//std::cout << "c" << std::endl;
					//std::cout << "radial criteria check " << std::endl;
					//std::cout << "structure sets size: " << structureSets.size() << std::endl;
					//testStruc->print();
					if (radialCriteria(*testStruc, radialCriteriaPercent))
					{
						//std::cout << "d" << std::endl;
						if (bondingRequirement(testStruc, radialCriteriaPercent, bondRequirement)) found = true;
						else found = false;
					}
					else found = false;
					//std::cout << "e" << std::endl;

					iter++;
				} while (!found && iter < criteriaIterations);

				//std::cout << "ending do loop for " << i << std::endl;

				//std::cout << "structure search while loop over with " << iter << std::endl;
				if (iter < criteriaIterations) {
					retValue.push_back(testStruc);
					//std::cout << "structure found" << std::endl;
					//testStruc->print();
				}
				else {
					delete testStruc;
					//std::cout << "symmetry failed moving on" << std::endl;
				}

			}
			elementDescriptorSets = {};
			for (int e = 0; e < composition.size(); e++) elementDescriptorSets.push_back({});
		}



		
	}


	if (variantSymmetry) {
		int numStructures = 1;
		for (int i = 0; i < elementDescriptorSets.size(); i++) {
			numStructures *= elementDescriptorSets[i].size();
		}
		std::cout << "Number of structure types from partial symmetry " << numStructures << std::endl;

		for (int i = 0; i < numStructures; i++)
		{
			//if (i % 1000 == 0) std::cout << "calculated so far: " << i << std::endl;
			//i = 84241;
			//std::cout << "WPI 1" << std::endl;
			structure* testStruc = nullptr;
			int iter = 0;
			bool found = false;
			//std::cout << "0 " << std::endl;
			std::vector<int> iterators = {};
			for (int e = 0; e < composition.size(); e++) iterators.push_back(0);
			//std::cout << "WPI 2" << std::endl;
			//std::cout << "1 " << std::endl;
			for (int it = 0; it < i; it++)
			{
				//std::cout << "WPI 3 " << it << std::endl;
				int e = 0;
				iterators[e] += 1;
				//std::cout << "it " << it << std::endl;
				while (iterators[e] >= elementDescriptorSets[e].size())
				{
					//std::cout << "moving up at " << it << std::endl;
					iterators[e] = 0;
					e++;
					iterators[e] += 1;
				}

			}



			//std::cout << "2 " << std::endl;
			//std::cout << "current iterators: ";
			//for (int it = 0; it < iterators.size(); it++) std::cout << iterators[it] << " ";
			//std::cout << std::endl;
			//std::cout << "starting do loop for " << i << std::endl;
			do {
				std::vector<std::vector<atom>> structureSets;

				//std::cout << "a" << std::endl;
				for (int e = 0; e < composition.size(); e++)
				{
					//std::cout << "loading PSD for " << elements[e] << std::endl;
					//elementDescriptorSets[e][iterators[e]].print();
					structureSets.push_back(partialSymmetrySet(elementDescriptorSets[e][iterators[e]], rmax));
				}
				//std::cout << "b" << std::endl;
				if (testStruc != nullptr) delete testStruc;
				testStruc = new structure(structureSets, elements);
				//std::cout << "c" << std::endl;
				//std::cout << "radial criteria check " << std::endl;
				//std::cout << "structure sets size: " << structureSets.size() << std::endl;
				//testStruc->print();
				if (radialCriteria(*testStruc, radialCriteriaPercent))
				{
					//std::cout << "d" << std::endl;
					if (bondingRequirement(testStruc, radialCriteriaPercent, bondRequirement)) found = true;
					else found = false;
				}
				else found = false;
				//std::cout << "e" << std::endl;

				iter++;
			} while (!found && iter < criteriaIterations);

			//std::cout << "ending do loop for " << i << std::endl;

			//std::cout << "structure search while loop over with " << iter << std::endl;
			if (iter < criteriaIterations) {
				retValue.push_back(testStruc);
				//testStruc->print();
			}
			else delete testStruc;

		}
	}

	/*
	COMBINES
	
	*/
	for (std::vector<bool*>::iterator pit = pointers.begin(); pit != pointers.end(); pit++) delete[] * pit;
	//std::cout << "internal pss size: " << retValue.size() << std::endl;
	return retValue;
}
std::vector<structure*> seedsWithSpecificSymmetry(bool extraHoles, bool alignHoles,int n, bool h, std::vector<int> composition, std::vector<std::string> elements, double rmax, double radialCriteriaPercent, int bondRequirement = -1, int criteriaIterations = 10000)
{
	std::vector<structure*> retValue;

	std::vector<bool*> pointers;
	/*
	Features to add in the future:
	-Combinations of elements into the same symmetry descriptor
	*/

	int Cmax = 0;
	for (int e = 0; e < composition.size(); e++) if (composition[e] > Cmax) Cmax = composition[e];
	std::vector<std::vector<partialSymmetryDescriptor>> elementDescriptorSets = {};
	for (int e = 0; e < composition.size(); e++) elementDescriptorSets.push_back({});
	
		

		//create a set for each element
		for (int e = 0; e < composition.size(); e++)
		{
			//the minimum number of holes must be n - composition.size() but can be higher if atoms are put on the CRA
			int holes = 0;
			if (n > composition[e]) holes += n - composition[e];
			std::vector<int> holeValues = { 0,holes };
			if (extraHoles) for (int i = 1; i <= composition[e]; i++) holeValues.push_back(holes + i);

			//std::vector<partialSymmetryDescriptor> descriptorSets = {};

			for (int i = 0; i < holeValues.size(); i++) {

				/*
				Todo, unfortunately we actually need to save all these parameters so that these sets can be continuously regenerated later :(

				*/
				std::vector<bool*> holeLocations = {};
				if (e == 0){
					if (!alignHoles)holeLocations = enumeration(n, holeValues[i],true);
					else holeLocations = alignedEnumeration(n, holeValues[i],true);
				}else{
					if (!alignHoles)holeLocations = enumeration(n, holeValues[i],false);
					else holeLocations = alignedEnumeration(n, holeValues[i],false);
				}
				std::cout << "hole locations size for n: " << n << ", holes: " << holeValues[i] << ":" << holeLocations.size() << std::endl;
				for (std::vector<bool*>::iterator pit = holeLocations.begin(); pit != holeLocations.end(); pit++) pointers.push_back(*pit);

					for (std::vector<bool*>::iterator it = holeLocations.begin(); it != holeLocations.end(); it++)
					{
						elementDescriptorSets[e].push_back(partialSymmetryDescriptor(composition[e], n, h, holeValues[i], *it));
					}
			}
			//elementDescriptorSets[e].push_back(descriptorSets);
		}
		/*
		COMPOUNDS

		todo, for now this is unnecessary.
		A future feature which will combine and divide elements into symmetries
		*/

			int numStructures = 1;
			for (int i = 0; i < elementDescriptorSets.size(); i++) {
				numStructures *= elementDescriptorSets[i].size();
			}
			std::cout << "Number of structure types from partial symmetry " << numStructures << std::endl;

			for (int i = 0; i < numStructures; i++)
			{
				if (i % 1000 == 0) std::cout << "calculated so far: " << i << std::endl;
				//i = 84241;
				//std::cout << "WPI 1" << std::endl;
				structure* testStruc = nullptr;
				int iter = 0;
				bool found = false;
				//std::cout << "0 " << std::endl;
				std::vector<int> iterators = {};
				for (int e = 0; e < composition.size(); e++) iterators.push_back(0);
				//std::cout << "WPI 2" << std::endl;
				//std::cout << "1 " << std::endl;
				for (int it = 0; it < i; it++)
				{
					//std::cout << "WPI 3 " << it << std::endl;
					int e = 0;
					iterators[e] += 1;
					//std::cout << "it " << it << std::endl;
					while (iterators[e] >= elementDescriptorSets[e].size())
					{
						//std::cout << "moving up at " << it << std::endl;
						iterators[e] = 0;
						e++;
						iterators[e] += 1;
					}

				}



				//std::cout << "2 " << std::endl;
				//std::cout << "current iterators: ";
				//for (int it = 0; it < iterators.size(); it++) std::cout << iterators[it] << " ";
				//std::cout << std::endl;
				//std::cout << "starting do loop for " << i << std::endl;
				do {
					std::vector<std::vector<atom>> structureSets;

					//std::cout << "a" << std::endl;
					for (int e = 0; e < composition.size(); e++)
					{
						//std::cout << "loading PSD for " << elements[e] << std::endl;
						//elementDescriptorSets[e][iterators[e]].print();
						/*
						the setOfE is nonsense
						. we need to check that the correct amount of atoms is being chosen
						*/
						std::vector<atom> setOfE = partialSymmetrySet(elementDescriptorSets[e][iterators[e]], rmax);
						//std::cout << "for " << elements[e] << ":" << setOfE.size();
						structureSets.push_back(setOfE);
					}
					//std::cout << "b" << std::endl;
					if (testStruc != nullptr) delete testStruc;
					/*
					What is being passed to structure is nonsense.
					*/

					testStruc = new structure(structureSets, elements);
					//std::cout << "c" << std::endl;
					//std::cout << "radial criteria check " << std::endl;
					//std::cout << "structure sets size: " << structureSets.size() << std::endl;
					//testStruc->print();
					if (radialCriteria(*testStruc, radialCriteriaPercent))
					{
						//std::cout << "d" << std::endl;
						if (bondingRequirement(testStruc, radialCriteriaPercent, bondRequirement)) found = true;
						else found = false;
					}
					else found = false;
					//std::cout << "e" << std::endl;
					//found = true;
					iter++;
				} while (!found && iter < criteriaIterations);

				//std::cout << "ending do loop for " << i << std::endl;

				//std::cout << "structure search while loop over with " << iter << std::endl;
				if (iter < criteriaIterations) {
					retValue.push_back(testStruc);
					//std::cout << "structure found" << std::endl;
					//testStruc->print();
				}
				else {
					delete testStruc;
					//std::cout << "symmetry failed moving on" << std::endl;
				}

			}
			elementDescriptorSets = {};
			for (int e = 0; e < composition.size(); e++) elementDescriptorSets.push_back({});

	/*
	COMBINES

	*/
	for (std::vector<bool*>::iterator pit = pointers.begin(); pit != pointers.end(); pit++) delete[] * pit;
	//std::cout << "internal pss size: " << retValue.size() << std::endl;
	return retValue;
}


structure* getSeedFromSymmetryMode(char mode, std::vector<int> composition, std::vector<std::string> elements, double rmax,int call, std::vector<std::pair<int, bool>> symmetries )
{
	structure* retValue = nullptr;
	switch (mode)
	{
	case 'i':
		//i for iterate through the possible symmetries
		//i for iterate through the possible symmetries
		retValue = seedFromPossibleSymmetries(symmetries, composition, elements, rmax, call);
		break;
	case 'r':
		//r for randomly choose from the possible symmetries
		retValue = seedFromPossibleSymmetries(symmetries, composition, elements, rmax, -1);
		break;
	case 'g':
		retValue = seedFromGroupTheory(composition, elements, rmax);
		//symetrical seed with random symmetry. g for group theory seed
		break;
	case 'n':
		retValue = randomSeedInRadius(rmax, composition, elements);
		//std::cout << "attempting to make seed" << std::endl;
		//retValue = randomSeedInRadiusFromRadialCriteria(rmax, composition, elements, 20);
		//std::cout << "seedmade" << std::endl;
		//retValue=	new structure(rmax, rmax, rmax, composition, elements);//causes cryptid pointer errors
		//std::cout << "PRINTING SEED" << std::endl;
		//retValue->print();
		//no symmetry
		break;
	}
	return retValue;
}

structure* getSymmetricaSeed(char mode, std::vector<int> composition, std::vector<std::string> elements, double rmax, int call = 0, std::vector<std::pair<int, bool>> symmetries = { std::pair<int,bool>(1,0) },int bondRequirement = 1, double radialCriteriaPercent = 20, int maxBIterations = 100, int maxRCIterations = 10000, bool Cmax= true, double nmultiplier = 1)
{
	structure* s = nullptr;

	if (mode != 'a')
	{
		/*
		The following 4 modes are accounted for in getSeedFromSymmetryMode:
		//procedural froms set
		//random symmetry from set
		//symetrical seed with random symmetry
		//no symmetry
		*/
		/*
		TODO: there may be something wrong with the following loops, i think it is possible to get stuck in an infinite loop if radial criteria is passing while bonding requirement is not but this can only happen on iteration 1. flip from do while to normal while to fix.
		*/
		do {
			//s = seedFromGroupTheory(composition, elements, rlimit);
			s = getSeedFromSymmetryMode(mode, composition, elements, rmax, call, symmetries);
			while (!radialCriteria(*s, radialCriteriaPercent))
			{
				delete s;
				s = getSeedFromSymmetryMode(mode, composition, elements, rmax, call, symmetries);
			}
		} while (!bondingRequirement(s, radialCriteriaPercent, bondRequirement));
	}else{
		//todo
		/*
		TODO
		*/
	}
	return s;
}

std::vector<bool*> enumeration(int size, int totalValue, bool initial)
{
	//base cases
	if (totalValue > size) return {};
	if (size == 1) {
		bool* b1 = new bool[1];
		if (totalValue == 1) b1[0] = 1;
		else b1[0] = 0;
		return { b1 };
	}
	if (totalValue == 0)
	{
		bool* b1 = new bool[size];
		for (int i = 0; i < size; i++) b1[i] = 0;
		return { b1 };
	}
	std::vector<bool*> retVal;

	if (initial)
	{
		std::vector<bool*> uninitialized = enumeration(size - 1, totalValue - 1, false);
		for (std::vector<bool*>::iterator it = uninitialized.begin(); it != uninitialized.end(); it++)
		{
			bool* newBool = new bool[size];
			for (int i = 0; i < size - 1; i++) newBool[i] = (*it)[i];
			delete[](*it);
			newBool[size - 1] = 1;
			retVal.push_back(newBool);
		}

	}
	else {

		//other cases
		std::vector<bool*> oneLess = enumeration(size - 1, totalValue - 1, false);
		std::vector<bool*> zeroLess = enumeration(size - 1, totalValue, false);


		for (std::vector<bool*>::iterator it = oneLess.begin(); it != oneLess.end(); it++)
		{
			bool* newBool = new bool[size];
			for (int i = 0; i < size - 1; i++) newBool[i] = (*it)[i];
			delete[](*it);
			newBool[size - 1] = 1;
			retVal.push_back(newBool);
		}
		for (std::vector<bool*>::iterator it = zeroLess.begin(); it != zeroLess.end(); it++)
		{
			bool* newBool = new bool[size];
			for (int i = 0; i < size - 1; i++) newBool[i] = (*it)[i];
			delete[](*it);
			newBool[size - 1] = 0;
			retVal.push_back(newBool);
		}
	}

	return retVal;
}

std::vector<bool*> alignedEnumeration(int size, int totalValue, bool initial)
{
	if (totalValue == 0)
	{
		bool* newSet = new bool[size];
		for (int i = 0; i < size; i++)
		{
			newSet[i] = 0;
		}
		return {newSet};
	}
	//this one just returns them all in a circle only
	

	std::vector<bool*> retVal;
	for (int i = 0; i < size; i++)
	{
		bool* newSet = new bool[size];
		for (int it = 0; it < size; it++) newSet[it] = 0;
		for (int j = 0; j < totalValue; j++) {
			if (i + j >= size) newSet[i + j - size] = 1;
			else newSet[i + j] = 1;
		}
		retVal.push_back(newSet);
		if (initial) i = size;
	}

	return retVal;
}


void partialSymmetryDescriptor::print()
{
	std::cout << "PSD: composition:" << composition << ", n:" << n << ", h:" << h << ", holes:" << holes << std::endl;
}


void writeToXyz(std::vector<structure*> structures,std::string outputName)
{
	std::ofstream outFile(outputName);
	if (!outFile.is_open()) {
		std::cout << "Error: Unable to open file " << outputName << " for writing." << std::endl;
		return;
	}

	for (int s= 0; s < structures.size();s++)
	{
		structure* cS = structures[s];
		int natoms = 0;
		//find center of mass
		double tx = 0;
		double ty = 0;
		double tz = 0;
		for (int i = 0; i < cS->set.size(); i++)
		{
			for (int j = 0; j < cS->set[i].size(); j++)
			{
				tx += cS->set[i][j].x;
				ty += cS->set[i][j].y;
				tz += cS->set[i][j].z;
				natoms++;
			}
		}
		tx = tx / natoms;
		ty = ty / natoms;
		tz = tz / natoms;

		/*

		*/

		//outFile << "# Structure OBJ File" << std::endl;

		outFile << std::to_string(natoms) << "\n" << "Frame " << std::to_string(s) <<  " ,Energy: " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << cS->energy << " , time: " << cS->time << "\n";
		for (int i = 0; i < cS->set.size(); i++)
		{
			for (int j = 0; j < cS->set[i].size(); j++)
			{
				outFile << cS->elements[i] << " " << std::to_string(cS->set[i][j].x - tx) << " " << std::to_string(cS->set[i][j].y - ty) << " " << std::to_string(cS->set[i][j].z - tz) << "\n";
			}
		}
	}
	
	

	

	outFile.close();


}

std::vector<structure*> getSeeds(std::vector<int> composition, std::vector<std::string> elements, double radius, double radialCriteriaPercent, int bondRequirement,int criteriaIterations,bool partialSymmetry, bool extraHoles, bool alignHoles, bool variantSymmetry, double nMultiplier, int nMin,std::vector<int> nValues)
{
	std::vector<structure*> seeds;
	std::vector<structure*> symmetrySeeds = {};
	std::vector<structure*> partialSymmetrySeeds = {};
	/*
	Pure symmetry
	*/
	int Cmax = 0;
	for (int i = 0; i < composition.size(); i++) if (composition[i] > Cmax) Cmax = composition[i];
	//in this smarter way to find just m, n and h


	std::vector<std::pair<int, bool>> theoreticalSymmetry = {};
	for (int n = Cmax; n > 0; n--)
	{
		theoreticalSymmetry.push_back(std::pair<int, bool>(n, 0));
		bool h = false;

		bool Hallowed = true;
		int nodd = 0;
		for (int i = 0; i < composition.size(); i++)
		{
			int CRAatoms = composition[i] % n;
			if (CRAatoms % 2 == 1) nodd++;
		}
		if (nodd > 0) Hallowed = false;
		if (Hallowed) theoreticalSymmetry.push_back(std::pair<int, bool>(n, 1));
	}

	std::vector<std::pair<int, bool>> reasonableSymmetry = {};

	for (std::vector<std::pair<int, bool>>::iterator it = theoreticalSymmetry.begin(); it != theoreticalSymmetry.end(); it++)
	{
		int nsym = it->first;
		int hsym = it->second;
		//std::cout << "n,h: " << nsym << " " << hsym << std::endl;
		structure* seed = nullptr;
		int iter = 0;
		bool found = false;
		do {
			if (seed != nullptr) delete seed;
			seed = seedWithSymmetry(nsym, hsym, composition, elements,radius);
			iter++;

			if (radialCriteria(*seed, radialCriteriaPercent))
			{
				if (bondingRequirement(seed, radialCriteriaPercent, bondRequirement)) found = true;
				else found = false;
			}
			else found = false;
			

		} while (!found && iter < criteriaIterations);
		//std::cout << subit << " " << iter << std::endl;
		if (iter < criteriaIterations) {
			symmetrySeeds.push_back(seed);
			//std::cout << "pushed " << nsym << " " << hsym << " with iter: " << iter << " and c: " << criteriaIterations << std::endl;
		}
		else delete seed;
	}

	/*
	Pure symmetry ove
	*/
	if (partialSymmetry) {
		if(nValues.size() == 0) partialSymmetrySeeds = seedsWithPartialSymmetry(extraHoles, alignHoles, variantSymmetry, nMultiplier, composition, elements, radius, radialCriteriaPercent, bondRequirement, criteriaIterations, nMin);
		else partialSymmetrySeeds = seedsWithPartialSymmetry(extraHoles, alignHoles, variantSymmetry, nValues, composition, elements, radius, radialCriteriaPercent, bondRequirement, criteriaIterations);
	}
	/*
	Add a no symmetry mode?
	*/

	for (std::vector<structure*>::iterator it = symmetrySeeds.begin(); it != symmetrySeeds.end(); it++) seeds.push_back(*it);
	for (std::vector<structure*>::iterator it = partialSymmetrySeeds.begin(); it != partialSymmetrySeeds.end(); it++) seeds.push_back(*it);
	return seeds;
}

bool structure::operator<(const structure& other)
{
	return this->energy < other.energy;
}


std::vector<structure*> structuresFromXYZ(std::string fileName)
{
	std::vector<structure*> retval;
	std::ifstream file(fileName);
	std::string in;
	while (std::getline(file, in))
	{
		int natoms = std::stoi(in);
		if (natoms > 0)
		{
			std::getline(file, in);//frame number
			//try read out the energy
			int spaces = 0;
			std::string energy = "";
			for (int i = 0; i < in.size(); i++)
			{
				if (in[i] == ' ') spaces++;
				else {
					if (spaces == 3) energy += in[i];
				}
			}
			double nrg = DBL_MAX;
			std::cout << energy << std::endl;
			if(energy.size() > 0)nrg = std::stod(energy);
			std::vector<std::pair<std::string, atom>> readAtoms;
			for (int i = 0; i < natoms; i++)
			{
				std::getline(file, in);
				//std::cout << in << std::endl;
				int endOfElement = 0;
				int endOfX = 0;
				int endOfY = 0;
				int endOfZ = 0;
				int endOfLine = 0;
				int whitespace = 0;
				for (int c = 0; c < in.size(); c++)
				{
					if (in[c] == ' ')
					{
						bool skipc = false;
						if (c != in.size() - 1) {
							if (in[c + 1] == ' ') {
								skipc = true;//more whitespace
								whitespace++;
							}
						}
						if (!skipc) {
							if (endOfElement > 0)
							{
								if (endOfX > 0)
								{
									if (endOfY > 0)
									{
										endOfZ = c;
									}
									else endOfY = c;
								}
								else endOfX = c;
							}
							else endOfElement = c - whitespace;
						}
					}
				}
				endOfZ = in.size();
				std::string element = in.substr(0, endOfElement);
				double x = std::stod(in.substr(endOfElement + 1, endOfX - endOfElement));
				double y = std::stod(in.substr(endOfX + 1, endOfY - endOfX));
				double z = std::stod(in.substr(endOfY +1, endOfZ - endOfY));
				std::pair<std::string, atom> eap = std::pair<std::string, atom>(element, atom(x, y, z));
				readAtoms.push_back(eap);
			}
			std::set<std::string> elementSet = {};
			for (std::vector<std::pair<std::string, atom>>::iterator it = readAtoms.begin(); it != readAtoms.end(); it++)
			{
				elementSet.insert(it->first);
			}
			std::vector<std::string> elements = {};
			std::vector<int> composition = {};
			std::vector<std::vector<atom>> atomSets = {};
			for (std::set<std::string>::iterator it = elementSet.begin(); it != elementSet.end(); it++)
			{
				elements.push_back(*it);
				composition.push_back(0);
				atomSets.push_back({});
			}
			for (std::vector<std::pair<std::string, atom>>::iterator it = readAtoms.begin(); it != readAtoms.end(); it++)
			{
				int e = 0;
				for (int i = 0; i < elements.size(); i++)
				{
					if (it->first == elements[i])
					{
						composition[i]++;
						atomSets[i].push_back(it->second);
					}
				}
			}
			structure* struc = new structure(atomSets, elements);
			struc->energy = nrg;
			retval.push_back(struc);
		}
	}
	return retval;

}

bool energyCompare(structure* a, structure* b) { return a->energy < b->energy; };


structure* randomSeedInRadius(double radius, std::vector<int> composition, std::vector<std::string> elements)
{
	std::vector<std::vector<atom>> set;
	for (int e = 0; e < composition.size(); e++)
	{
		std::vector<atom > element;
		for (int a = 0; a < composition[e]; a++)
		{
			double r = radius*double(rand()) / double(RAND_MAX);
			double phi = 360.0 * double(rand()) / double(RAND_MAX);
			double theta = 360.0 * double(rand()) / double(RAND_MAX);
			double x = radius * sin(theta) * cos(phi);
			double y = radius * sin(theta) * sin(phi);
			double z = radius * cos(theta);
			element.push_back(atom(x, y, z));
		}
		set.push_back(element);
	}
	structure* retValue = new structure(set, elements);
	return retValue;
}

structure* randomSeedInRadiusFromRadialCriteria(double radius, std::vector<int> composition, std::vector<std::string> elements,double radialCriteriaPercent)
{
	std::vector<std::vector<atom>> set;
	std::vector<atom> allatoms;
	for (int e = 0; e < composition.size(); e++)
	{
		std::vector<atom > element;
		for (int a = 0; a < composition[e]; a++)
		{
			bool totalaccept = false;
			atom catom;
			while (!totalaccept) {
				totalaccept = true;
				double r = radius * double(rand()) / double(RAND_MAX);
				double phi = 360.0 * double(rand()) / double(RAND_MAX);
				double theta = 360.0 * double(rand()) / double(RAND_MAX);
				double x = radius * sin(theta) * cos(phi);
				double y = radius * sin(theta) * sin(phi);
				double z = radius * cos(theta);

				/*
				Check that the atom will fit in
				*/
				catom = atom(x, y, z);
				
				int compositionIterator = 0;
				int atomIterator = 0;
				std::string celement = elements[0];
				for (int b = 0; b < allatoms.size(); b++)
				{
					

					double diameterVDW = picometersToAngstrom(vanDerWaalsRadii(elements[e])) + picometersToAngstrom(vanDerWaalsRadii(celement));
					double diameterR1 = picometersToAngstrom(covalentRadii(elements[e], 1)) + picometersToAngstrom(covalentRadii(celement, 1));
					double diameterR2 = picometersToAngstrom(covalentRadii(elements[e], 2)) + picometersToAngstrom(covalentRadii(celement, 2));
					if (picometersToAngstrom(covalentRadii(elements[e], 2)) == 0 || picometersToAngstrom(covalentRadii(celement, 2)) == 0) diameterR2 = 0;
					double diameterR3 = picometersToAngstrom(covalentRadii(elements[e], 3)) + picometersToAngstrom(covalentRadii(celement, 3));
					if (picometersToAngstrom(covalentRadii(elements[e], 3)) == 0 || picometersToAngstrom(covalentRadii(celement, 3)) == 0) diameterR2 = 0;
					double boundR1 = diameterR1 * ((double)radialCriteriaPercent / (double)100);
					double boundR2 = diameterR2 * ((double)radialCriteriaPercent / (double)100);
					double boundR3 = diameterR3 * ((double)radialCriteriaPercent / (double)100);
					//do not perform on the same atom
					double dist = euclideanDistance(allatoms[b], catom);
					bool accept = false;

					if (dist < diameterVDW)
					{

						if (dist > diameterR1 - boundR1 && dist < diameterR1 + boundR1) accept = true;
						else {
							if (diameterR2 != 0)
							{
								if (dist > diameterR2 - boundR2 && dist < diameterR2 + boundR2) accept = true;
								else {
									if (diameterR3 != 0)
									{
										if (dist > diameterR3 - boundR3 && dist < diameterR3 + boundR3) accept = true;
									}
								}
							}
						}
					}
					if (!accept) totalaccept = false;
					//else if (bugging) std::cout << "outside of van der waals for " << s.elements[i] << a << " and " << s.elements[j] << b << " being longer than the van der waals radii," << picometersToAngstrom(vanDerWaalsRadii(s.elements[i])) << ", of " << s.elements[i] << a << " at " << dist << std::endl;


					atomIterator++;
					if (atomIterator >= composition[compositionIterator])
					{
						atomIterator = 0;
						compositionIterator++;
						celement = elements[compositionIterator];
					}

				}
			}

			element.push_back(catom);
			allatoms.push_back(catom);
		}
		set.push_back(element);
	}
	structure* retValue = new structure(set, elements);
	return retValue;
}

std::vector<int *> multinomialCoefficient(int n, std::vector<int> k)
{
	
		//base cases
	if (n == 1)
	{
		int remainder = 0;
		for (int i = 0; i < k.size(); i++) if (k[i] == 1) remainder = i;
		int* arrangement = new int[1];
		arrangement[0] = remainder;
		return {arrangement};
	}
	std::vector<int*> retVal = {};
	for (int ki = 0; ki < k.size(); ki++)
	{
		if (k[ki] > 0)
		{
			std::vector<int> kprev = k;
			kprev[ki]--;
			std::vector<int*> prev = multinomialCoefficient(n - 1, kprev);
			for (std::vector<int*>::iterator pit = prev.begin(); pit != prev.end(); pit++)
			{
				int* arrangement = new int[n];
				for (int i = 0; i < n - 1; i++) arrangement[i] = (*pit)[i];
				delete[](*pit);
				arrangement[n - 1] = ki;
				retVal.push_back(arrangement);
			}
		}
	}
	return retVal;


}

bool connectednessRequirement(structure* s, int percent, bool numBonds)
{
	//unlike bonding requirement, this ensures atoms are within covalent criterias upper bound, not lower bound

	int natoms = 0;
	for (int i = 0; i < s->set.size(); i++) natoms += s->set[i].size();
	//std::cout << "natoms: " << natoms << std::endl;
	if (natoms == 0) {
		std::cout << "cryptid pointer error detected from createXnew loop in BHC algo" << std::endl;
		return 0; 
	}
	if (natoms == 1) return 1;
	if (numBonds == 0) {
		return 1;
		
	}
	bool** graph = new bool* [natoms];
	for (int i = 0; i < natoms; i++) graph[i] = new bool[natoms];
	int index_a = 0;
	for (int i = 0; i < s->set.size(); i++)
	{
		for (int a = 0; a < s->set[i].size(); a++)
		{
			double radiusA1 = picometersToAngstrom(covalentRadii(s->elements[i], 1));
			int index_b = 0;
			for (int j = 0; j < s->set.size(); j++)
			{
				for (int b = 0; b < s->set[j].size(); b++)
				{
					if (i == j && b == a) {
						graph[index_a][index_b] = 0;//same atom
					}
					else {
						double radiusB1 = picometersToAngstrom(covalentRadii(s->elements[j], 1));
						

						double diameterR1 = radiusA1 + radiusB1;
						double boundR1 = diameterR1 * ((double)percent / (double)100);
						double dist = euclideanDistance(s->set[i][a], s->set[j][b]);

						bool bond = false;

						if (dist < diameterR1 + boundR1)
						{
							//std::cout << "dist: " << dist << " diameterR1: " << diameterR1 << " boindR1: " << boundR1 << std::endl;
							bond = true;
						}
						
						graph[index_a][index_b] = bond;
						//if(!bond) std::cout << "a: " << index_a << " b: " << index_b << std::endl;
					}
					index_b++;
				}
			}

			index_a++;
		}
	}
	//step 2: ensure the graph is connected
	bool* visited = new bool[natoms];
	for (int b = 0; b < natoms; b++) visited[b] = false;
	std::queue<int> atomQueue;
	atomQueue.push(0);
	while (!atomQueue.empty())
	{
		int a = atomQueue.front();
		visited[a] = true;
		//std::cout << "visited" << std::endl;
		for (int b = 0; b < natoms; b++)
		{
			if (graph[a][b] && !visited[b]) atomQueue.push(b);
		}


		//remove the current atom
		atomQueue.pop();
	}
	//ensuring all atoms were visited will ensure the graph is connected
	bool connected = true;
	for (int b = 0; b < natoms; b++) if (!visited[b]) {
		//std::cout << "proof of check 1" << std::endl;
		connected = false;
		//get information about the disconnection:
		/*std::cout << "atom: " << b << " not connected";
		int tempb = b + 1;
		for (int i = 0; i < s->composition.size(); i++)
		{
			if (tempb > s->composition[i])
			{
				tempb -= s->composition[i];
			}
			else {
				std::cout << ", element: " << s->elements[i] << std::endl;
				i = s->composition.size();
			}
		}
		*/
	}
	//step 3: number of bonds requirement
	bool bondRequirement = true;
	for (int a = 0; a < natoms; a++)
	{
		//std::cout << "proof of check 2" << std::endl;
		int bondNo = 0;
		for (int b = 0; b < natoms; b++)
		{
			if (graph[a][b]) bondNo++;
		}
		if (bondNo < numBonds) bondRequirement = false;
	}
	//step 4: cleanup
	for (int i = 0; i < natoms; i++) delete[] graph[i];
	delete[] graph;
	delete[] visited;
	//step 5: return
	//std::cout << connected << bondRequirement << std::endl;
	return (connected && bondRequirement);
}

std::vector<bool> randomNchooseC(int n, int c)
{
	//std::cout << "parameters given: " << n << " " << c << std::endl;
	std::vector<bool>  retVal = {};
	int rn = n;
	int rc = c;
	for (int i = 0; i < n; i++)
	{
		retVal.push_back(false);
		//std::cout << i << std::endl;
		int rval = rand() % rn;
		if (rval < rc) {
			rc--;
			retVal[i] = true;
		}
		rn--;
		
	}
	return retVal;
}
template<class T>
std::vector<std::vector<T>> combineSets(std::vector<std::vector<T>> setA, std::vector<std::vector<T>> setB)
{
	std::vector<std::vector<T>> retVal = {};
	for (int i = 0; i < setA.size(); i++)
	{
		//std::cout << "CS--" << i << std::endl;
		std::vector<T> element = {};
		for (int j = 0; j < setA[i].size(); j++) element.push_back(setA[i][j]);
		//std::cout << "CS-2-" << i << std::endl;
		for (int j = 0; j < setB[i].size(); j++) element.push_back(setB[i][j]);
		retVal.push_back(element);
	}
	return retVal;
}
void moveSet(std::vector<std::vector<atom>>& set, double x, double y, double z)
{
	for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++) for (std::vector<atom>::iterator itt = it->begin(); itt != it->end(); itt++)
	{
		itt->x += x;
		itt->y += y;
		itt->z += z;
	}
}

std::vector<structure*> proceduralSeeds(std::vector<int> composition, std::vector<std::string> elements, int numRings, std::vector<int> nValues, double covalentCriteriaPercent,std::vector<int> allowedHoles,int keepPerCombo, int keepPerRingCombination, double achange)
{
	//std::cout << "elements : ";
	//for (int eee = 0; eee < elements.size(); eee++) std::cout << elements[eee] << " ";
	//std::cout << std::endl;
	/*
	It is reasonable to expected to use this function for multiple ring values, but not acceptable to use it for the same ring value with different nValues as this will lead to duplicate seeds.
	nvalues cannot contain duplicates
	*/

	
	bool dOut = true;//debug output
	int debugging = 0;
	//int debu = true;
	int natoms = 0;
	for (int i = 0; i < composition.size(); i++) natoms += composition[i];
	//if (maxHoles == INT_MAX) maxHoles = INT_MAX - natoms;//avoid integer overflow
	/*
	Determine n values for the rings, this will also determine the number of holes.
	*/
	std::vector < std::vector<int>> ringCombinations = {};//the inner vector is length numRings, and the integer is the iterator of nValues for the nvalue of each ring. the total vector contains the number of possible combinations
	int nrc = pow(nValues.size(),numRings);//number of combinations of ring number distributions for each element
	//this variable will be useful later but only needs to be found once
	double totalCovalentRadiiHorizontal = 0;
	double totalCovalentRadiiVertical = 0;
	for (int i = 0; i < elements.size(); i++) totalCovalentRadiiHorizontal += double(covalentRadii(elements[i], 1,true))/100.0 * composition[i];//divide by 100 for pm to A
	double averageCovalentRadiiHorizontal = totalCovalentRadiiHorizontal / (double)natoms;
	//std::cout << averageCovalentRadiiHorizontal << " acr,n " << natoms <<" "<<  totalCovalentRadiiHorizontal << elements[0] << covalentRadii(elements[0], 1,true) << std::endl;
	for (int i = 0; i < elements.size(); i++) totalCovalentRadiiVertical += double(covalentRadii(elements[i], 1,false)) / 100.0 * composition[i];//divide by 100 for pm to A
	double averageCovalentRadiiVertical = totalCovalentRadiiVertical / (double)natoms;
	//std::cout << averageCovalentRadiiVertical << " acr,n " << natoms << " " << totalCovalentRadiiVertical << elements[0] << covalentRadii(elements[0], 1,false) << std::endl;
	std::vector<int> ringCombination = {};
	for (int rc = 0; rc < numRings; rc++) ringCombination.push_back(0);
	std::cout << "num rings: " << numRings << " ring combo size: " << ringCombination.size();
	for (int rc = 0; rc < nrc; rc++) {
		int ring = 0;
		if(rc != 0)ringCombination[ring]++;
		while (ringCombination[ring] >= nValues.size())
		{
			ringCombination[ring] = 0;
			ringCombination[++ring]++;
			if (debugging > 1) std::cout << "d ringCombination[" << ring << "]: " << ringCombination[ring] << std::endl;
		}
		//Check if this ring combination is valid
		int maxAtoms = 0;
		for (std::vector<int>::iterator it = ringCombination.begin(); it != ringCombination.end(); it++ ) maxAtoms += nValues[(*it)];//the maximum number of atoms that can be fit in this combination is the sum of the nvalues of the ring.
		if (allowedHoles.size() == 0) ringCombinations.push_back(ringCombination);
		else {
			bool RCallowed = false;
			for (std::vector<int>::iterator it = allowedHoles.begin(); it != allowedHoles.end(); it++) if (maxAtoms >= natoms && maxAtoms - *it == natoms) {
				RCallowed = true;
				std::cout <<"P: " << maxAtoms << ":maxAtoms,natoms:" << natoms << ",numHoles:" << *it << std::endl;
			}
			else {
				std::cout << "DNP: " << maxAtoms << ":maxAtoms,natoms:" << natoms << ",numHoles:" << *it << std::endl;
			}
			if(RCallowed) ringCombinations.push_back(ringCombination);
		}
		//if (maxAtoms >= natoms && maxAtoms - maxHoles <= natoms) //utilizes the copy constructor to get this instanc of ring combination, without the further changes
		//if (maxAtoms - maxHoles <= natoms)std::cout << "MA: " << maxAtoms << ", mH: " << maxHoles << ", natoms: " << natoms << std::endl;
		//else std::cout << "passed" << std::endl;
	}
	/*
	Replace the indeces in ring combinations with the number of atoms in the rings
	*/
	for (int rco = 0; rco < ringCombinations.size(); rco++) for (int rcoi = 0; rcoi < ringCombinations[rco].size(); rcoi++) ringCombinations[rco][rcoi] = nValues[ringCombinations[rco][rcoi]];
	if (debugging || true)
	{
		std::cout << "ring combinations:" << ringCombinations.size() << std::endl;
		for (std::vector<std::vector<int>>::iterator ccit = ringCombinations.begin(); ccit != ringCombinations.end(); ccit++)
		{
			for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	/*
	remove ring combinations that are exact reverses
	*/
	std::vector < std::vector<int>> filteredRingCombinations = {};
	bool* keep = new bool[ringCombinations.size()];
	for (int rC1 = 0; rC1 < ringCombinations.size(); rC1++)
	{
		keep[rC1] = true;
		for (int rC2 = rC1; rC2 < ringCombinations.size(); rC2++)
		{
			if (rC1 != rC2)
			{
				bool identical = true;
				for (int r = 0; r < numRings; r++)
				{
					if (ringCombinations[rC1][r] != ringCombinations[rC1][numRings - (r + 1)]) identical = false;
				}
				if (identical) keep[rC1] = false;
			}
		}
	}
	for (int rC1 = 0; rC1 < ringCombinations.size(); rC1++) if (keep[rC1]) filteredRingCombinations.push_back(ringCombinations[rC1]);
	
	ringCombinations = filteredRingCombinations;
	//filter done

	/*
				For each ring with smaller rings as the next ring we need combinations for combining these rings.
				For example for the rings of size 3,2,1,4 we want a sitation with 3+2, 1, 4 and 3+2+1, 4, and 3, 2+1, 4
				we do not want to combine rings if smaller rings proceed because 1,2,3,4 is considered a different combination, and 3+2+1 will be duplicate.
				Additionally in 3, 1, 2, 4, 3, 2+1, 4 would be duplicate
				and 1+2+3+4 would be duplucate from 4,3,2,1 and 3,2,1,4
	*/
	std::vector < std::vector< std::vector<int>>> combinedRingCombinations = {};
	if (debugging)
	{
		std::cout << "fitlered ring combinations:" << filteredRingCombinations.size() << std::endl;
		for (std::vector<std::vector<int>>::iterator ccit = filteredRingCombinations.begin(); ccit != filteredRingCombinations.end(); ccit++)
		{
			for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	for (std::vector<std::vector<int>>::iterator ringCombination = ringCombinations.begin(); ringCombination != ringCombinations.end(); ringCombination++)
	{
		std::vector<std::vector<std::vector<int>>> currentRingCombination = {};//for each ring, this contains the possibility of smaller rings
		for (int ring = 0; ring < ringCombination->size(); ring++)
		{
			std::vector<std::vector<int>> currentRingCombinations = {};//each of these represents one ring, the vector inside represents combinations of this ring with smaller rings
			std::vector<int> thisRing = { (*ringCombination)[ring] };
			currentRingCombinations.push_back(thisRing);
			for (int nextRing = ring+1; nextRing < ringCombination->size(); nextRing++)
			{
				if ((*ringCombination)[nextRing-1] > (*ringCombination)[nextRing])
				{
					//create a combination with this smaller ring inside and the rings inbetween
					thisRing.push_back((*ringCombination)[nextRing]);
					currentRingCombinations.push_back(thisRing);
				}
				else {
					nextRing = ringCombination->size();
				}

			}
			currentRingCombination.push_back(currentRingCombinations);
		}
		//somethign to completely enumerate all fused ring combos that are compatible. ie 4+3 , 2+1 and not 4+3+2, 2+1

		/*
		Now we have a vector (for each ring) containing a vector of all their possible sub rings represented as a vector
		Now we must turn these into possible ring combinations which is just a vector of a ring which is a vector
		*/
		std::vector<std::vector<std::vector<int>>> validRingCombinations = { };
		if (debugging > 1)
		{
			std::cout << "currentRingCombination: " << currentRingCombination.size() << std::endl;
			for (std::vector < std::vector< std::vector<int>>>::iterator cit = currentRingCombination.begin(); cit != currentRingCombination.end(); cit++)
			{
				for (std::vector<std::vector<int>>::iterator ccit = cit->begin(); ccit != cit->end(); ccit++)
				{
					std::cout << "{";
					for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
					std::cout << "} ";
				}
				std::cout << std::endl;
			}
		}

		for (int ring = 0; ring < currentRingCombination.size(); ring++)
		{
			std::vector<std::vector<std::vector<int>>> currentLevelValidRingCombinations = { };
			if (ring == 0)
			{
				for (int ringCombo = 0; ringCombo < currentRingCombination[ring].size(); ringCombo++) currentLevelValidRingCombinations.push_back({ currentRingCombination[ring][ringCombo] });
				if (debugging) {
					std::cout << "current level valid ring combinations size: " << currentLevelValidRingCombinations.size() << std::endl;
					for (std::vector < std::vector< std::vector<int>>>::iterator cit = currentLevelValidRingCombinations.begin(); cit != currentLevelValidRingCombinations.end(); cit++)
					{
						for (std::vector<std::vector<int>>::iterator ccit = cit->begin(); ccit != cit->end(); ccit++)
						{
							std::cout << "{";
							for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
							std::cout << "} ";
						}
						std::cout << std::endl;
					}
				}
				validRingCombinations = currentLevelValidRingCombinations;
				if (debugging) std::cout << "moving on from ring 0, valid ring combinations size: " << validRingCombinations.size() << std::endl;
			}
			else {
				if (debugging) std::cout << "not ring 0, valid ring combinations before adding next ring: " << validRingCombinations.size() << std::endl;
				//for each previous valid ring combination
				for (int pvrc = 0; pvrc < validRingCombinations.size(); pvrc++)
				{

					//get the current number of rings inside
					int includedRings = 0;
					for (int pvrci = 0; pvrci < validRingCombinations[pvrc].size(); pvrci++) includedRings += validRingCombinations[pvrc][pvrci].size();
					if (debugging) std::cout << "included rings[index]:" << includedRings << std::endl;
					if (includedRings-1 < ring)
					{
						
						for (int ringCombo = 0; ringCombo < currentRingCombination[ring].size(); ringCombo++)
						{
							std::vector<std::vector<int>> currentLevelValidRingCombination = validRingCombinations[pvrc];
							currentLevelValidRingCombination.push_back(currentRingCombination[ring][ringCombo]);
							currentLevelValidRingCombinations.push_back(currentLevelValidRingCombination);
						}
					}
					else {
						//this ring is already included, dont add a new ring
						std::vector<std::vector<int>> currentLevelValidRingCombination = validRingCombinations[pvrc];
						currentLevelValidRingCombinations.push_back(currentLevelValidRingCombination);
					}
				}
				if (debugging) {
					std::cout << "current level valid ring combinations size: " << currentLevelValidRingCombinations.size() << std::endl;
					for (std::vector < std::vector< std::vector<int>>>::iterator cit = currentLevelValidRingCombinations.begin(); cit != currentLevelValidRingCombinations.end(); cit++)
					{
						for (std::vector<std::vector<int>>::iterator ccit = cit->begin(); ccit != cit->end(); ccit++)
						{
							std::cout << "{";
							for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
							std::cout << "} ";
						}
						std::cout << std::endl;
					}
				}
			}
			validRingCombinations = currentLevelValidRingCombinations;
		}
		for (std::vector<std::vector<std::vector<int>>>::iterator vit = validRingCombinations.begin(); vit != validRingCombinations.end(); vit++) combinedRingCombinations.push_back((*vit));
	}


	if (dOut) std::cout << "# of ring combinations: " << combinedRingCombinations.size() << std::endl;
	if (dOut)
	{
		for (std::vector < std::vector< std::vector<int>>>::iterator cit = combinedRingCombinations.begin(); cit != combinedRingCombinations.end(); cit++)
		{
			for (std::vector<std::vector<int>>::iterator ccit = cit->begin(); ccit != cit->end(); ccit++)
			{
				std::cout << "{";
				for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
				std::cout << "} ";
			}
			std::cout << std::endl;
		}
	}
	
	std::vector<std::vector<std::pair<std::vector<std::vector<atom>>,std::vector<std::vector<int>>>>> allSets;
	
	/*
	Tabulate total number of combinations
	*/
	int totalCombos = 0;
	for (std::vector<std::vector<std::vector<int>>>::iterator ringCombination = combinedRingCombinations.begin(); ringCombination != combinedRingCombinations.end(); ringCombination++)
	{
		int ringSpots = 0;
		for (std::vector<std::vector<int>>::iterator it = ringCombination->begin(); it != ringCombination->end(); it++) for (std::vector<int>::iterator rit = it->begin(); rit != it->end(); rit++) ringSpots += (*rit);//the maximum number of atoms that can be fit in this combination is the sum of the nvalues of the ring.
		int nHoles = ringSpots - natoms;
		int nOnes = 0;
		std::vector<bool*> holeLocations = enumeration(ringSpots - nOnes, natoms - nOnes, true);//the boolean at the end indicates if one location will always be held constant. this should be done to reduce hole locations that are identical
		std::vector<int*> elementDistributions = multinomialCoefficient(natoms, composition);
		for (std::vector<bool*>::iterator holeLocation = holeLocations.begin(); holeLocation != holeLocations.end(); holeLocation++)
		{
			for (std::vector<int*>::iterator elementDistribution = elementDistributions.begin(); elementDistribution != elementDistributions.end(); elementDistribution++)
			{
				totalCombos++;
			}
		}
	}
	std::cout << "Total combos: " << totalCombos << " reducing to at most " << keepPerCombo << std::endl;
	std::vector<bool> keepVector;
	if(keepPerCombo > 0) keepVector = randomNchooseC(totalCombos, keepPerCombo);
	else keepVector = randomNchooseC(totalCombos, totalCombos);
	int comboIterator = 0;

	for (std::vector<std::vector<std::vector<int>>>::iterator ringCombination = combinedRingCombinations.begin(); ringCombination != combinedRingCombinations.end(); ringCombination++)
	{
		std::vector<std::pair<std::vector<std::vector<atom>>,std::vector<std::vector<int>>>> notAllSets;
		/*if (debu) {
			std::cout << "D0 ring combination";
			for (std::vector<std::vector<int>>::iterator ccit = ringCombination->begin(); ccit != ringCombination->end(); ccit++)
			{
				std::cout << "{";
				for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
				std::cout << "} ";
			}
			std::cout << std::endl;

		}*/
		/*//for each ring combination now we need to iterate combinations of possible holes.
		we must be careful to not create duplicate structures, or produce as few as possible. 
		Therefore because in later steps it will be possible to place rings inside of eachother, if the ring size is 1 a hole will not be allocated.

		*/
		int ringSpots = 0;
		for (std::vector<std::vector<int>>::iterator it = ringCombination->begin(); it != ringCombination->end(); it++) for(std::vector<int>::iterator rit = it->begin(); rit!= it->end(); rit++) ringSpots += (*rit);//the maximum number of atoms that can be fit in this combination is the sum of the nvalues of the ring.
		int nHoles = ringSpots - natoms;
		//if(true)std::cout << "nHoles: " << nHoles << std::endl;
		int nOnes = 0;//the number of rings with n = 1 , these cannot be holes!
		/*
		{1} cannot be a hole, but can {x,1} have a hole at 1? It is a goode idea to allow for {x,1} to have a hole at the one spot, allowing for wider rings which are unrealistic on their own but realistic if there are rings above and below
		I changed my mind, no ones will have holes
		*/
		for (std::vector<std::vector<int>>::iterator it = ringCombination->begin(); it != ringCombination->end(); it++) for(std::vector<int>::iterator oit = it->begin(); oit != it->end(); oit++) if(*oit == 1) nOnes++;

		std::vector<bool*> holeLocations = enumeration(ringSpots - nOnes, natoms - nOnes, true);//the boolean at the end indicates if one location will always be held constant. this should be done to reduce hole locations that are identical
		//now that we have the hole locations we also want the distribution of elements in each spot
		std::vector<int*> elementDistributions = multinomialCoefficient(natoms, composition);
		
		//combinations of hole locations and element distributions
		if (dOut) std::cout << "\t# of hole location and element distribution combinations: " << holeLocations.size()*elementDistributions.size() << std::endl;

		
		

		for (std::vector<bool*>::iterator holeLocation = holeLocations.begin(); holeLocation != holeLocations.end(); holeLocation++)
		{
			for (std::vector<int*>::iterator elementDistribution = elementDistributions.begin(); elementDistribution != elementDistributions.end(); elementDistribution++)
			{
				if (keepVector[comboIterator])
				{


					bool atomPrintout = false;
					if (false)
					{
						//print details to see if the mistake has already been made
						std::cout << "printout\n---------------------\nhole locations: ";
						for (int i = 0; i < ringSpots - nOnes; i++) std::cout << (*holeLocation)[i] << " ";
						std::cout << "std::endl";
						std::cout << "EDS: " << elementDistributions.size() << std::endl;
						std::cout << "element distribution:";
						for (int i = 0; i < natoms; i++) std::cout << (*elementDistribution)[i] << " ";
						std::cout << "std::endl";
						std::cout << "ring combinaiton: ";
						for (std::vector<std::vector<int>>::iterator ccit = ringCombination->begin(); ccit != ringCombination->end(); ccit++)
						{
							std::cout << "{";
							for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
							std::cout << "} ";
						}
						std::cout << std::endl;
						std::cout << "--------------------\n";
					}

					//if (debu) std::cout << "D1: inner loop" << std::endl;
					//std::cout << "creating new set" << std::endl;
					std::vector<std::vector<atom>> set;
					std::vector<std::vector<int>> assignment;
					int assign = 0;

					for (int e = 0; e < elements.size(); e++) {
						set.push_back({});
						assignment.push_back({});
					}
					if (set.size() == 0) std::cout << "WARNING SET SIZE 0 AFTER CREATION" << std::endl;



					if (achange > 0)
					{

						/*
					Finally we can create the structure!
					*/
					/*
					For each ring we must determine a center radii (for each atom type) that allows atoms in the ring to be spaced out well.
					If there are center atoms, then we will choose whichever radii is larger, the one from the distance from the center atom or between atoms.
					As a ring may have many types of atoms, each will be farther from the center and this will be determined for each element type individually as if they were the only one in the ring.
					*/
					//std::cout << "RM-- 1" << std::endl;
						std::vector<std::vector<atom>> growingSet;
						std::vector<std::vector<int>> growingAssignment;
						int gassign = 0;
						for (int e = 0; e < elements.size(); e++) {
							growingSet.push_back({});
							growingAssignment.push_back({});
						}
						//std::cout << "RM-- 2" << std::endl;
						int atomsPrior = 0;//counts the atoms we already visited, for access from element distribution
						int atomsPriorHoles = 0;//counts atoms we alreadt visited - ones, to access hole distribution
						int atomIterator = 0;
						int atomIteratorHoles = 0;
						double lastZ = -achange;//whats the best starting value?

						bool validCombination = true;

						//std::cout << "RM-- 3" << std::endl;
						for (int ring = 0; ring < ringCombination->size() && validCombination; ring++)
						{

							//std::cout << "RM-- 4:"<< ring << std::endl;
							std::vector<std::vector<atom>> ringOfAtoms;
							std::vector<std::vector<int>> ringAssignment;
							int rassign = 0;
							for (int e = 0; e < elements.size(); e++) {
								ringOfAtoms.push_back({});
								ringAssignment.push_back({});
							}
							double lastRingAvgRadii = 0;
							//std::cout << "RM-- 5" << std::endl;
							double lastRadius = 0;
							for (int subring = (*ringCombination)[ring].size() - 1; subring >= 0 && validCombination; subring--)
							{
								//std::cout << "RM-- 6:" << subring << std::endl;
								double radius = lastRadius;
								bool subRingApproved = false;
								int tatomsPrior = atomsPrior;//counts the atoms we already visited, for access from element distribution
								int tatomsPriorHoles = atomsPriorHoles;//counts atoms we alreadt visited - ones, to access hole distribution
								int tatomIterator = atomIterator;
								int tatomIteratorHoles = atomIteratorHoles;

								int radiusSteps = 0;
								int maxRadiusSteps = (*ringCombination)[ring][subring] * (double)5 / achange;//not sure about this for large rings

								

								double minRadius = 0;
								double maxRadius = 0;
								bool minFound = false;
								bool maxFound = false;

								

								while ((!minFound || !maxFound) && validCombination) {
									int snatoms = 0;
									tatomsPrior = atomsPrior;//counts the atoms we already visited, for access from element distribution
									tatomsPriorHoles = atomsPriorHoles;//counts atoms we alreadt visited - ones, to access hole distribution
									tatomIterator = atomIterator;
									tatomIteratorHoles = atomIteratorHoles;
									//std::cout << "stuck R? ?" << radius << std::endl;

									//std::cout << "RM-- 7 looping for subring approval" << std::endl;
									std::vector<std::vector<atom>> subringOfAtoms;
									std::vector<std::vector<int>> subRingAssignment;
									int sassign = assign;//this used to be 0, but that reset the assignment every ring. we need sassign as we may destroy this structure many times and dont want assignment to increase forever.
									
									for (int e = 0; e < elements.size(); e++) {
										subringOfAtoms.push_back({});
										subRingAssignment.push_back({});
									}
									//std::cout << "RM-- 8" << std::endl;
									for (int a = 0; a < (*ringCombination)[ring][subring]; a++)
									{
										//std::cout << "RM-- 9:" << a << std::endl;
										bool thisHole = true;
										if ((*ringCombination)[ring][subring] != 1) {
											if (atomPrintout) std::cout << "H" << tatomIteratorHoles << ',' << (*holeLocation)[tatomIteratorHoles] << ':';
											thisHole = (*holeLocation)[tatomIteratorHoles];//somehow 1 is for no hole now
											tatomIteratorHoles++;
										}

										if (thisHole)
										{
											snatoms++;
											int thisElement = (*elementDistribution)[tatomIterator];

											tatomIterator++;

											double angle = a * 2 * M_PI / (*ringCombination)[ring][subring];

											double x = radius * cos(angle);
											double y = radius * sin(angle);
											if (atomPrintout) std::cout << elements[thisElement] << " ";

											if ((*ringCombination)[ring][subring] == 1)
											{
												x = 0;
												y = 0;
											}

											if (elements[thisElement].size() > 2)
											{

												std::vector<structure* > ssss = structuresFromXYZ(elements[thisElement]);
												moveSet(ssss[0]->set, x, y, 0);
												for (int es = 0; es < ssss[0]->set.size(); es++)
												{
													std::string celement = ssss[0]->elements[es];
													int thisElementE = 0;
													for (int ee = 0; ee < elements.size(); ee++) if (celement == elements[ee]) {
														thisElementE = ee;
													}
													for (int as = 0; as < ssss[0]->set[es].size(); as++)
													{
														ssss[0]->set[es][as].fixedGroup = thisElement;
														subringOfAtoms[thisElementE].push_back(ssss[0]->set[es][as]);
														subRingAssignment[thisElementE].push_back(sassign);//used to be sassign, but its not clear why. this each ring to start with assignment 0
													}
												}
												delete ssss[0];
												sassign++;
												
											}
											else {
												//std::cout << "pushing element to set in " << elements[thisElement] << std::endl;
												subringOfAtoms[thisElement].push_back(atom(x, y, 0));
												subRingAssignment[thisElement].push_back(sassign++/*sassign++*/);//used to be sassign, but its not clear why. this each ring to start with assignment 0
											}
											//if (debu) std::cout << "D10" << std::endl;
											//needs to be added to a set by element type

										}
										else {

										}
										//std::cout << "RM-- 10" << std::endl;
										//if (debu) std::cout << "D4: atoms survived" << std::endl;
									}
									//std::cout << "RM-- 11" << std::endl;
									ringOfAtoms[0];
									//std::cout << subringOfAtoms[0].size();
									subringOfAtoms[0];
									//std::cout << "RM-- 11.05" << std::endl;
									std::vector<std::vector<atom>> combinedSubRing = combineSets(ringOfAtoms, subringOfAtoms);
									//std::cout << "RM-- 11.1" << std::endl;
									std::vector<std::vector<int>> combinedAssignment = combineSets(ringAssignment, subRingAssignment);
									//std::cout << "RM-- 11.2" << std::endl;
									structure checkS(combinedSubRing, elements, combinedAssignment);
									//std::cout << "RM-- 11.3" << std::endl;
									//std::cout << "min found: " << minFound << std::endl;
									//std::cout << "max found: " << maxFound << std::endl;
									//checkS.print();
									if (!minFound) {
										//std::cout << "Min search" << std::endl;
										minFound = covalentCriteria(&checkS, covalentCriteriaPercent);
										minRadius = radius;
									}
									else {
										//std::cout << "max search" << radius << std::endl;
										//std::cout << bondingRequirement(&checkS, covalentCriteriaPercent,1);
										//maxFound = !(covalentCriteria(&checkS, covalentCriteriaPercent) && bondingRequirement(&checkS,covalentCriteriaPercent,1));
										maxFound = !connectednessRequirement(&checkS, covalentCriteriaPercent, 1);
										//maxFound = !bondingRequirement(&checkS, covalentCriteriaPercent, 1);
										if (!maxFound) maxRadius = radius;
									}
									//std::cout << "miF:" << minFound << "maF:" << maxFound << std::endl;
									radius += achange;
									/*
									subRingApproved = covalentCriteria(&checkS, covalentCriteriaPercent);//connectednessRequirement(&checkS, covalentCriteriaPercent, 1) && 
									//std::cout << "RM-- 11.4" << std::endl;

									if (subRingApproved)
									{
										//std::cout << "sub ring approved, radius: " << radius << std::endl;
										ringAssignment = combinedAssignment;
										ringOfAtoms = combinedSubRing;
									}
									else {
										//std::cout << "sub ring not approved, radius: " << radius << std::endl;
										radius += achange;
									}
									//std::cout << "RM-- 12" << std::endl;
									lastRadius = radius;
									*/
									if (radiusSteps > maxRadiusSteps) {
										validCombination = false;
										std::cout << "failed due to radius" << std::endl;
									}
									if (snatoms == 1)
									{
										std::cout << "snatoms = 1" << std::endl;
										minFound = true;
										maxFound = true;
										
									}
								}
								if (validCombination)
								{
									//std::cout << "max rad:" << maxRadius << " minRadius " << minRadius;
									double chosenRadius = minRadius / 2 + maxRadius / 2;
									//std::cout << "Chosen radius: " << chosenRadius << std::endl;
									radius = chosenRadius;
									//std::cout << "using radius: " << radius << std::endl;
									tatomsPrior = atomsPrior;//counts the atoms we already visited, for access from element distribution
									tatomsPriorHoles = atomsPriorHoles;//counts atoms we alreadt visited - ones, to access hole distribution
									tatomIterator = atomIterator;
									tatomIteratorHoles = atomIteratorHoles;


									//std::cout << "RM-- 7 looping for subring approval" << std::endl;
									std::vector<std::vector<atom>> subringOfAtoms;
									std::vector<std::vector<int>> subRingAssignment;
									int sassign = assign;//this assignment actually gets used. the one in the loop appears to get erased. we must start off 
									//regardless this is not in a while loop and will only occur once per ring so we can use assign
									for (int e = 0; e < elements.size(); e++) {
										subringOfAtoms.push_back({});
										subRingAssignment.push_back({});
									}
									//std::cout << "RM-- 8" << std::endl;
									for (int a = 0; a < (*ringCombination)[ring][subring]; a++)
									{
										//std::cout << "RM-- 9:" << a << std::endl;
										bool thisHole = true;
										if ((*ringCombination)[ring][subring] != 1) {
											if (atomPrintout) std::cout << "H" << tatomIteratorHoles << ',' << (*holeLocation)[tatomIteratorHoles] << ':';
											thisHole = (*holeLocation)[tatomIteratorHoles];//somehow 1 is for no hole now
											tatomIteratorHoles++;
										}

										if (thisHole)
										{
											int thisElement = (*elementDistribution)[tatomIterator];

											tatomIterator++;

											double angle = a * 2 * M_PI / (*ringCombination)[ring][subring];

											double x = radius * cos(angle);
											double y = radius * sin(angle);
											if (atomPrintout) std::cout << elements[thisElement] << " ";

											if ((*ringCombination)[ring][subring] == 1)
											{
												x = 0;
												y = 0;
											}

											if (elements[thisElement].size() > 2)
											{
												//this code block and ones like it are for when our element is a group of atoms and not a single atom. this can be extremely confusing to understand without that knowledge.
												std::vector<structure* > ssss = structuresFromXYZ(elements[thisElement]);
												moveSet(ssss[0]->set, x, y, 0);
												for (int es = 0; es < ssss[0]->set.size(); es++)
												{
													std::string celement = ssss[0]->elements[es];
													int thisElementE = 0;
													for (int ee = 0; ee < elements.size(); ee++) if (celement == elements[ee]) {
														thisElementE = ee;
													}
													for (int as = 0; as < ssss[0]->set[es].size(); as++)
													{
														ssss[0]->set[es][as].fixedGroup = thisElement;
														subringOfAtoms[thisElementE].push_back(ssss[0]->set[es][as]);
														subRingAssignment[thisElementE].push_back(assign);//this used to be sassign,  but it reset everything to 0. this is not while looped so its not clear why we cant use assign
													}
												}
												delete ssss[0];
												assign++; 
											}
											else {
												//std::cout << "pushing element to set in " << elements[thisElement] << std::endl;
												subringOfAtoms[thisElement].push_back(atom(x, y, 0));
												subRingAssignment[thisElement].push_back(assign++);
											}
											//if (debu) std::cout << "D10" << std::endl;
											//needs to be added to a set by element type

										}
										else {

										}
										//std::cout << "RM-- 10" << std::endl;
										//if (debu) std::cout << "D4: atoms survived" << std::endl;
									}
									//std::cout << "RM-- 11" << std::endl;
									ringOfAtoms[0];
									//std::cout << subringOfAtoms[0].size();
									subringOfAtoms[0];
									//std::cout << "RM-- 11.05" << std::endl;
									std::vector<std::vector<atom>> combinedSubRing = combineSets(ringOfAtoms, subringOfAtoms);
									//std::cout << "RM-- 11.1" << std::endl;
									std::vector<std::vector<int>> combinedAssignment = combineSets(ringAssignment, subRingAssignment);
									//std::cout << "RM-- 11.2" << std::endl;
									structure checkS(combinedSubRing, elements, combinedAssignment);

									ringAssignment = combinedAssignment;
									ringOfAtoms = combinedSubRing;
									lastRadius = radius;
								}

								atomsPrior = tatomsPrior;//counts the atoms we already visited, for access from element distribution
								atomsPriorHoles = tatomsPriorHoles;//counts atoms we alreadt visited - ones, to access hole distribution
								atomIterator = tatomIterator;
								atomIteratorHoles = tatomIteratorHoles;




								//if (atomPrintout) std::cout << ") ";
								//if (debu) std::cout << "D3: subring SURVIVED" << std::endl;
								//std::cout << std::endl << "LRAR now: " << lastRingAvgRadii << std::endl;
								/*if (lastRingAvgRadii > 1000)
								{
									std::cout << "lastRingAvgRadii fail: " << lastRingAvgRadii << ", totalradii: " << totalRadii << ", ratoms: " << ratoms << std::endl;
								}*/

								//std::cout << "LRAR after: " << lastRingAvgRadii << ", totalRadii: " << totalRadii << " ratoms: " << ratoms << std::endl;

							}
							bool ringApproved = false;
							double z = lastZ + achange;

							double minZ = 0;
							double maxZ = 0;
							bool minZfound = false;
							bool maxZfound = false;
							//std::cout << "RM-- 13" << std::endl;
							moveSet(ringOfAtoms, 0, 0, z);//move to the last z to start
							int zsteps = 0;
							int zstepMAX = double(5) / achange;
							while (!minZfound && !maxZfound && validCombination)
							{
								//std::cout << "stuck ? Z" << std::endl;
								//std::cout << "RM-- 14: ring approval" << std::endl;
								moveSet(ringOfAtoms, 0, 0, achange);
								z += achange;
								std::vector<std::vector<atom>> combinedRing = combineSets(growingSet, ringOfAtoms);
								std::vector<std::vector<int>> combinedAssignment = combineSets(growingAssignment, ringAssignment);
								structure checkS(combinedRing, elements, combinedAssignment);
								if (!minZfound) {
									minZfound = covalentCriteria(&checkS, covalentCriteriaPercent);//&& connectednessRequirement(&checkS, covalentCriteriaPercent, 1) &&
									minZ = z;
									//std::cout << "minzfoudn" << std::endl;
								}
								else {
									maxZfound = !connectednessRequirement(&checkS, covalentCriteriaPercent, 1);
									
									//maxZfound = covalentCriteria(&checkS, covalentCriteriaPercent);//&& connectednessRequirement(&checkS, covalentCriteriaPercent, 1) &&
									if(!maxZfound) maxZ = z;
									//std::cout << "maxzfoudn" << std::endl;

								}
								/*
								if (ringApproved)
								{
									//std::cout << "ring approved" << std::endl;
									growingAssignment = combinedAssignment;
									growingSet = combinedRing;
									lastZ = z;
									//std::cout << "ring approved Z: " << z << std::endl;
								}//else std::cout << "ring not approved, Z:" << z << std::endl;
								//if (debu) std::cout << "D2: ring SURVIVED" << std::endl;
								*/
								zsteps++;
								if (zsteps > zstepMAX) {
									validCombination = false;
									//std::cout << "failed due to z" << std::endl;
								}
							}
							if (validCombination)
							{
								z = minZ / 2 + maxZ / 2;
								moveSet(ringOfAtoms, 0, 0, z - maxZ -achange);//currently at maxZ + achange, move back to z this way
								std::vector<std::vector<atom>> combinedRing = combineSets(growingSet, ringOfAtoms);
								std::vector<std::vector<int>> combinedAssignment = combineSets(growingAssignment, ringAssignment);
								structure checkS(combinedRing, elements, combinedAssignment);
								growingAssignment = combinedAssignment;
								growingSet = combinedRing;
								lastZ = z;
							}
							//std::cout << "RM-- 15" << std::endl;


						}
						if (validCombination) {
							set = growingSet;
							assignment = growingAssignment;
							std::cout << "accepted valid combination" << std::endl;
						}
						else {
							std::cout << "threw out invalid comination" << std::endl;
						}
					}
					else {
						/*
						Finally we can create the structure!
						*/
						/*
						For each ring we must determine a center radii (for each atom type) that allows atoms in the ring to be spaced out well.
						If there are center atoms, then we will choose whichever radii is larger, the one from the distance from the center atom or between atoms.
						As a ring may have many types of atoms, each will be farther from the center and this will be determined for each element type individually as if they were the only one in the ring.
						*/
						int atomsPrior = 0;//counts the atoms we already visited, for access from element distribution
						int atomsPriorHoles = 0;//counts atoms we alreadt visited - ones, to access hole distribution
						int atomIterator = 0;
						int atomIteratorHoles = 0;
						double lastZ = 0;
						for (int ring = 0; ring < ringCombination->size(); ring++)
						{
							//if (debu) std::cout << "D2 rings" << std::endl;
							//for each subring
							std::vector < std::vector<atom>> ringOfAtoms;
							if (atomPrintout) std::cout << "{";
							//determine z distance from last ring
							double z = 0;
							if (ring != 0)
							{
								//find the two closest n values
								//use the nvalues from the centermost ring of this ring and the previous, again this is just an estimate.
								double radii1 = averageCovalentRadiiVertical;
								//std::cout << "radii1: " << radii1 << std::endl;
								double lastRadii2 = averageCovalentRadiiVertical / sin(M_PI / ((*ringCombination)[ring - 1][0]));//the radius where the distance between atoms is the covalent diameter, the angle between atoms must be divided by 2 and thereofre does not have 2 M_PI
								//std::cout << "lastradii2: " << lastRadii2 << std::endl;
								double Radii2 = averageCovalentRadiiVertical / sin(M_PI / ((*ringCombination)[ring][0]));//the radius where the distance between atoms is the covalent diameter, the angle between atoms must be divided by 2 and thereofre does not have 2 M_PI
								//std::cout << "Radii2: " << Radii2 << std::endl;
								//std::cout << "Angle for radii2: " << M_PI / ((*ringCombination)[ring][0]) << ", number: " << ((*ringCombination)[ring][0]) << " for ring: " << ring << std::endl;
								double lastR = 0;
								double thisR = 0;
								if (Radii2 > radii1) thisR = Radii2;
								else thisR = radii1;
								if (lastRadii2 > radii1) lastR = lastRadii2;
								else lastR = radii1;
								if ((*ringCombination)[ring][0]) thisR = 0;
								if ((*ringCombination)[ring - 1][0]) lastR = 0;
								//std::cout << "asin: " << (thisR - lastR) / averageCovalentRadii << " , " << thisR << " " << lastR << " " << averageCovalentRadii << std::endl;
								double theta = asin((thisR - lastR) / averageCovalentRadiiVertical);
								if (std::isnan(theta))
								{
									std::cout << "isnan theta:, " << (thisR - lastR) / averageCovalentRadiiVertical << ",thisR: " << thisR << ",lastR " << lastR << std::endl;
								}
								//std::cout << "theta: " << theta << std::endl;
								z = 2 * averageCovalentRadiiVertical * cos(theta) + lastZ;
								//std::cout << "z: " << z << std::endl;
								if (z < 0) z = -z;
								if (std::isnan(z)) {
									std::cout << "nan detected:" << std::endl;
									std::cout << "acr" << averageCovalentRadiiVertical << ",cos: " << cos(theta) << ",theta: " << theta << ", " << lastZ << std::endl;
								}
							}
							double lastRingAvgRadii = 0;

							for (int subring = (*ringCombination)[ring].size() - 1; subring >= 0; subring--)
							{
								// std::cout << " (S:" << (*ringCombination)[ring][subring] << ";";
								//if (debu) std::cout << "D3: subring " << std::endl;
								//calculate the radii for each element type
								std::vector<double> elementalRadii = {};
								for (std::vector<std::string>::iterator e = elements.begin(); e != elements.end(); e++)
								{
									double radii;
									//std::cout << "last ring avg radii: " << lastRingAvgRadii << std::endl << "_______________" << std::endl;
									/*if (subring == 0) */ radii = covalentRadii(*e, 1, true) + 2 * lastRingAvgRadii * 100;//100 for pm conversion, 2 to make it a diameter
									//else radii = 2 * covalentRadii(*e, 1);//the atoms in the ring below may not be the same type, but this is just an estimate anyways
									double radii2 = covalentRadii(*e, 1, true) / sin(M_PI / ((*ringCombination)[ring][subring]));//the radius where the distance between atoms is the covalent diameter, the angle between atoms must be divided by 2 and thereofre does not have 2 M_PI
									if ((*ringCombination)[ring][subring] == 1) radii2 = 0;
									/*if (radii2 > 1000)
									{
										std::cout << "problem starter: " << radii2 << ", from ring: " << (*ringCombination)[ring][subring];
									}*/
									if (radii > radii2) elementalRadii.push_back(radii);
									else elementalRadii.push_back(radii2);
									//std::cout << "elemental radii: "<< elementalRadii[elementalRadii.size() - 1] << std::endl;
								}
								int ratoms = 0;
								double totalRadii = 0;
								for (int a = 0; a < (*ringCombination)[ring][subring]; a++)
								{

									//if (debu) std::cout << "D4: atoms" << std::endl;

									//if (debu) std::cout << "element distribution size/natoms:" << natoms << " atom iterator: " << atomIterator << std::endl;
									bool thisHole = true;
									if ((*ringCombination)[ring][subring] != 1) {
										if (atomPrintout) std::cout << "H" << atomIteratorHoles << ',' << (*holeLocation)[atomIteratorHoles] << ':';
										thisHole = (*holeLocation)[atomIteratorHoles];//somehow 1 is for no hole now
										atomIteratorHoles++;
									}
									else if (atomPrintout) std::cout << "Hskip:";

									//if (debu) std::cout << "D5" << std::endl;

									if (thisHole)
									{

										ratoms++;
										int thisElement = (*elementDistribution)[atomIterator];
										totalRadii += elementalRadii[thisElement] / 100.0;

										atomIterator++;
										//if (debu) std::cout << "D6" << std::endl;
										double angle = a * 2 * M_PI / (*ringCombination)[ring][subring];
										//if (debu) std::cout << "D7:" << thisElement << " in " << elementalRadii.size() << std::endl;

										/*
										NEW Change 8/17
										The distance between two elements next to eachother should be the priority rather than the distance between the oposite side of the polygon
										polygon radius = 2(elementalRadius)/2sin(angle/2)
										
										*/
										double polyR = elementalRadii[thisElement] / sin(angle / 2);
										//double x = elementalRadii[thisElement] / 100.0 * cos(angle);
										double x = polyR / 100.0 * cos(angle);

										//if (debu) std::cout << "D8" << std::endl;
										
										//double y = elementalRadii[thisElement] / 100.0 * sin(angle);
										
										double y = polyR / 100.0 * sin(angle);
										std::cout << "chord length: " << elementalRadii[thisElement] << std::endl;

										//std::cout << "ANGLE,x,y: " << angle << " " << x << " " << y <<std::endl;
										//if (debu) std::cout << "D9" << std::endl;
										if (atomPrintout) std::cout << elements[thisElement] << " ";
										//std::cout << "added an atom at ring: " << ring << ", subring: " << subring << std::endl;

										if ((*ringCombination)[ring][subring] == 1)
										{
											x = 0;
											y = 0;
										}
										/*if ((x > 1000 && x > 0) || (x < -1000 && x < 0) || (y > 1000 && y > 0) || (y < -1000 && y < 0))
										{
											std::cout << "strange numbers found." << std::endl;
											std::cout << "z: " << z << std::endl;
											std::cout << "elemental radii: " << elementalRadii[thisElement] << std::endl;
											std::cout << "angle: " << angle << std::endl;
											std::cout << "a: " << a << "ring#:  " << (*ringCombination)[ring][subring] << std::endl;
											std::cout << "cos: " << cos(angle) << ", sin: " << sin(angle) << std::endl;
											std::cout << "x: " << x << ",y: " << y << std::endl << std::endl;
											std::cout << "lastringavgradii: " << lastRingAvgRadii << std::endl;
										}*/

										if (elements[thisElement].size() > 2)
										{

											std::vector<structure* > ssss = structuresFromXYZ(elements[thisElement]);
											moveSet(ssss[0]->set, x, y, z);
											//std::cout << "checking the sub structure elements: " << std::endl;
											//for (int ec = 0; ec < ssss[0]->elements.size(); ec++) std::cout << ssss[0]->elements[ec] << "" << ssss[0]->set[ec].size() <<" ";
											//std::cout << std::endl;
											//std::cout << " here 5.1?" << ssss.size() << " " << elements[thisElement] << std::endl;
											for (int es = 0; es < ssss[0]->set.size(); es++)
											{

												//figure out the actual this element
												std::string celement = ssss[0]->elements[es];
												int thisElementE = 0;
												//std::cout << "starting element: " << celement << " for set index: " << es << std::endl;
												for (int ee = 0; ee < elements.size(); ee++) if (celement == elements[ee]) {
													thisElementE = ee;
													//std::cout << "element identified in set: " << celement << " " << elements[ee] << " for index " << ee << std::endl;
												}
												//else {
												//	std::cout << "celement:" << celement << ", did not match:" << elements[ee] << ", indeces:" << es << "," << ee << std::endl;
												//}
												for (int as = 0; as < ssss[0]->set[es].size(); as++)
												{

													ssss[0]->set[es][as].fixedGroup = thisElement;
													//std::cout << "pushing element to set in " << elements[thisElement] << " with index " << thisElementE << std::endl;
													set[thisElementE].push_back(ssss[0]->set[es][as]);
													assignment[thisElementE].push_back(assign);

												}

											}
											delete ssss[0];
											assign++;
										}
										else {
											//std::cout << "pushing element to set in " << elements[thisElement] << std::endl;
											set[thisElement].push_back(atom(x, y, z));
											assignment[thisElement].push_back(assign++);
										}
										//if (debu) std::cout << "D10" << std::endl;
										//needs to be added to a set by element type

									}

									//if (debu) std::cout << "D4: atoms survived" << std::endl;
								}
								//if (atomPrintout) std::cout << ") ";
								//if (debu) std::cout << "D3: subring SURVIVED" << std::endl;
								//std::cout << std::endl << "LRAR now: " << lastRingAvgRadii << std::endl;

								lastRingAvgRadii += totalRadii / (double)ratoms;
								/*if (lastRingAvgRadii > 1000)
								{
									std::cout << "lastRingAvgRadii fail: " << lastRingAvgRadii << ", totalradii: " << totalRadii << ", ratoms: " << ratoms << std::endl;
								}*/

								//std::cout << "LRAR after: " << lastRingAvgRadii << ", totalRadii: " << totalRadii << " ratoms: " << ratoms << std::endl;

							}
							if (atomPrintout) std::cout << "} ";
							//if (debu) std::cout << "D2: ring SURVIVED" << std::endl;
							lastZ = z;
						}
					}
					//if (debu) std::cout << "D1: inner loop SURVIVED" << std::endl;

					//debugging
					//count atomsin set and see whats wrong
					/*int satoms = 0;
					for (int s = 0; s < set.size(); s++) satoms += set[s].size();
					if (satoms != natoms)
					{
						std::cout << "wrong atoms: " << satoms << " vs. " << natoms << std::endl;
						//print details to see if the mistake has already been made
						std::cout << "printout\n---------------------\nhole locations: ";
						for (int i = 0; i < ringSpots - nOnes; i++) std::cout << (*holeLocation)[i] << " ";
						std::cout << "std::endl";
						std::cout << "EDS: " << elementDistributions.size() << std::endl;
						std::cout << "element distribution:";
						for (int i = 0; i < natoms; i++) std::cout << (*elementDistribution)[i] << " ";
						std::cout << "std::endl";
						std::cout << "ring combinaiton: ";
						for (std::vector<std::vector<int>>::iterator ccit = ringCombination->begin(); ccit != ringCombination->end(); ccit++)
						{
							std::cout << "{";
							for (std::vector<int>::iterator cccit = ccit->begin(); cccit != ccit->end(); cccit++) std::cout << *cccit << " ";
							std::cout << "} ";
						}
						std::cout << std::endl;
						std::cout << "--------------------\n";
					}
					else std::cout << "right atoms: " <<"-------" << natoms << std::endl;*/
					/*std::cout << "adding set with atomic distribution by element: ";
					for (std::vector<std::vector<atom>>::iterator at = set.begin(); at != set.end(); at++) std::cout << at->size() << " ";
					std::cout << std::endl;*/
					
					int setsize = 0;
					for (std::vector<std::vector<atom>>::iterator it = set.begin(); it != set.end(); it++) setsize += it->size();
					if(setsize != 0)  notAllSets.push_back(std::pair<std::vector<std::vector<atom>>, std::vector<std::vector<int>>>(set, assignment));
					//std::cout << "ADDED MOLECULE" << std::endl;
					if (atomPrintout) std::cout << std::endl;
				}
				comboIterator++;
			}
		}
		allSets.push_back(notAllSets);

		/*
		clean holeLocations and multinomial coefficient
		*/
		for (std::vector<bool*>::iterator hit = holeLocations.begin(); hit != holeLocations.end(); hit++) delete[](*hit);
		for (std::vector<int*>::iterator eit = elementDistributions.begin(); eit != elementDistributions.end(); eit++) delete[](*eit);

	
	}

	//make structures from all sets
	std::vector<std::vector<structure*>> allPreCriteriaStructures = {};
	//std::cout << "all sets size: " << allSets.size() << std::endl;
	for (int is = 0; is < allSets.size(); is++) {
		std::vector<structure*> preCriteriaStructures = {};
		for (std::vector<std::pair<std::vector<std::vector<atom>>,std::vector<std::vector<int>>>>::iterator sit = allSets[is].begin(); sit != allSets[is].end(); sit++)
		{
			//std::cout << "elements: ";
			//for (int eee = 0; eee < elements.size(); eee++) std::cout << elements[eee] << " ";
			//std::cout << std::endl;
			//std::cout << "set size: " << sit->first.size() << " ";
			//for (int sss= 0; sss < sit->first.size(); sss++) std::cout << sit->first[sss].size() << " ";
			//std::cout << std::endl;
			structure* struc = new structure(sit->first, elements, sit->second);
			/*std::cout << "adding STRUCTURE with atomic distribution by element: ";
			for (std::vector<std::vector<atom>>::iterator at = struc->set.begin(); at != struc->set.end(); at++) std::cout << at->size() << " ";
			std::cout << std::endl;*/
			//struc->assignment = sit->second;
			
			//std::cout << std::endl;
			//std::cout << "NEXT" << std::endl;
			//struc->print();
			preCriteriaStructures.push_back(struc);

		}
		allPreCriteriaStructures.push_back(preCriteriaStructures);
	}

	//check the bonding criteria only, just to remove ridiculous stuff
	std::vector<structure*> criteriaStructures = {};
	//std::cout << "pre criteria structures: " << preCriteriaStructures.size() << std::endl;
	//std::cout << "all pre criteria structures size: " << allPreCriteriaStructures.size() << std::endl;
	for (int is = 0; is < allPreCriteriaStructures.size(); is++) {
		std::vector<structure*> theseCriteriaStructures = {};
		//std::cout << "pre criteria structures size: " << allPreCriteriaStructures[is].size() << std::endl;
		for (std::vector<structure*>::iterator sit = allPreCriteriaStructures[is].begin(); sit != allPreCriteriaStructures[is].end(); sit++)
		{
			if (true || connectednessRequirement((*sit), covalentCriteriaPercent, 1))//TODO eventually i want to implement bond requirements for each element type
			{
				/*std::cout << "PASSING structure with atomic distribution by element: ";
				for (std::vector<std::vector<atom>>::iterator at = (*sit)->set.begin(); at != (*sit)->set.end(); at++) std::cout << at->size() << " ";
				std::cout << std::endl;*/
				theseCriteriaStructures.push_back(*sit);

			}
			else {
				delete* sit;
				std::cout << "failed criteria" << std::endl;
			}
		}

		if (keepPerRingCombination == -1 || keepPerRingCombination >= theseCriteriaStructures.size())
		{
			//std::cout << "keep per ring combination: "<< keepPerRingCombination << " vs. These criteria structures size : " << theseCriteriaStructures.size() << std::endl;
			for (std::vector<structure*>::iterator csit = theseCriteriaStructures.begin(); csit != theseCriteriaStructures.end(); csit++) criteriaStructures.push_back(*csit);
		}
		else {
			//std::cout << "These criteria structures size: " << theseCriteriaStructures.size() << std::endl;
			std::vector<bool> skeep = randomNchooseC(theseCriteriaStructures.size(), keepPerRingCombination);
			for (int sit = 0; sit < theseCriteriaStructures.size(); sit++)
			{
				if (skeep[sit]) criteriaStructures.push_back(theseCriteriaStructures[sit]);
				else delete theseCriteriaStructures[sit];
			}
		}
		int added = 0;
	}
	if (dOut) std::cout << "total structures made:" << criteriaStructures.size() << std::endl;
	return criteriaStructures;
}
std::vector<structure*> proceduralSeeds(std::vector<int> composition, std::vector<std::string> elements, std::vector<std::tuple<std::vector<int>, std::vector<int>,std::vector<int>>> nRnVaH, double covalentCriteriaPercent,int keepPerCombo, int keepPerRingCombination,double achange){
	std::vector<structure*> retVal;
	for (int i = 0; i < nRnVaH.size(); i++)
	{

		std::cout << "mrnvah combination " << i << std::endl;
		//std::vector<structure*> tempVal = proceduralSeeds(composition, elements, std::get<0>(nRnVaH[i]), std::get<1>(nRnVaH[i]), covalentCriteriaPercent, std::get<2>(nRnVaH[i]), keepPerRingCombination);
		std::vector<int> numRings = std::get<0>(nRnVaH[i]);
		std::vector<int> nValues = std::get<1>(nRnVaH[i]);
		std::vector<int> allowedHoles = std::get<2>(nRnVaH[i]);
		std::vector<structure*> tempTempVal;
		//std::cout << "starting loop " << std::endl;
		for (int j = 0; j < numRings.size(); j++)
		{
			//std::cout << "running seed generation #" << j + j*i + 1 << " of " << numRings.size()*nRnVaH.size() << std::endl;
			//std::cout << "j " << j << std::endl;
			std::vector<structure*> tempVal = proceduralSeeds(composition, elements, numRings[j], nValues, covalentCriteriaPercent, allowedHoles, keepPerCombo,keepPerRingCombination,achange);
			//std::cout << "temp val size: " << tempVal.size() << std::endl;
			/*
			if (keepPerCombo == -1 || keepPerCombo >= tempVal.size()) {
				for (std::vector<structure*>::iterator it = tempVal.begin(); it != tempVal.end(); it++) retVal.push_back(*it);
			}
			else {
				std::vector<bool> skeep = randomNchooseC(tempVal.size(), keepPerCombo);

				/*
				std::cout << "skeep parameters: " << tempVal.size() << " , " << keepPerCombo << std::endl;
				for (int sk = 0; sk < skeep.size(); sk++) std::cout << skeep[sk] << " ";
				std::cout << std::endl;
				

				for (int sit = 0; sit < tempVal.size(); sit++)
				{
					//std::cout << "sit " << sit << std::endl;
					if (skeep[sit]) retVal.push_back(tempVal[sit]);
					else delete tempVal[sit];
				}
			}*/
			for (int sit = 0; sit < tempVal.size(); sit++)
			{
				retVal.push_back(tempVal[sit]);
			}
			
		}
		//std::cout << "finished loop " << std::endl;
	}
	return retVal;
}


std::vector<structure*> combineSurface(structure* surface, std::vector<structure*> seeds, double covalentCriteriaPercent, double achange) {
	std::vector<structure*> newSeeds = {};
	double sx = 0.0;
	double sy = 0.0;
	double sz = 0.0;
	double szmax = -1;
	double szmin = -1;
	for (std::vector<structure*>::iterator it = seeds.begin(); it != seeds.end(); it++)
	{
		bool works = false;
		while (!works) {
			sz += achange;
			std::vector<std::vector<atom>> set;
			std::vector<std::string> elements;
			int elIt = 0;
			std::vector<std::vector<int>> assignment;
			int assignIt = 0;
			for (std::vector<std::vector<atom>>::iterator sit = (*it)->set.begin(); sit != (*it)->set.end(); sit++)
			{
				std::vector<atom> celement = {};
				std::vector<int> assi = {};
				for (std::vector<atom>::iterator ait = sit->begin(); ait != sit->end(); ait++)
				{
					celement.push_back(atom(ait->x, ait->y, ait->z));
					assi.push_back((*it)->assignment[assignIt++]);
				}
				elements.push_back(elements[elIt++]);
				set.push_back(celement);
				assignment.push_back(assi);

			}

			elIt = 0;
			for (std::vector<std::vector<atom>>::iterator sit = surface->set.begin(); sit != surface->set.end(); sit++)
			{
				//check the element
				int eIndex = -1;
				for (int e = 0; e < (*it)->elements.size(); e++)
				{
					if ((*it)->elements[e] == surface->elements[elIt])
					{
						eIndex = e;
					}
				}
				if (eIndex == -1)
				{
					elements.push_back(surface->elements[elIt]);
					eIndex = elements.size() - 1;
					set.push_back({});
					assignment.push_back({});
				}

				for (std::vector<atom>::iterator ait = sit->begin(); ait != sit->end(); ait++)
				{
					set[eIndex].push_back(atom(ait->x + sx, ait->y + sy, ait->z + sz));
					assignment[eIndex].push_back(-1);
				}
				elIt++;
			}


			structure* cstruct = new structure(set, elements, assignment, false);
			works = covalentCriteria(cstruct, covalentCriteriaPercent) && bondingRequirement(cstruct, covalentCriteriaPercent, 1);

		}

	}
	return {};

	
}
/*
std::vector<int>* copyAssignment(std::vector<int>* oldAs)
{
	std::vector<int>* newAssignment = new std::vector<int>;
	for (int i = 0; i < oldAs->size(); i++)
	{
		newAssignment->push_back((*oldAs)[i]);
	}
	return newAssignment;
}
*/