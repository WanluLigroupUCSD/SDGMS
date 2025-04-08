#include "calculators.hpp"
const int DEBUG = 0;
int getBasisSetElectrons(std::string element, int complexity)
{
	if (element == "H") {
		return 1;
	}
	else if (element == "He") {
		return 2;
	}
	else if (element == "Li") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 3;
	}
	else if (element == "Be") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 4;
	}
	else if (element == "B") {
		return 3;
	}
	else if (element == "C") {
		return 4;
	}
	else if (element == "N") {
		return 5;
	}
	else if (element == "O") {
		return 6;
	}
	else if (element == "F") {
		return 7;
	}
	else if (element == "Ne") {
		return 8;
	}
	else if (element == "Na") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 9;
	}
	else if (element == "Mg") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 10;
	}
	else if (element == "Al") {
		return 3;
	}
	else if (element == "Si") {
		return 4;
	}
	else if (element == "P") {
		return 5;
	}
	else if (element == "S") {
		return 6;
	}
	else if (element == "Cl") {
		return 7;
	}
	else if (element == "Ar") {
		return 8;
	}
	else if (element == "K") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 9;
	}
	else if (element == "Ca") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 10;
	}
	else if (element == "Sc") {
		if (complexity == 1) return 3;
		else if (complexity == 2) return 11;
	}
	else if (element == "Ti") {
		if (complexity == 1) return 4;
		else if (complexity == 2) return 12;
	}
	else if (element == "V") {
		if (complexity == 1) return 5;
		else if (complexity == 2) return 13;
	}
	else if (element == "Cr") {
		if (complexity == 1) return 6;
		else if (complexity == 2) return 14;
	}
	else if (element == "Mn") {
		if (complexity == 1) return 7;
		else if (complexity == 2) return 15;
	}
	else if (element == "Fe") {
		if (complexity == 1) return 8;
		else if (complexity == 2) return 16;
	}
	else if (element == "Co") {
		if (complexity == 1) return 9;
		else if (complexity == 2) return 17;
	}
	else if (element == "Ni") {
		if (complexity == 1) return 10;
		else if (complexity == 2) return 18;
	}
	else if (element == "Cu") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 11;
		else if (complexity == 3) return 19;
	}
	else if (element == "Zn") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 12;
		else if (complexity == 3) return 20;
	}
	else if (element == "Ga") {
		if (complexity == 1) return 3;
		else if (complexity == 2) return 13;
	}
	else if (element == "Ge") {
		return 4;
	}
	else if (element == "As") {
		return 5;
	}
	else if (element == "Se") {
		return 6;
	}
	else if (element == "Br") {
		return 7;
	}
	else if (element == "Kr") {
		return 8;
	}
	else if (element == "Rb") {
		return 1;
	}
	else if (element == "Sr") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 10;
	}
	else if (element == "Y") {
		if (complexity == 1) return 3;
		else if (complexity == 2) return 11;
	}
	else if (element == "Zr") {
		if (complexity == 1) return 4;
		else if (complexity == 2) return 12;
	}
	else if (element == "Nb") {
		if (complexity == 1) return 5;
		else if (complexity == 2) return 13;
	}
	else if (element == "Mo") {
		if (complexity == 1) return 6;
		else if (complexity == 2) return 14;
	}
	else if (element == "Tc") {
		if (complexity == 1) return 7;
		else if (complexity == 2) return 15;
	}
	else if (element == "Ru") {
		if (complexity == 1) return 8;
		else if (complexity == 2) return 16;
	}
	else if (element == "Rh") {
		if (complexity == 1) return 9;
		else if (complexity == 2) return 17;
	}
	else if (element == "Pd") {
		if (complexity == 1) return 10;
		else if (complexity == 2) return 18;
	}
	else if (element == "Ag") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 11;
		else if (complexity == 3) return 19;
	}
	else if (element == "Cd") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 12;
	}
	else if (element == "In") {
		if (complexity == 1) return 3;
		else if (complexity == 2) return 13;
		else if (complexity == 3) return 146;
	}
	else if (element == "Sn") {
		return 4;
	}
	else if (element == "Sb") {
		return 5;
	}
	else if (element == "Te") {
		return 6;
	}
	else if (element == "I") {
		if (complexity == 1) return 7;
		else if (complexity == 2) return 129;
		else if (complexity == 3) return 125;
	}
	else if (element == "Xe") {
		return 8;
	}
	else if (element == "Cs") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 9;
	}
	else if (element == "Ba") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 10;
	}
	else if (element == "La") {
		return 11;
	}
	else if (element == "Ce") {
		return 12;
	}
	else if (element == "Pr") {
		return 13;
	}
	else if (element == "Nd") {
		return 14;
	}
	else if (element == "Pm") {
		return 15;
	}
	else if (element == "Sm") {
		return 16;
	}
	else if (element == "Eu") {
		return 17;
	}
	else if (element == "Gd") {
		return 18;
	}
	else if (element == "Tb") {
		return 19;
	}
	else if (element == "Dy") {
		return 20;
	}
	else if (element == "Ho") {
		return 21;
	}
	else if (element == "Er") {
		return 22;
	}
	else if (element == "Tm") {
		return 23;
	}
	else if (element == "Yd") {
		return 24;
	}
	else if (element == "Lu") {
		return 25;
	}
	else if (element == "Hf") {
		return 12;
	}
	else if (element == "Ta") {
		if (complexity == 1) return 5;
		else if (complexity == 2) return 13;
	}
	else if (element == "W") {
		if (complexity == 1) return 6;
		else if (complexity == 2) return 14;
	}
	else if (element == "Re") {
		if (complexity == 1) return 7;
		else if (complexity == 2) return 15;
	}
	else if (element == "Os") {
		if (complexity == 1) return 8;
		else if (complexity == 2) return 16;
	}
	else if (element == "Ir") {
		if (complexity == 1) return 9;
		else if (complexity == 2) return 17;
	}
	else if (element == "Pt") {
		if (complexity == 1) return 10;
		else if (complexity == 2) return 18;
	}
	else if (element == "Au") {
		if (complexity == 1) return 1;
		else if (complexity == 2) return 11;
		else if (complexity == 3) return 19;
	}
	else if (element == "Hg") {
		if (complexity == 1) return 2;
		else if (complexity == 2) return 12;
	}
	else if (element == "Tl") {
		if (complexity == 1) return 3;
		else if (complexity == 2) return 13;
	}
	else if (element == "Pb") {
		return 4;
	}
	else if (element == "Bi") {
		return 5;
	}
	else if (element == "Po") {
		return 6;
	}
	else if (element == "At") {
		return 7;
	}
	else if (element == "Rn") {
		return 8;
	}
	else if (element == "Fr") {
		//no basis?
	}
	else if (element == "Ra") {
		//no basis?
	}
	else if (element == "Ac") {
		if (complexity == 1) return 11;
		else if (complexity == 2) return 29;
	}
	else if (element == "Th") {
		if (complexity == 1) return 12;
		else if (complexity == 2) return 30;
	}
	else if (element == "Pa") {
		if (complexity == 1) return 13;
		else if (complexity == 2) return 31;
	}
	else if (element == "U") {
		if (complexity == 1) return 14;
		else if (complexity == 2) return 32;
	}
	else if (element == "Np") {
		if (complexity == 1) return 15;
		else if (complexity == 2) return 33;
	}
	else if (element == "Pu") {
		if (complexity == 1) return 16;
		else if (complexity == 2) return 34;
	}
	else if (element == "Am") {
		if (complexity == 1) return 17;
		else if (complexity == 2) return 35;
	}
	else if (element == "Cm") {
		if (complexity == 1) return 18;
		else if (complexity == 2) return 36;
	}
	else if (element == "Bk") {
		if (complexity == 1) return 20;
		else if (complexity == 2) return 38;
	}
	else if (element == "Cf") {
		//no basis?
	}
	else if (element == "Es") {
		if (complexity == 1) return 21;
		else if (complexity == 2) return 39;
	}
	else if (element == "Fm") {
		if (complexity == 1) return 22;
		else if (complexity == 2) return 40;
	}
	else if (element == "Md") {
		if (complexity == 1) return 23;
		else if (complexity == 2) return 41;
	}
	else if (element == "No") {
		if (complexity == 1) return 24;
		else if (complexity == 2) return 42;
	}
	else if (element == "Lr") {
		if (complexity == 1) return 25;
		else if (complexity == 2) return 43;
	}
	else if (element == "Rf") {
		//no basis
	}
	else if (element == "Db") {
		//no basis
	}
	else if (element == "Sg") {
		//no basis
	}
	else if (element == "Bh") {
		//no basis
	}
	else if (element == "Hs") {
		//no basis
	}
	else if (element == "Mt") {
		//no basis
	}
	else if (element == "Ds") {
		//no basis
	}
	else if (element == "Rg") {
		//no basis
	}
	else if (element == "112") {
		//no basis
	}
	else if (element == "113") {
		//no basis
	}
	else if (element == "114") {
		//no basis
	}
	else if (element == "115") {
		//no basis
	}
	else if (element == "116") {
		//no basis
	}
	else if (element == "117") {
		//no basis
	}
	else if (element == "118") {
		//no basis
	}
	return 0;
}
std::string getBasisSet(std::string element, int complexity)
{
	return "GTH-PADE-q" + std::to_string(getBasisSetElectrons(element, complexity));
}

void writeToCP2K(structure& s, std::string task, std::string filename, std::string basis_file, std::string potential_file, double a, double b, double c)
{
	std::string output;
	std::ofstream file(filename);

	std::string global = "&GLOBAL\n\tPROJECT " + task + "\n\tRUN_TYPE ENERGY_FORCE\n\tPRINT_LEVEL LOW\n&END GLOBAL";

	std::string DFT = "\n\t&DFT"
		"\n\t\tBASIS_SET_FILE_NAME " + basis_file;
	"\n\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS " + potential_file;
	"\n\t\t&QS"
		"\n\t\t\tEPS_DEFAULT 1.0E-10"
		"\n\t\t&END QS"
		"\n\t\t&MGRID"
		"\n\t\t\tNGRIDS 4"
		"\n\t\t\tCUTOFF 300"
		"\n\t\t\tREL_CUTOFF 60"
		"\n\t\t&END MGRID"
		"\n\t\t&XC"
		"\n\t\t\t&XC_FUNCTIONAL PADE"
		"\n\t\t\t&END XC_FUNCTIONAL"
		"\n\t\t&END XC"
		"\n\t\t&SCF"
		"\n\t\t\tSCF_GUESS ATOMIC"
		"\n\t\t\tEPS_SCF 1.0E-7"
		"\n\t\t\tMAX_SCF 300"
		"\n\t\t\t&DIAGONALIZATION ON"
		"\n\t\t\t\tALGORITHM STANDARD"
		"\n\t\t\t&END DIAGONALIZATION"
		"\n\t\t\t&MIXING T"
		"\n\t\t\t\tMETHOD BROYDEN_MIXING"
		"\n\t\t\t\tALPHA 0.4"
		"\n\t\t\t\tNBROYDEN 8"
		"\n\t\t\t&END MIXING"
		"\n\t\t&END SCF"
		"\n\t&END DFT";

	std::string subsys = "\n\t&SUBSYS";
	for (int e = 0; e < s.elements.size(); e++)
	{
		subsys += "\n\t\t&KIND " + s.elements[e];
		subsys += "\n\t\t\tELEMENT\t" + s.elements[e];
		subsys += "\n\t\t\tBASIS_SET DZVP-GTH-PADE"
			"\n\t\t\tPOTENTIAL " + getBasisSet(s.elements[e], 1) +
			"\n\t\t&END KIND";
	}

	/*
	We need to know the cell shape. From a materials science perspective it makes sense for this to be the smallest shape that will fit all the atoms in side. Determine this.
	*/
	if (a == 0 && b == 0 && c == 0)
	{
		double maxa = 0;
		double maxb = 0;
		double maxc = 0;
		double mina = DBL_MAX;
		double minb = DBL_MAX;
		double minc = DBL_MAX;

		for (int e = 0; e < s.set.size(); e++)
		{
			for (int atom = 0; atom < s.set[e].size(); atom++)
			{
				if (s.set[e][atom].x > maxa) {
					maxa = s.set[e][atom].x;
				}
				else if (s.set[e][atom].x < mina)
				{
					mina = s.set[e][atom].x;
				}

				if (s.set[e][atom].y > maxb) {
					maxb = s.set[e][atom].y;
				}
				else if (s.set[e][atom].y < minb)
				{
					minb = s.set[e][atom].y;
				}

				if (s.set[e][atom].z > maxc) {
					maxc = s.set[e][atom].z;
				}
				else if (s.set[e][atom].z < minc)
				{
					minc = s.set[e][atom].z;
				}

			}
		}
		//double the cell shape to account for the negative numbers?
		a = maxa - mina;
		b = maxb - minb;
		c = maxc - minc;
	}

	subsys += "\n\t\t&CELL"
		"\n\t\t\tABC " + std::to_string(a) + " " + std::to_string(b) + " " + std::to_string(c) +
		"\n\t\t&END CELL"
		//"\n\t\t# " + std::to_string(numMolecules) + " " + formula + "(TIP5P," + std::to_string(pressureBars) + "bar, " + std::to_string(temperature) + "K) a = 9.8528"
		"\n\t\t&COORD";
	for (int e = 0; e < s.set.size(); e++)
	{
		for (int atom = 0; atom < s.set[e].size(); atom++)
		{
			subsys += "\n\t" + s.elements[e] + " " + std::to_string(s.set[e][atom].x) + " " + std::to_string(s.set[e][atom].y) + " " + std::to_string(s.set[e][atom].z);
		}
	}
	subsys += "\n\t\t&END COORD"
		"\n\t&END SUBSYS";

	std::string print = "\n\t&PRINT"
		"\n\t\t&FORCES ON"
		"\n\t\t&END FORCES"
		"\n\t&END PRINT";

	std::string forceeval = "&FORCE_EVAL\n\tMETHOD Quickstep" + DFT + subsys + print + "\n&END FORCE_EVAL";
	file << global << "\n" << forceeval;
	file.close();
}

void writeToGuassian(structure& s, std::string task, std::string filename, std::string basis_name, std::string method, int charge, int state, int maxCycles, int maxSCF, std::string customBasis)
{
	//max cycles = 0 implies this is just energy and no optimization
	//a,b,c are useless for guassian? so is potential name. the variables are left here for compatibility.
	/*

	# HF/6-31G(d)	Route section

water energy	Title section

0   1	Molecule specification
O  -0.464   0.177   0.0
H  -0.464   1.137   0.0
H   0.441  -0.143   0.0
	*/
	//route section
	std::string outString = "# " + method + "/" + basis_name;
	if (maxCycles > 0) outString += " Opt=(MaxCycles=" + std::to_string(maxCycles) + ")";
	if (maxSCF > 0) outString += " SCF = (MaxCycle = " + std::to_string(maxSCF) + ")";
	outString += "\n\n";
	//title section
	outString += task + "\n\n";
	//molecule specification, charge, state singlet doublet triplet
	outString += std::to_string(charge) + " " + std::to_string(state);
	//molecules
	for (int e = 0; e < s.set.size(); e++)
	{
		for (int atom = 0; atom < s.set[e].size(); atom++)
		{
			outString += "\n" + s.elements[e] + " " + std::to_string(s.set[e][atom].x) + " " + std::to_string(s.set[e][atom].y) + " " + std::to_string(s.set[e][atom].z);
		}
	}
	outString += "\n\n";//blank line at the end of the structure is required
	std::ofstream file(filename);
	file << outString;
	if (customBasis != "empty")
	{
		//add a custom basis file to the end
		std::ifstream cB(customBasis);
		std::string cbbline;
		while (std::getline(cB, cbbline))
		{
			file << cbbline;
		}
	}
	file.close();



}
void writeToADF(structure& s, std::string task, std::string filename, std::string basis_name, std::string method, int maxSCF, int charge, int state, std::string converge, bool unrestricted)
{

	//task is either GeometryOptimization or Energy??
	std::string out = "#!/bin/sh\n\n\"$AMSBIN/ams\" << eor\n\n";
	out += "Task " + task + "\n";
	out += "System\n";
	out += "\tCharge " + std::to_string(charge) + "\n";
	out += "\tAtoms";


	for (int e = 0; e < s.set.size(); e++)
	{
		for (int atom = 0; atom < s.set[e].size(); atom++)
		{
			out += "\n\t\t" + s.elements[e];
			if (s.set[e][atom].x > 0) out += "  " + std::to_string(s.set[e][atom].x).substr(0, 7);
			else out += " " + std::to_string(s.set[e][atom].x).substr(0, 8);
			if (s.set[e][atom].y > 0) out += "  " + std::to_string(s.set[e][atom].y).substr(0, 7);
			else out += " " + std::to_string(s.set[e][atom].y).substr(0, 8);
			if (s.set[e][atom].z > 0) out += "  " + std::to_string(s.set[e][atom].z).substr(0, 7);
			else out += " " + std::to_string(s.set[e][atom].z).substr(0, 8);
		}
	}
	out += "\n\tEnd\nEnd\n";
	/*
	if (task == "GeometryOptimization")
	{
		out += "GeometryOptimization\n";
		out += "\toptim all\n";
		out += "\titerations 300\n";
		out += "\tstep rad=0.15 angle=10.0\n";
		out += "\thessupd BFGS\n";
		out += "\tConvergence Gradient=1.0e-4\n";
		out += "\tConvergence e=1.0e-4 grad=1.0e-4 rad=1.0e-2 angle=0.5\n";
		out += "End\n";
		out += "\t\n";
		out += "SCF\n";
		out += "\titerations 300\n";
		out += "\tconverge 1.0e-6 1.0e-6\n";
		out += "\tmixing 0.2\n";
		out += "\tlshift 0.0\n";
		out += "\tDIIS\n\t\tn 500\n\t\tok 0.00001\n\t\tcyc 500\n\tEnd\n";// cx = 5.0 cxx = 10.0\n";
		out += "END\n";
	}*/

	//Engine block
	out += "Engine ADF\n";//should be ADF
	out += "\tUnrestricted ";
	if (unrestricted) out += "Yes\n";
	else out += "No\n";
	out += "\tspinpolarization " + std::to_string(state) + "\n";
	out += "\tBasis\n\t\tType " + basis_name + "\n\t\tCore Large\n";
	out += "\tEnd\n";
	out += "\trelativity\n";
	out += "\t\tlevel scalar\n";
	out += "\t\tformalism ZORA\n";
	out += "\tEnd\n";
	out += "\tSCF\n";
	out += "\t\titerations " + std::to_string(maxSCF) + "\n";
	out += "\t\tconverge " + converge + " " + converge + "\n";
	out += "\t\tmixing 0.2\n";
	out += "\t\tlshift 0.0\n";
	out += "\t\tdiis\n\t\t\tn 500\n\t\t\tok 0.00001\n\t\t\tcyc 500\n\t\tEnd\n";// cx = 5.0 cxx = 10.0\n";
	out += "\tEND\n";

	out += "\txc\n";
	out += "\t\tgga " + method + "\n";
	out += "\tEnd\n";
	out += "EndEngine\n";
	out += "eor";
	std::ofstream file(filename + ".run");
	file << out;
	file.close();



}

structure* optimizationADF(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, std::string AMSHOME, int amscalculation, int maxSCF, std::string converge, double timeoutminutes, double timeLeftMinutes, bool keepFiles, bool unrestricted)
{

	if (amscalculation != -1)
	{
		//we are going to make separate folders for each ams

	}
	else
	{
		char* ac;
		std::string deleteFiles1 = "rm -r ams.results\n";
		std::string deleteFiles2 = "rm " + filename + ".run\n";
		std::string deleteFiles3 = "rm " + filename + ".out\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


		ac = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			ac[i] = deleteFiles2[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

		ac = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			ac[i] = deleteFiles3[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

	}
	/*

	*/
	std::string dirname = "calc" + std::to_string(amscalculation);
	char* ac;
	if (amscalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}

	writeToADF(s, "GeometryOptimization", dirname + "/" + filename, basis_name, method, maxSCF, charge, state, converge, unrestricted);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + filename + ".run");

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + filename + ".run");
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	//commands += "export AMSHOME=" + AMSHOME + "\n";//AMSHOME = /home/jpburkhardt/scm/ams2023.104
	//commands += "export AMSBIN=$AMSHOME/bin\nexport AMSRESOURCES=$AMSHOME/atomicdata\nexport SCMLICENSE=$AMSHOME/license.txt\n";
	//commands += "AMS_CALCULATION=" + filename + "\n" + "dos2unix $AMS_CALCULATION.run\n" + "chmod u+x $AMS_CALCULATION.run\n";
	//commands += "./$AMS_CALCULATION.run > " + filename + ".out\n";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + filename + ".run\n" + "chmod u+x " + filename + ".run\n";
	commands += "nohup ./" + filename + ".run > " + filename + ".out &\n";

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;


	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;


	//reading file input
	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + filename + ".out");


	//stall for cp2k to finish
	std::cout << "checking for guassian completion" << std::endl;



	while (!getline(*aOutFile, checker))
	{
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");
	}
	std::cout << "ADF started" << std::endl;

	//the file is writing complete?
	bool finishedOutput_energy = false;
	bool finishedOutput_geometry = false;
	std::vector<std::string> optStrings = {};
	std::string outLine;
	std::string energyString;

	bool timeOut = false;
	bool currentGeometry = false;//checks if there is at least some geometry optimized from input
	bool currentEnergy = false;
	std::string currentEnergyString = "";

	int natoms = 0;
	for (int i = 0; i < s.composition.size(); i++) natoms += s.composition[i];
	bool failed = false;
	while (!timeOut && (!finishedOutput_energy || !finishedOutput_geometry))
	{
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				//kill the task so we dont waste cpu resource
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				timeOut = true;
				//return nullptr;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				timeOut = true;
				//return nullptr;
			}

		}

		optStrings = {};
		//reopen file to get new updates
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");

		//look for the energy line

		while (getline(*aOutFile, outLine) && !failed) {

			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = "Energy (hartree)";
			std::string match2 = "Optimized geometry:";
			std::string match3 = "Geometry optimization failed";
			std::string match4 = "BAD CORE INTEGRAL";
			std::string match5 = "Process received SIGTERM";

			std::string someGeometry = "current energy";
			std::string someOptimization = "--------";
			std::string geometry = "Geometry";
			if (outLine.size() >= someGeometry.size())
			{
				bool matching = true;
				for (int m = 0; m < someGeometry.size(); m++)
				{
					if (outLine[m] != someGeometry[m]) matching = false;
				}
				if (matching) {
					currentEnergy = true;//cannot be made false
					currentEnergyString = outLine;
				}
			}
			bool intermediateGeomCheck = false;
			if (outLine.size() >= someOptimization.size())
			{
				bool matching = true;
				for (int m = 0; m < someOptimization.size(); m++)
				{
					if (outLine[m] != someOptimization[m]) matching = false;
				}
				if (matching) {
					intermediateGeomCheck = true;//cannot be made false
				}
			}
			if (intermediateGeomCheck)
			{
				bool getOpt = false;
				if (!getline(*aOutFile, outLine))
				{


					if (outLine.size() >= geometry.size())
					{
						bool matching = true;
						for (int m = 0; m < geometry.size(); m++)
						{
							if (outLine[m] != geometry[m]) matching = false;
						}
						if (matching) {
							getOpt = true;//cannot be made false
						}
					}
				}
				if (getOpt)
				{
					std::vector<std::string> tempOptStrings = {};

					bool stop = false;
					for (int i = 0; i < 4; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
					for (int e = 0; e < s.elements.size(); e++)  for (int a = 0; a < s.set[e].size(); a++) {
						if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
						tempOptStrings.push_back(outLine);
					}
					if (tempOptStrings.size() == natoms) {
						optStrings = tempOptStrings;
						currentGeometry = true;
					}
				}
			}



			if (outLine.size() >= match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) {
					finishedOutput_energy = true;//cannot be made false
					energyString = outLine;
				}
			}

			bool getopt = false;

			if (outLine.size() >= match2.size())
			{
				bool matching = true;
				for (int m = 0; m < match2.size(); m++)
				{
					if (outLine[m] != match2[m]) matching = false;
				}
				if (matching) {
					getopt = true;//cannot be made false
					finishedOutput_geometry = true;//geometry needs its own bool now because there are cases where energy calculates but geometry does not and crash line geometry optimization failed hasnt been printed yet
				}
			}
			if (getopt)
			{
				std::vector<std::string> tempOptStrings = {};
				bool stop = false;
				for (int i = 0; i < 7; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
				for (int e = 0; e < s.elements.size(); e++)  for (int a = 0; a < s.set[e].size(); a++) {
					if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
					tempOptStrings.push_back(outLine);
				}
				if (tempOptStrings.size() == natoms) optStrings = tempOptStrings;
			}
			//when the AMS calculation is successful, energy is always printed after the geometry, and therefore its fine to get the geometry here.

			failed = false;
			if (outLine.size() >= match3.size())
			{
				bool matching = true;
				for (int m = 0; m < match3.size(); m++)
				{
					if (outLine[m] != match3[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (outLine.size() >= match4.size())
			{
				bool matching = true;
				for (int m = 0; m < match4.size(); m++)
				{
					if (outLine[m] != match4[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (outLine.size() >= match5.size())
			{
				bool matching = true;
				for (int m = 0; m < match5.size(); m++)
				{
					if (outLine[m] != match5[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			/*if (failed)
			{
				return nullptr;
			}*/
		}



	}
	std::cout << "finished output: " << finishedOutput_energy << std::endl;

	if (timeOut && (!currentEnergy || !currentGeometry))
	{
		failed = true;
		//return nullptr;
	}
	if (failed && optStrings.size() != natoms)
	{
		//got no structure; complete failure
		return nullptr;
	}

	double seconds_since_start = difftime(time(0), start);

	std::string getout = "cd ../\n";
	if (amscalculation != -1) {
		ac = new char[getout.length() + 1];
		for (int i = 0; i < getout.length(); i++)
		{
			ac[i] = getout[i];
		}

		std::system(ac);
		std::cout << "command: " << getout << std::endl;
		delete[] ac;
	}
	std::cout << "ADF time: " << seconds_since_start << std::endl;

	structure* st = nullptr;



	std::vector<std::string> elements = {};
	std::vector<std::vector<atom>> set = {};

	for (int i = 0; i < optStrings.size(); i++)
	{
		std::cout << optStrings[i] << std::endl;

		std::string cline = optStrings[i];

		std::string element = "";
		if (cline[12] == ' ') element = cline.substr(13, 1);//it goes the other way ....
		else element = cline.substr(12, 2);
		//x
		double x;
		double y;
		double z;
		bool sign = false;
		if (cline[18] == '-') sign = true;
		if (sign) x = std::stod(cline.substr(18, 11));
		else x = std::stod(cline.substr(19, 10));
		if (cline[33] == '-') sign = true;
		else sign = false;
		if (sign) y = std::stod(cline.substr(33, 11));
		else y = std::stod(cline.substr(34, 10));
		if (cline[48] == '-') sign = true;
		else sign = false;
		if (sign) z = std::stod(cline.substr(48, 11));
		else z = std::stod(cline.substr(49, 10));



		//figure out which set this belongs to.
		int index = -1;
		for (int i = 0; i < elements.size(); i++)
		{
			if (elements[i] == element)
			{
				index = i;
			}
		}
		//if (buggo) std::cout << "making atom for elemental index " << index << " with coords " << x << " " << y << " " << z << std::endl;
		if (index == -1)
		{
			//if (buggo) std::cout << "not yet indexed" << std::endl;
			//atom type not yet indexed
			elements.push_back(element);
			std::vector < atom> alist = { atom(x,y,z) };
			set.push_back(alist);
		}
		else {
			//if (buggo) std::cout << "pushing back" << std::endl;
			set[index].push_back(atom(x, y, z));
		}
	}
	//if (buggo) std::cout << "this structure is not the last one starting over with next" << std::endl;

	if (st != nullptr) delete st;
	st = new structure(set, elements);


	std::cout << "found energy" << std::endl;

	std::string matchString;
	if (!timeOut) matchString = "Energy (hartree)";
	else {
		matchString = "current energy";
		energyString = currentEnergyString;
	}
	std::string num = "";
	for (int i = matchString.length(); i < energyString.length(); i++)
	{
		if (energyString[i] != ' ') num += energyString[i];
	}
	double energy = std::stod(num);
	st->energy = energy;
	// Close the file
	//if there was no convergence then nullptr will be returned
	//if (buggo) std::cout << "closing" << std::endl;
	aOutFile->close();
	delete aOutFile;


	//remove files
	if (!keepFiles)
	{
		if (amscalculation == -1)
		{

			std::string deleteFiles1 = "rm -r ams.results\n";
			std::string deleteFiles2 = "rm " + filename + ".run\n";
			std::string deleteFiles3 = "rm " + filename + ".out\n";

			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;


			ac = new char[deleteFiles2.length() + 1];
			for (int i = 0; i < deleteFiles2.length(); i++)
			{
				ac[i] = deleteFiles2[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;

			ac = new char[deleteFiles3.length() + 1];
			for (int i = 0; i < deleteFiles3.length(); i++)
			{
				ac[i] = deleteFiles3[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;
		}
		else {
			std::string deleteFiles1 = "rm -r " + dirname + "\n";
			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;

		}

	}

	return st;
	//return energyHartree;
}

double energyADF(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, std::string AMSHOME, int adfcalculation, std::string converge, int timeoutminutes, int timeLeftMinutes, bool keepFiles, bool unrestricted)
{
	if (adfcalculation == -1)//always delete with ams, make sure its gone to prevent crashes
	{
		char* ac;
		std::string deleteFiles1 = "rm -r ams.results\n";
		std::string deleteFiles2 = "rm " + filename + ".run\n";
		std::string deleteFiles3 = "rm " + filename + ".out\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


		ac = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			ac[i] = deleteFiles2[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

		ac = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			ac[i] = deleteFiles3[i];
		}
		std::system(ac);
		std::cout << "command: " << ac << std::endl;
		delete[] ac;

	}
	std::string dirname = "calc" + std::to_string(adfcalculation);
	char* ac;
	if (adfcalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\n";//cd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;

	}

	writeToADF(s, "SinglePoint", dirname + "/" + filename, basis_name, method, 300, charge, state, converge, unrestricted);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + filename + ".run");

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + filename + ".run");
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	//commands += "export AMSHOME=" + AMSHOME + "\n";//AMSHOME = /home/jpburkhardt/scm/ams2023.104
	//commands += "export AMSBIN=$AMSHOME/bin\nexport AMSRESOURCES=$AMSHOME/atomicdata\nexport SCMLICENSE=$AMSHOME/license.txt\n";
	//commands += "AMS_CALCULATION=" + filename + "\n";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + filename + ".run\n" + "chmod u+x " + filename + ".run\n";
	commands += "nohup ./" + filename + ".run > " + filename + ".out &\n";

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;


	//reading file input
	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + filename + ".out");


	//stall for cp2k to finish
	std::cout << "checking for adf completion" << std::endl;




	while (!getline(*aOutFile, checker))
	{
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");
	}
	std::cout << "ADF started" << std::endl;

	//the file is writing complete?
	bool finishedOutput = false;
	std::string outLine;
	std::string energyString;
	while (!finishedOutput)
	{
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				//kill the task so we dont waste cpu resource
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill ams\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;
			}

		}


		//reopen file to get new updates
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + filename + ".out");

		//look for the energy line

		while (getline(*aOutFile, outLine)) {
			//come back here plz
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = "Energy (hartree)";
			std::string match4 = "BAD CORE INTEGRAL";
			std::string match5 = "Process received SIGTERM";
			if (outLine.size() >= match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) {
					finishedOutput = true;//cannot be made false
					energyString = outLine;
				}
			}
			bool failed = false;
			if (outLine.size() >= match4.size())
			{
				bool matching = true;
				for (int m = 0; m < match4.size(); m++)
				{
					if (outLine[m] != match4[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (outLine.size() >= match5.size())
			{
				bool matching = true;
				for (int m = 0; m < match5.size(); m++)
				{
					if (outLine[m] != match5[m]) matching = false;
				}
				if (matching) failed = true;//cannot be made false
			}
			if (failed)
			{
				return DBL_MAX;
			}
		}



	}
	std::cout << "finished output: " << finishedOutput << std::endl;


	double seconds_since_start = difftime(time(0), start);
	if (adfcalculation != -1) {
		std::string getout = "cd ../\n";

		ac = new char[getout.length() + 1];
		for (int i = 0; i < getout.length(); i++)
		{
			ac[i] = getout[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}
	std::cout << "ADF time: " << seconds_since_start << std::endl;


	std::cout << "found energy" << std::endl;
	std::string matchString = "Energy (hartree)";
	std::string num = "";
	for (int i = matchString.length(); i < energyString.length(); i++)
	{
		if (energyString[i] != ' ') num += energyString[i];
	}
	double energy = std::stod(num);
	// Close the file
	//if there was no convergence then nullptr will be returned
	//if (buggo) std::cout << "closing" << std::endl;
	aOutFile->close();
	delete aOutFile;


	//remove files
	if (!keepFiles)
	{
		if (adfcalculation == -1)
		{

			std::string deleteFiles1 = "rm -r ams.results\n";
			std::string deleteFiles2 = "rm " + filename + ".run\n";
			std::string deleteFiles3 = "rm " + filename + ".out\n";

			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;


			ac = new char[deleteFiles2.length() + 1];
			for (int i = 0; i < deleteFiles2.length(); i++)
			{
				ac[i] = deleteFiles2[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;

			ac = new char[deleteFiles3.length() + 1];
			for (int i = 0; i < deleteFiles3.length(); i++)
			{
				ac[i] = deleteFiles3[i];
			}
			std::system(ac);
			std::cout << "command: " << ac << std::endl;
			delete[] ac;
		}
		else {
			std::string deleteFiles1 = "rm -r " + dirname + "\n";
			ac = new char[deleteFiles1.length() + 1];
			for (int i = 0; i < deleteFiles1.length(); i++)
			{
				ac[i] = deleteFiles1[i];
			}
			std::cout << "command: " << ac << std::endl;
			std::system(ac);
			delete[] ac;

		}

	}

	return energy;
	//return energyHartree;
}

double basinHoppingEnergy(structure& s, int ittt, std::string task, double a, double b, double c)
{
	if (DEBUG > 1)
	{
		//scannign structure
		for (std::vector<std::vector<atom>>::iterator it = s.set.begin(); it != s.set.end(); it++)
		{
			std::cout << "element set" << std::endl;
			for (int i = 0; i < it->size(); i++)
			{
				std::cout << "atom: " << (*it)[i].x << " " << (*it)[i].y << " " << (*it)[i].z << std::endl;
			}
		}
		std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
		std::string inFileName = fileName + "_in.txt";
		writeToCP2K(s, task, inFileName, "/home/jpburkhardt/GMS/data/BASIS_SET", "/home/jpburkhardt/GMS/data/GTH_POTENTIALS", a, b, c);

	}

	/*
	New basin hopping energy iterator
	*/
	double energy = 0;
	double sigma = 1;
	double epsilon = 1;
	if (DEBUG > 1) std::cout << "starting energy calculation" << std::endl;

	//the iterators i made for structure only work with do while loops
	structure::Iterator it = s.begin();
	do {
		if (DEBUG > 3) std::cout << "new jt loop about to begin" << std::endl;
		structure::Iterator jt = s.begin();
		do {
			if (DEBUG > 3) std::cout << "comparing" << "pointer " << it.m_ptr << " to " << "pointer " << jt.m_ptr << std::endl;
			if (it != jt)
			{
				//if (DEBUG) std::cout << "doing some math" << std::endl;

				if (DEBUG > 2) std::cout << "atom " << it.m_ptr << ": " << it->x << " " << it->y << " " << it->z << " vs atom2 " << jt.m_ptr << ": " << jt->x << " " << jt->y << " " << jt->z << std::endl;
				double r = euclideanDistance((*it), (*jt));
				//if (DEBUG) std::cout << "is r 0: " << r << std::endl;
				//std::cout << "r: " << r << std::endl;
				double e = pow(sigma / r, 12) - pow(sigma / r, 6);
				if (DEBUG > 1) std::cout << it.indexI << "," << it.indexJ << " vs. " << jt.indexI << "," << jt.indexJ << " e: " << e << std::endl;
				energy += e;
			}
			else {
				if (DEBUG > 1) std::cout << "did not calculate an atom against itself" << std::endl;
			}
			if (DEBUG > 3) std::cout << "about to iterate jt" << std::endl;
			jt++;
		} while (jt != s.end());

		it++;
	} while (it != s.end());



	/*
	for (std::vector<std::vector<atom>>::iterator it = s.set.begin(); it != s.set.end(); it++)
	{
		if (DEBUG) std::cout << "iterating outer set" << std::endl;
		for (int i = 0; i < it->size(); i++)
		{
			if (DEBUG) std::cout << "inner loop 1: all inner sets" << std::endl;
			for (std::vector<std::vector<atom>>::iterator itt = s.set.begin(); itt != s.set.end(); itt++)
			{
				if (DEBUG) std::cout << "inner loop 2: inner sets size: " << itt->size() << std::endl;
				for (int j = 0; j < itt->size(); j++)
				{
					if (DEBUG) std::cout << "about to do math" << std::endl;
					if (&(*it)[i] != &(*itt)[j])
					{
						if (DEBUG) std::cout << "doing some math" << std::endl;
						if (DEBUG) std::cout << "it size:" << it->size() << " " << i << std::endl;
						if (DEBUG) std::cout << "itt size:" << itt->size() << " " << j << std::endl;
						double r = euclideanDistance((*it)[i], (*itt)[j]);
						//if (DEBUG) std::cout << "is r 0: " << r << std::endl;
						//std::cout << "r: " << r << std::endl;
						double e = pow(sigma / r, 12) - pow(sigma / r, 6);
						//std::cout << "e: " << e << std::endl;
						energy += e;
					}
					else {
						if (DEBUG) std::cout << "did not calculate duplicate" << std::endl;
					}
					if (DEBUG) std::cout << "math done" << std::endl;
				}
				if (DEBUG) std::cout << "here 2" << std::endl;
			}
			if (DEBUG) std::cout << "here 1" << std::endl;
		}
		if (DEBUG) std::cout << "outer loop complete" << std::endl;
	}startModifiedBasinHopping(init, x, y, z, percentChangeI, percentChangeF, composition, elements, charge, state, dft, method, basis, optimizationCycles, coordinationSteps, HFsteps, time, localMinimaTrappingSteps, radialCriteriaPercent, outputFileName, taskName, a, b, c, blindGradientDescent, directionExponent, proportionExponent, deltaEThr);

	*/
	if (DEBUG > 1) std::cout << "returning" << std::endl;
	return 0.5 * epsilon * energy;
};

structure* optimizationBasinHopping(structure& s, int ittt, std::string task, int maxCycles, double x, double y, double z)
{
	//perform max cylces random steps
	std::vector<double> movementVector = {};
	double lowestEnergy = s.energy;
	structure* st = nullptr;
	for (int i = 0; i < s.elements.size(); i++)
	{
		for (int j = 0; j < s.set[i].size(); j++)
		{
			movementVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * x));
			movementVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * y));
			movementVector.push_back(((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * z));
		}

	}
	for (int i = 0; i < maxCycles; i++)
	{
		for (int i = 0; i < s.elements.size(); i++)
		{
			for (int j = 0; j < s.set[i].size(); j++)
			{
				movementVector[i * 3] = (((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * x));
				movementVector[i * 3 + 1] = (((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * y));
				movementVector[i * 3 + 2] = (((2.0 * (rand() % RAND_MAX) / RAND_MAX - 1) * z));
			}

		}
		structure* sn = nullptr;
		if (st != nullptr) sn = new structure(*st, movementVector);
		else sn = new structure(s, movementVector);
		double snergy = basinHoppingEnergy(*sn, 0, "task", x, y, z);
		sn->energy = snergy;
		if (snergy < lowestEnergy)
		{
			if (st != nullptr) delete st;
			st = sn;
			lowestEnergy = snergy;
		}
	}
	return st;
}



double energyCP2K(structure& s, int ittt, std::string task, double a, double b, double c)
{

	//these a b c are for cell shape
	std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
	std::string inFileName = fileName + "_in.txt";
	std::string outFileName = fileName + "_out.txt";
	writeToCP2K(s, task, inFileName, "/home/jpburkhardt/GMS/data/BASIS_SET", "/home/jpburkhardt/GMS/data/GTH_POTENTIALS", a, b, c);

	//ensure file has been written
	std::string checker;
	std::ifstream* cp2kInput = new std::ifstream(inFileName);

	while (!getline(*cp2kInput, checker))
	{
		cp2kInput->close();
		delete cp2kInput;
		cp2kInput = new std::ifstream(inFileName);
	}
	cp2kInput->close();
	delete cp2kInput;
	std::cout << "input file written" << std::endl;

	std::string cp2kcommand = "cp2k.psmp -i " + inFileName + " - o " + outFileName;

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	char* cp2kc = new char[cp2kcommand.length() + 1];
	for (int i = 0; i < cp2kcommand.length(); i++)
	{
		cp2kc[i] = cp2kcommand[i];
	}
	cp2kc[cp2kcommand.length()] = '\0';
	std::cout << "command: " << cp2kc << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "command: " << cp2kc << std::endl;
	std::system(cp2kc);
	std::cout << "command written" << std::endl;

	std::ifstream* cp2koutFile = new std::ifstream(outFileName + ".txt");


	//stall for cp2k to finish
	std::cout << "checking for cp2k completion" << std::endl;
	time_t start = time(0);


	while (!getline(*cp2koutFile, checker))
	{
		cp2koutFile->close();
		delete cp2koutFile;
		cp2koutFile = new std::ifstream(outFileName);
	}
	std::cout << "cp2k started" << std::endl;
	//the file is located but energy may not be calculated yet
	bool energyCalculated = false;
	std::string energyLine;
	while (!energyCalculated)
	{
		//reopen file to get new updates
		cp2koutFile->close();
		delete cp2koutFile;
		cp2koutFile = new std::ifstream(outFileName);

		//look for the energy line
		while (getline(*cp2koutFile, energyLine)) {
			std::string match = "ENERGY| Total";
			if (energyLine.size() > 6)
			{
				if (energyLine[1] == 'E' && energyLine[2] == 'N' && energyLine[3] == 'E' && energyLine[4] == 'R' && energyLine[5] == 'G' && energyLine[6] == 'Y')
				{
					energyCalculated = true;
				}
			}
		}
	}

	double seconds_since_start = difftime(time(0), start);
	std::cout << "cp2k time: " << seconds_since_start << std::endl;
	std::cout << "reading file" << std::endl;

	double value = 0.0;
	// Use a while loop together with the getline() function to read the file line by line

	double enrg = 0.0;
	int count = 0;
	int j = energyLine.size() - 1;
	while (energyLine[j - 1] != ' ')
	{
		j--;
	}

	std::string val;
	while (j != energyLine.size())
	{
		val += energyLine[j];
		j++;
	}
	std::cout << val << std::endl;
	value = std::stod(val);


	// Close the file
	cp2koutFile->close();
	delete cp2koutFile;
	std::cout << "energy: " << value << std::endl;
	return value;
}


double energyGuassian(structure& s, int ittt, std::string task, std::string basis, std::string method, int charge, int state, int maxSCF, bool keepFiles, std::string customBasis)
{
	//outfile name should end in log, not sure what the input should end in

	//these a b c are for cell shape
	std::string guassianTaskName = task + " energy";

	std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
	std::string inFileName = fileName + ".com";
	std::string outFileName = fileName + ".log";

	writeToGuassian(s, guassianTaskName, inFileName, basis, method, charge, state, 0, maxSCF, customBasis);

	//ensure file has been written
	std::string checker;
	std::ifstream* guassianInput = new std::ifstream(inFileName);

	while (!getline(*guassianInput, checker))
	{
		guassianInput->close();
		delete guassianInput;
		guassianInput = new std::ifstream(inFileName);
	}
	guassianInput->close();
	delete guassianInput;
	std::cout << "input file written" << std::endl;

	std::string guassiancommand = "g16 <" + inFileName + " >" + outFileName;

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
	//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	char* gc = new char[guassiancommand.length() + 1];
	for (int i = 0; i < guassiancommand.length(); i++)
	{
		gc[i] = guassiancommand[i];
	}
	gc[guassiancommand.length()] = '\0';
	std::cout << "command: " << gc << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "command: " << gc << std::endl;
	time_t start = time(0);
	std::system(gc);
	delete[]gc;
	//std::cout << "command written" << std::endl;

	std::ifstream* gOutFile = new std::ifstream(outFileName);


	//stall for cp2k to finish
	std::cout << "checking for guassian completion" << std::endl;



	while (!getline(*gOutFile, checker))
	{
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);
	}
	std::cout << "guassian started" << std::endl;
	//the file is located but energy may not be calculated yet
	bool finishedOutput = false;

	std::string outLine;
	while (!finishedOutput)
	{
		//reopen file to get new updates
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);

		//look for the energy line
		while (getline(*gOutFile, outLine)) {
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = " Job cpu time:";

			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
			}
		}

	}
	double seconds_since_start = difftime(time(0), start);
	std::cout << "gaussian time: " << seconds_since_start << std::endl;

	std::string energyLine;
	int lastLines = 0;
	int clines = 0;
	//reopen file to get new updates
	gOutFile->close();
	delete gOutFile;
	gOutFile = new std::ifstream(outFileName);
	bool energyCalculated = false;
	//look for the energy line
	while (!energyCalculated) {
		//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line
		if (getline(*gOutFile, energyLine))
		{
			std::string match = " SCF Done:  E";

			if (energyLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (energyLine[m] != match[m]) matching = false;
				}
				energyCalculated = matching;
			}
		}
		else {
			break;
		}

	}
	if (!energyCalculated) return DBL_MAX;//this might mess with gradient descent


	std::cout << "reading file" << std::endl;

	//find the start of the energy
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));




	// Close the file
	gOutFile->close();
	delete gOutFile;
	std::cout << "energy: " << energyHartree << std::endl;

	//remove files
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm " + inFileName + "\n";
		std::string deleteFiles2 = "rm " + outFileName + "\n";
		std::string deleteFiles3 = "rm -r core.*\n";

		gc = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			gc[i] = deleteFiles1[i];
		}
		std::cout << "command: " << gc << std::endl;
		std::system(gc);
		delete[] gc;


		gc = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			gc[i] = deleteFiles2[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

		gc = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			gc[i] = deleteFiles3[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

	}

	return energyHartree;
}

structure* optimizationGaussian(structure& s, int ittt, std::string task, std::string basis, std::string method, int maxCycles, int charge, int state, int maxSCF, bool keepFiles, std::string customBasis)
{
	//the method should be HF instead of DFT for better convergence
	//outfile name should end in log, not sure what the input should end in


	std::string guassianTaskName = task + " opt";

	std::string fileName = task + "_" + std::to_string(ittt) + "_" + std::to_string(rand() / 10);
	std::string inFileName = fileName + ".com";
	std::string outFileName = fileName + ".log";

	writeToGuassian(s, guassianTaskName, inFileName, basis, method, charge, state, maxCycles, maxSCF, customBasis);

	//ensure file has been written
	std::string checker;
	std::ifstream* guassianInput = new std::ifstream(inFileName);

	while (!getline(*guassianInput, checker))
	{
		guassianInput->close();
		delete guassianInput;
		guassianInput = new std::ifstream(inFileName);
	}
	guassianInput->close();
	delete guassianInput;
	std::cout << "input file written" << std::endl;

	std::string guassiancommand = "g16 <" + inFileName + " >" + outFileName;

	//"cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i testingCH4_1_929_in.txt -o " + outFileName + ".txt";
		//std::string cp2kcommand = "cd \"\"" + cp2kpath + "\"\"\nstartcp2k -i " + inFileName + ".txt -o " + outFileName + ".txt\n";
	std::cout << "printing input" << std::endl;

	char* gc = new char[guassiancommand.length() + 1];
	for (int i = 0; i < guassiancommand.length(); i++)
	{
		gc[i] = guassiancommand[i];
	}
	gc[guassiancommand.length()] = '\0';
	std::cout << "command: " << gc << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "command: " << gc << std::endl;
	time_t start = time(0);
	std::system(gc);
	delete[] gc;
	//std::cout << "command written" << std::endl;

	std::ifstream* gOutFile = new std::ifstream(outFileName);


	//stall for cp2k to finish
	std::cout << "checking for guassian completion" << std::endl;



	while (!getline(*gOutFile, checker))
	{
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);
	}
	std::cout << "guassian started" << std::endl;

	//the file is writing complete?
	bool finishedOutput = false;

	std::string outLine;
	while (!finishedOutput)
	{
		//reopen file to get new updates
		gOutFile->close();
		delete gOutFile;
		gOutFile = new std::ifstream(outFileName);

		//look for the energy line
		while (getline(*gOutFile, outLine)) {
			//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

			std::string match = " Job cpu time:";

			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
			}
		}

	}


	double seconds_since_start = difftime(time(0), start);
	std::cout << "gaussian time: " << seconds_since_start << std::endl;
	std::cout << "reading file" << std::endl;

	//find the start of the optimization
	//this will take the orientation everytime but iguess thats ok
	//restart the file just incase
	gOutFile->close();
	delete gOutFile;
	gOutFile = new std::ifstream(outFileName);
	std::string cline;
	structure* st = nullptr;

	//bool buggo = true;
	//if (buggo) std::cout << "opening line one" << std::endl;
	std::string matchString = "Input orientation:";
	std::string eMatch = " SCF Done:  E";
	std::string energyLine = "";



	std::cout << "reading file" << std::endl;

	//find the start of the energy

	while (getline(*gOutFile, cline)) {
		//if (buggo) std::cout << "line one opened" << std::endl;
		if (cline.length() > matchString.length())
		{

			bool matching = false;
			int matchIt = 0;
			//if (buggo) std::cout << "checking matches" << std::endl;

			for (int l = 0; l < cline.length() - matchString.length(); l++)
			{
				//if (buggo) std::cout << "l=" << l << " cline length: " << cline.length() << " , length - match: " << cline.length() - matchString.length() << " match string: " << matchString << std::endl;
				if (cline[l] == matchString[matchIt])
				{
					//if (buggo) std::cout << "match" << std::endl;
					matchIt++;
					matching = true;
				}
				else if (matchIt < 18) {
					//if (buggo) std::cout << "no match" << std::endl;
					//if (matchIt > 10 && buggo) std::cout << "closed early on: " << cline << " with a matchit of " << matchIt << " and a l of " << l << " and a cline length of " << cline.length() << " and a matchit length of " << matchString.length() << std::endl;
					matchIt = 0;
					matching = false;
				}
				else {
					//the strings are already matched let the match go through
				}
			}
			if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
			if (matching && matchIt >= matchString.length())
			{
				//if (buggo) std::cout << "matched" << std::endl;
				//a geometry has been identified

				//create set of atoms for those inside
				std::vector<std::vector<atom>> set;
				std::vector<std::string> elements;
				//the next 4 lines do not contain atom info
				//if (buggo) std::cout << "skipping 4 lines" << std::endl;
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				char space = ' ';


				//if (buggo) std::cout << "reading an atom" << std::endl;

				while (getline(*gOutFile, cline) && cline[1] != '-')
				{
					std::string cnumber = "";

					//if (buggo) std::cout << cline << std::endl;
					int num = 0;
					std::string element;
					double x = 0;
					double y = 0;
					double z = 0;
					std::cout << "cline: " << cline << std::endl;
					for (int i = 0; i < cline.length(); i++)
					{
						if (cline[i] == space)
						{
							std::cout << "current cnumber: " << cnumber << std::endl;
							if (cnumber.length() > 0)
							{
								//if (buggo) std::cout << "cnumber larger than one. it is " << cnumber << std::endl;
								num++;
								switch (num)
								{
								case 1:
									//center number, not sure what this is for. do nothing
									break;
								case 2:
									//atomic number
									element = atomicNumber(std::stoi(cnumber));
									std::cout << "elemenent: " << element << std::endl;
									break;
								case 3:
									//atomic type, not sure what this is. do nothing.
									break;
								case 4:
									//x coordinate
									std::cout << " x: " << x << std::endl;
									x = std::stof(cnumber);
									break;
								case 5:
									//y
									std::cout << " y: " << y << std::endl;
									y = std::stof(cnumber);
									break;
								case 6:
									//z
									std::cout << " z: " << z << std::endl;
									z = std::stof(cnumber);
									break;
								default:
									break;

								}
								cnumber = "";
							}
							//else {
								//if (buggo) std::cout << "current char not large than one" << std::endl;
							//}

						}
						else {
							//std::cout << "adding " << std::to_string(cline[i]) << " or " << cline[i] << "? to cnumber" << std::endl;
							cnumber += cline[i];
						}
					}
					//figure out which set this belongs to.
					int index = -1;
					for (int i = 0; i < elements.size(); i++)
					{
						if (elements[i] == element)
						{
							index = i;
						}
					}
					//if (buggo) std::cout << "making atom for elemental index " << index << " with coords " << x << " " << y << " " << z << std::endl;
					if (index == -1)
					{
						//if (buggo) std::cout << "not yet indexed" << std::endl;
						//atom type not yet indexed
						elements.push_back(element);
						std::vector < atom> alist = { atom(x,y,z) };
						set.push_back(alist);
					}
					else {
						//if (buggo) std::cout << "pushing back" << std::endl;
						set[index].push_back(atom(x, y, z));
					}
				}
				//if (buggo) std::cout << "this structure is not the last one starting over with next" << std::endl;

				if (st != nullptr) delete st;
				st = new structure(set, elements);


			}
			//else if (buggo) std::cout << "no match" << std::endl;

		}
		if (cline.size() > eMatch.size())
		{
			bool matching = true;
			for (int m = 0; m < eMatch.size(); m++)
			{
				if (cline[m] != eMatch[m]) matching = false;
			}
			if (matching) energyLine = cline;
		}
		//structure(std::vector<std::vector<atom>> set, std::vector<std::string> elements, bool sort = false)

	}
	// Close the file
	//if there was no convergence then nullptr will be returned
	//if (buggo) std::cout << "closing" << std::endl;
	gOutFile->close();
	delete gOutFile;
	////
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));
	st->energy = energyHartree;
	/////
	//remove files
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm " + inFileName + "\n";
		std::string deleteFiles2 = "rm " + outFileName + "\n";
		std::string deleteFiles3 = "rm -r core.*\n";

		gc = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			gc[i] = deleteFiles1[i];
		}
		std::cout << "command: " << gc << std::endl;
		std::system(gc);
		delete[] gc;


		gc = new char[deleteFiles2.length() + 1];
		for (int i = 0; i < deleteFiles2.length(); i++)
		{
			gc[i] = deleteFiles2[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

		gc = new char[deleteFiles3.length() + 1];
		for (int i = 0; i < deleteFiles3.length(); i++)
		{
			gc[i] = deleteFiles3[i];
		}
		std::system(gc);
		std::cout << "command: " << gc << std::endl;
		delete[] gc;

	}

	return st;
}

structure* newOptimizationGaussian(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, int amscalculation, int maxSCF, int maxCycles, double timeoutminutes, double timeLeftMinutes, bool keepFiles, std::string customBasis)
{
	std::string dirname = "calc" + std::to_string(amscalculation);
	char* ac;
	if (amscalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}

	std::string guassianTaskName = filename + " opt";

	std::string inFileName = filename + ".com";
	std::string outFileName = filename + ".log";

	writeToGuassian(s, guassianTaskName, dirname + "/" + inFileName, basis_name, method, charge, state, maxCycles, maxSCF, customBasis);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + inFileName);

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + inFileName);
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + inFileName + "\n";
	commands += "nohup g16 <" + inFileName + " >" + outFileName + " &\n";
	std::cout << "printing input" << std::endl;


	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;

	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + outFileName);


	std::cout << "checking for gaussian completion" << std::endl;



	bool failedToStart = false;
	while (!getline(*aOutFile, checker) && !failedToStart)
	{

		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);

		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				failedToStart = true;
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return nullptr;

			}

		}
		if (timeoutminutes != -1)
		{
			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				failedToStart = true;
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return nullptr;
			}

		}
	}
	std::cout << "Gaussian started" << std::endl;

	bool finishedOutput = false;
	bool finishedEnergy = false;
	bool finishedStructure = false;
	std::vector<std::string> optStrings = {};
	std::string outLine;
	std::string energyString;

	bool timedOut = false;


	while (!timedOut && (!finishedOutput || !finishedEnergy || !finishedStructure))
	{
		if (finishedOutput && !(finishedEnergy && finishedStructure))
		{
			commands = "pkill g16\n";
			ac = new char[commands.length() + 1];
			for (int i = 0; i < commands.length(); i++)
			{
				ac[i] = commands[i];
			}
			ac[commands.length()] = '\0';
			std::cout << "commands: " << ac << std::endl;
			std::cout << "writing command" << std::endl;
			std::system(ac);
			delete[] ac;
			//timedOut = true;
			return nullptr;

		}
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				timedOut = true;
				//return nullptr;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				timedOut = true;
				//return nullptr;
			}

		}

		optStrings = {};
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);

		std::string match = " Job cpu time:";
		std::string matchString = "Input orientation:";
		std::string eMatch = " SCF Done:  E";
		while (getline(*aOutFile, outLine)) {
			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
			}
			if (outLine.size() > matchString.size())
			{
				bool matching = true;
				for (int m = 0; m < matchString.size(); m++)
				{
					if (outLine[m] != matchString[m]) matching = false;
				}
				if (matching) finishedStructure = true;//cannot be made false
			}
			if (outLine.size() > eMatch.size())
			{
				bool matching = true;
				for (int m = 0; m < eMatch.size(); m++)
				{
					if (outLine[m] != eMatch[m]) matching = false;
				}
				if (matching) finishedEnergy = true;//cannot be made false
			}
		}





	}
	std::cout << "finished output: " << finishedOutput << std::endl;
	double seconds_since_start = difftime(time(0), start);

	if (!finishedEnergy || !finishedOutput)
	{
		return nullptr;//otherwise we at least can get an energy and geometry
	}

	aOutFile->close();
	delete aOutFile;
	std::ifstream* gOutFile = new std::ifstream(dirname + "/" + outFileName);
	std::string cline;
	structure* st = nullptr;

	std::string matchString = "Input orientation:";
	std::string eMatch = " SCF Done:  E";
	std::string energyLine = "";

	std::cout << "reading file" << std::endl;

	while (getline(*gOutFile, cline)) {
		if (cline.length() > matchString.length())
		{

			bool matching = false;
			int matchIt = 0;
			for (int l = 0; l < cline.length() - matchString.length(); l++)
			{
				if (cline[l] == matchString[matchIt])
				{
					matchIt++;
					matching = true;
				}
				else if (matchIt < 18) {
					matchIt = 0;
					matching = false;
				}
				else {
					//the strings are already matched let the match go through
				}
			}
			if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
			if (matching && matchIt >= matchString.length())
			{
				std::vector<std::vector<atom>> set;
				std::vector<std::string> elements;
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				getline(*gOutFile, cline);
				char space = ' ';
				while (getline(*gOutFile, cline) && cline[1] != '-')
				{
					std::string cnumber = "";
					int num = 0;
					std::string element;
					double x = 0;
					double y = 0;
					double z = 0;
					for (int i = 0; i < cline.length(); i++)
					{
						if (cline[i] == space)
						{
							if (cnumber.length() > 0)
							{
								num++;
								switch (num)
								{
								case 1:
									//center number, not sure what this is for. do nothing
									break;
								case 2:
									//atomic number
									element = atomicNumber(std::stoi(cnumber));
									break;
								case 3:
									//atomic type, not sure what this is. do nothing.
									break;
								case 4:
									//x coordinate
									x = std::stof(cnumber);
									break;
								case 5:
									//y
									y = std::stof(cnumber);
									break;
								case 6:
									//z
									z = std::stof(cnumber);
									break;
								default:
									break;

								}
								cnumber = "";
							}
						}
						else {
							cnumber += cline[i];
						}
					}
					int index = -1;
					for (int i = 0; i < elements.size(); i++)
					{
						if (elements[i] == element)
						{
							index = i;
						}
					}
					if (index == -1)
					{
						elements.push_back(element);
						std::vector < atom> alist = { atom(x,y,z) };
						set.push_back(alist);
					}
					else {
						set[index].push_back(atom(x, y, z));
					}
				}
				if (st != nullptr) delete st;
				st = new structure(set, elements);


			}
		}
		if (cline.size() > eMatch.size())
		{
			bool matching = true;
			for (int m = 0; m < eMatch.size(); m++)
			{
				if (cline[m] != eMatch[m]) matching = false;
			}
			if (matching) energyLine = cline;
		}
	}
	gOutFile->close();
	delete gOutFile;

	//find the start of the energy
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));
	st->energy = energyHartree;

	std::cout << "Gaussian time: " << seconds_since_start << std::endl;
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm -r " + dirname + "\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


	}
	return st;
}
double newEnergyGaussian(structure& s, std::string filename, std::string basis_name, std::string method, int charge, int state, int amscalculation, int maxSCF, int maxCycles, double timeoutminutes, double timeLeftMinutes, bool keepFiles, std::string customBasis)
{
	std::string dirname = "calc" + std::to_string(amscalculation);
	char* ac;
	if (amscalculation != -1) {
		std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";

		ac = new char[mkdir.length() + 1];
		for (int i = 0; i < mkdir.length(); i++)
		{
			ac[i] = mkdir[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;
	}

	std::string guassianTaskName = filename + " energy";

	std::string inFileName = filename + ".com";
	std::string outFileName = filename + ".log";

	writeToGuassian(s, guassianTaskName, dirname + "/" + inFileName, basis_name, method, charge, state, 0, maxSCF, customBasis);

	std::string checker;
	std::ifstream* adfInput = new std::ifstream(dirname + "/" + inFileName);

	while (!getline(*adfInput, checker))
	{
		adfInput->close();
		delete adfInput;
		adfInput = new std::ifstream(dirname + "/" + inFileName);
	}
	adfInput->close();
	delete adfInput;

	std::cout << "input file written" << std::endl;

	std::string commands = "";
	commands += "cd " + dirname + "\n";
	commands += "dos2unix " + inFileName + "\n";
	commands += "nohup g16 <" + inFileName + " >" + outFileName + " &\n";
	std::cout << "printing input" << std::endl;


	ac = new char[commands.length() + 1];
	for (int i = 0; i < commands.length(); i++)
	{
		ac[i] = commands[i];
	}
	ac[commands.length()] = '\0';
	std::cout << "commands: " << ac << std::endl;

	std::cout << "writing command" << std::endl;
	std::cout << "commands: " << ac << std::endl;
	time_t start = time(0);
	std::system(ac);
	delete[] ac;

	std::ifstream* aOutFile = new std::ifstream(dirname + "/" + outFileName);


	std::cout << "checking for gaussian completion" << std::endl;


	bool failedToStart = false;
	while (!getline(*aOutFile, checker) && !failedToStart)
	{

		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);

		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				failedToStart = true;
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				failedToStart = true;
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;
			}

		}
	}
	//faiif (failedToStart) return DBL_MAX;
	std::cout << "Gaussian started" << std::endl;

	bool finishedOutput = false;
	bool finishedEnergy = false;
	std::vector<std::string> optStrings = {};
	std::string outLine;
	std::string energyString;



	std::string energyLine = "";
	while (!finishedOutput || !finishedEnergy)
	{
		if (finishedOutput && !(finishedEnergy))
		{

			commands = "pkill g16\n";
			ac = new char[commands.length() + 1];
			for (int i = 0; i < commands.length(); i++)
			{
				ac[i] = commands[i];
			}
			ac[commands.length()] = '\0';
			std::cout << "commands: " << ac << std::endl;
			std::cout << "writing command" << std::endl;
			std::system(ac);
			delete[] ac;
			return DBL_MAX;

		}
		if (timeLeftMinutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

			}

		}
		if (timeoutminutes != -1)
		{

			auto seconds_since_start = difftime(time(0), start);
			if (seconds_since_start > timeoutminutes * 60)
			{
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				time_t start = time(0);
				std::system(ac);
				delete[] ac;
				return DBL_MAX;
			}

		}

		optStrings = {};
		aOutFile->close();
		delete aOutFile;
		aOutFile = new std::ifstream(dirname + "/" + outFileName);

		std::string match = " Job cpu time:";
		std::string eMatch = " SCF Done:  E";
		while (getline(*aOutFile, outLine)) {
			if (outLine.size() > match.size())
			{
				bool matching = true;
				for (int m = 0; m < match.size(); m++)
				{
					if (outLine[m] != match[m]) matching = false;
				}
				if (matching) finishedOutput = true;//cannot be made false
				else {
					//std::cout << match << "\n vs \n" << outLine << std::endl;
				}
			}
			if (outLine.size() > eMatch.size())
			{
				bool matching = true;
				for (int m = 0; m < eMatch.size(); m++)
				{
					if (outLine[m] != eMatch[m]) matching = false;
				}
				if (matching) {
					finishedEnergy = true;//cannot be made false
					energyLine = outLine;
				}
			}
		}


	}
	std::cout << "finished output: " << finishedOutput << std::endl;
	double seconds_since_start = difftime(time(0), start);


	aOutFile->close();
	delete aOutFile;

	//find the start of the energy
	int eStart = 0;
	int eFin = 0;
	bool equalPassed = false;
	for (int i = 0; i < energyLine.size(); i++)
	{
		if (eFin) {}
		else if (eStart) {
			if (energyLine[i] == ' ')
			{
				eFin = i;
			}
		}
		else if (equalPassed) {
			if (energyLine[i] == ' ') {}
			else eStart = i;
		}
		else if (energyLine[i] == '=') equalPassed = true;

	}
	//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
	double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));

	std::cout << "Gaussian time: " << seconds_since_start << std::endl;
	if (!keepFiles)
	{

		std::string deleteFiles1 = "rm -r " + dirname + "\n";

		ac = new char[deleteFiles1.length() + 1];
		for (int i = 0; i < deleteFiles1.length(); i++)
		{
			ac[i] = deleteFiles1[i];
		}
		std::cout << "command: " << ac << std::endl;
		std::system(ac);
		delete[] ac;


	}
	return energyHartree;
}


std::vector<double> direction(structure& a, structure& b)
{
	//returns the direction between two structures
	std::vector<double> dir;
	for (int i = 0; i < a.elements.size(); i++)
	{
		for (int j = 0; j < a.set[i].size(); j++)
		{
			dir.push_back(b.set[i][j].x - a.set[i][j].x);
			dir.push_back(b.set[i][j].y - a.set[i][j].y);
			dir.push_back(b.set[i][j].z - a.set[i][j].z);
		}
	}

	return dir;
}

/*double energyCalculation(structure* s, std::string basis, std::string method, int charge, int state, char calculator, std::string taskName, int energyCall, bool keepFiles, std::string converge, double timeLeftMinutes, double timeoutMinutes, bool unrestricted, std::string customBasis, std::string AMSHOME) {
	double retValue = 0;
	switch (calculator)
	{
	case 'a':
		//adf/scm energy
		retValue = energyADF(*s, taskName, basis, method, charge, state, AMSHOME, energyCall, converge, timeoutMinutes, timeLeftMinutes, keepFiles, unrestricted);
		break;
	case 'b':
		//basin hopping energy
		retValue = basinHoppingEnergy(*s, energyCall++, taskName, 0, 0, 0);
		break;
	case 'g':
		//gaussain energy
		retValue = newEnergyGaussian(*s, taskName, basis, method, charge, state, energyCall, 300, 0, timeoutMinutes, timeLeftMinutes, keepFiles, customBasis);
		//retValue = energyGuassian(*s, energyCall++, taskName, basis, method, charge, state, keepFiles);
		break;
	default:
		retValue = energyGuassian(*s, energyCall++, taskName, basis, method, charge, state, keepFiles);
		break;
	}
	return retValue;
}

structure* optimizationCalculation(char calculator, structure* Xnew, int step, std::string taskName, int optimizationIterations, int maxSCF, double x, double y, double z, std::vector<double>& optimizationMovement, std::vector<structure*>& structures, std::string basis, std::string method, int charge, int state, bool keepFiles, std::string converge, double timeLeftMinutes, double timeoutMinutes, bool unrestricted, std::string customBasis, std::string AMSHOME )
{
	structure* optimizedStructure;
	switch (calculator)
	{
	case 'a':
	{
		optimizedStructure = optimizationADF(*Xnew, taskName, basis, method, charge, state, AMSHOME, step, maxSCF, converge, timeoutMinutes, timeLeftMinutes, keepFiles, unrestricted);
		break;
	}
	case 'b':
	{
		optimizedStructure = optimizationBasinHopping(*Xnew, step, taskName, optimizationIterations, x, y, z);
		break;
	}
	case 'g':
	{
		optimizedStructure = newOptimizationGaussian(*Xnew, taskName, basis, method, charge, state, step, maxSCF, optimizationIterations, timeoutMinutes, timeLeftMinutes, keepFiles, customBasis);
		//optimizedStructure = optimizationGaussian(*Xnew, step, taskName, basis, method, optimizationIterations, charge, state, maxSCF, keepFiles);
		break;
	}
	default:
	{
		optimizedStructure = optimizationGaussian(*Xnew, step, taskName, basis, method, optimizationIterations, charge, state, maxSCF, keepFiles);
		break;
	}
	}
	if (optimizedStructure == nullptr)
	{
		if (DEBUG) std::cout << "failed to optimize" << std::endl;
		//optimizedStructure = Xnew;//give it back Xnew with a terrible energy
		//optimizedStructure->energy = DBL_MAX;
	}
	else {
		optimizationMovement = direction(*Xnew, *optimizedStructure);
	}
	return optimizedStructure;
}
*/

/*
Calculators from premade dft files
*/
structure* optimization(std::string calculator, structure* x, std::string dftFileName,std::string pythonFile, double timeoutminutes, double timeLeftMinutes, std::string xyzfilename, int calculationNo)
{
	
	/*

	*/

	std::string dirname = "calc" + std::to_string(calculationNo);
	char* ac;
	std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";
	ac = new char[mkdir.length() + 1];
	for (int i = 0; i < mkdir.length(); i++)
	{
		ac[i] = mkdir[i];
	}
	std::cout << "command: " << ac << std::endl;
	std::system(ac);
	delete[] ac;
	writeToXyz({ x }, dirname + "/" + xyzfilename + ".xyz");

	

	std::string iFT = "";
	std::string oFT = "";

	if (calculator == "ADF") {
		iFT = ".run"; oFT = ".out";
	}
	else if (calculator == "Gaussian")
	{
		iFT = ".com"; oFT = ".log";
	}
	

	std::ifstream dftFile(dftFileName + iFT);

	std::ofstream newDFTFile(dirname + "/" + dftFileName + iFT);

	if (calculator == "Gaussian")
	{
		//gaussian does not allow xyz input.
		std::string gline = "";
		while (getline(dftFile, gline))
		{
			if (gline.substr(0, std::string("GEOM").size()) == "GEOM")
			{
				std::string addline = "";
				for (int e = 0; e < x->elements.size(); e++)
				{
					
					for (int a = 0; a < x->set[e].size(); a++) {
						newDFTFile << x->elements[e] << " " << std::to_string(x->set[e][a].x) << " " << std::to_string(x->set[e][a].y) << " " << std::to_string(x->set[e][a].z) << std::endl;

					}

				}
			}
			else {
				newDFTFile << gline << std::endl;
			}
		}
		dftFile.close();
		newDFTFile.close();
	}
	else {


	newDFTFile << dftFile.rdbuf();
	dftFile.close();
	newDFTFile.close();
	}

	if (calculator == "P")
	{
		std::string commands = "";
		//commands += "cd " + dirname + "\n";
		//commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "nohup python " + pythonFile + " O " + dirname + " &\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + "opt.xyz");
		std::cout << "checking for completion " << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timedOut = false;
		while (!getline(*aOutFile, checker) && !timedOut)
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + "opt.xyz");


			/*
			Timeout code
			*/
			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1)
			{
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}
				//else {
				//	std::cout << "time not passed sss: " << seconds_since_start << " !> " << timeLeftMinutes * 60 - 5 * 60 << std::endl;
				//}
			}
			if (timeoutminutes != -1)
			{
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}
			}
		}
		std::cout << "Optimization started writing" << std::endl;

		bool finished = false;
		bool failure = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;
		double energyHartree;
		int natoms = 0;
		for (int nn = 0; nn < x->composition.size(); nn++)
		{
			natoms += x->composition[nn];
		}
		while (!timedOut && (!finished))
		{
			
			if (timeLeftMinutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;

				}

			}
			if (timeoutminutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeoutminutes * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}

			}

			optStrings = {};
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + "opt.xyz");

			getline(*aOutFile, outLine);
			//std::cout << "OL:" << outLine << std::endl;
			if (outLine == "-1")
			{
				std::cout << "optimization failed." << std::endl;
				finished = true;
				failure = true;
			}
			else {
				int nlines = 0;
				while (getline(*aOutFile, outLine)) nlines++;
				if (nlines == natoms + 1)
				{
					finished = true;

				}
			}
			





		}
		
		double seconds_since_start = difftime(time(0), start);

		if (failure == true || finished == false)
		{
			return nullptr;
		}

		
		structure* st = structuresFromXYZ(dirname + "/" + "opt.xyz")[0];
		st->assignment = x->assignment;
		return st;
	}

	if (calculator == "W")
	{
		std::string commands = "";
		//commands += "cd " + dirname + "\n";
		//commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "Python " + pythonFile + " O " + dirname + "\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + "opt.xyz");
		std::cout << "checking for completion " << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timedOut = false;
		while (!getline(*aOutFile, checker) && !timedOut)
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + "opt.xyz");


			/*
			Timeout code
			*/
			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1)
			{
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}
				//else {
				//	std::cout << "time not passed sss: " << seconds_since_start << " !> " << timeLeftMinutes * 60 - 5 * 60 << std::endl;
				//}
			}
			if (timeoutminutes != -1)
			{
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}
			}
		}
		std::cout << "Optimization started writing" << std::endl;

		bool finished = false;
		bool failure = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;
		double energyHartree;
		int natoms = 0;
		for (int nn = 0; nn < x->composition.size(); nn++)
		{
			natoms += x->composition[nn];
		}
		while (!timedOut && (!finished))
		{

			if (timeLeftMinutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;

				}

			}
			if (timeoutminutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeoutminutes * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}

			}

			optStrings = {};
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + "opt.xyz");

			getline(*aOutFile, outLine);
			//std::cout << "OL:" << outLine << std::endl;
			if (outLine == "-1")
			{
				std::cout << "optimization failed." << std::endl;
				finished = true;
				failure = true;
			}
			else {
				int nlines = 0;
				while (getline(*aOutFile, outLine)) nlines++;
				if (nlines == natoms + 1)
				{
					finished = true;

				}
			}






		}

		double seconds_since_start = difftime(time(0), start);

		if (failure == true || finished == false)
		{
			return nullptr;
		}


		structure* st = structuresFromXYZ(dirname + "/" + "opt.xyz")[0];
		st->assignment = x->assignment;
		return st;
	}

	if (calculator == "ADF") {
		std::string commands = "";
		commands += "cd " + dirname + "\n";
		commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "nohup ./" + dftFileName + iFT + " > " + dftFileName + oFT + " &\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::cout << "checking for ADF completion" << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timeOut = false;
		while (!getline(*aOutFile, checker) && !timeOut)
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);


			/*
			Timeout code
			*/
			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1)
			{
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: Skipping structure optimization due to error with computer" << std::endl;
					timeOut = true;
					//return nullptr;
				}
			}
			if (seconds_since_start > timeoutminutes * 60)
			{
				std::cout << "WARNING: ADF did not start" << std::endl;
				timeOut = true;
				//return nullptr;
			}
			
		}
		std::cout << "ADF started" << std::endl;

		//the file is writing complete?
		bool finishedOutput_energy = false;
		bool finishedOutput_geometry = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;
		
		bool currentGeometry = false;//checks if there is at least some geometry optimized from input
		bool currentEnergy = false;
		std::string currentEnergyString = "";

		int natoms = 0;
		for (int i = 0; i < x->composition.size(); i++) natoms += x->composition[i];
		//std::cout << "RML natoms: " << natoms << std::endl;
		bool failed = false;
		while (!timeOut && (!finishedOutput_energy || !finishedOutput_geometry) && !failed)
		{
			if (timeLeftMinutes != -1)
			{
				
				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					//kill the task so we dont waste cpu resource
					commands = "pkill ams\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timeOut = true;
					//return nullptr;
					std::cout << "WARNING: Skipping structure optimization due to error with calculator" << std::endl;

				}

			}
			if (timeoutminutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeoutminutes * 60)
				{
					commands = "pkill ams\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timeOut = true;
					//return nullptr;
				}

			}

			optStrings = {};
			//reopen file to get new updates
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + ".out");

			//look for the energy line

			while (getline(*aOutFile, outLine) && !failed) {

				//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

				std::string match = "Energy (hartree)";
				std::string match2 = "Optimized geometry:";
				std::string match3 = "Geometry optimization failed";
				std::string match4 = "BAD CORE INTEGRAL";
				std::string match5 = "Process received SIGTERM";

				std::string someGeometry = "current energy";
				std::string someOptimization = "--------";
				std::string geometry = "Geometry";
				if (outLine.size() >= someGeometry.size())
				{
					bool matching = true;
					for (int m = 0; m < someGeometry.size(); m++)
					{
						if (outLine[m] != someGeometry[m]) matching = false;
					}
					if (matching) {
						currentEnergy = true;//cannot be made false
						currentEnergyString = outLine;
					}
				}
				bool intermediateGeomCheck = false;
				if (outLine.size() >= someOptimization.size())
				{
					bool matching = true;
					for (int m = 0; m < someOptimization.size(); m++)
					{
						if (outLine[m] != someOptimization[m]) matching = false;
					}
					if (matching) {
						intermediateGeomCheck = true;//cannot be made false
					}
				}
				if (intermediateGeomCheck)
				{
					bool getOpt = false;
					if (!getline(*aOutFile, outLine))
					{


						if (outLine.size() >= geometry.size())
						{
							bool matching = true;
							for (int m = 0; m < geometry.size(); m++)
							{
								if (outLine[m] != geometry[m]) matching = false;
							}
							if (matching) {
								getOpt = true;//cannot be made false
							}
						}
					}
					if (getOpt)
					{
						std::vector<std::string> tempOptStrings = {};

						bool stop = false;
						for (int i = 0; i < 4; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
						for (int e = 0; e < x->elements.size(); e++)  for (int a = 0; a < x->set[e].size(); a++) {
							if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
							tempOptStrings.push_back(outLine);
						}
						if (tempOptStrings.size() == natoms) {
							optStrings = tempOptStrings;
							currentGeometry = true;
						}
					}
				}



				if (outLine.size() >= match.size())
				{
					bool matching = true;
					for (int m = 0; m < match.size(); m++)
					{
						if (outLine[m] != match[m]) matching = false;
					}
					if (matching) {
						finishedOutput_energy = true;//cannot be made false
						energyString = outLine;
					}
				}

				bool getopt = false;

				if (outLine.size() >= match2.size())
				{
					bool matching = true;
					for (int m = 0; m < match2.size(); m++)
					{
						if (outLine[m] != match2[m]) matching = false;
					}
					if (matching) {
						getopt = true;//cannot be made false
						finishedOutput_geometry = true;//geometry needs its own bool now because there are cases where energy calculates but geometry does not and crash line geometry optimization failed hasnt been printed yet
					}
				}
				if (getopt)
				{
					std::vector<std::string> tempOptStrings = {};
					bool stop = false;
					for (int i = 0; i < 7; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
					for (int e = 0; e < x->elements.size(); e++)  for (int a = 0; a < x->set[e].size(); a++) {
						if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
						tempOptStrings.push_back(outLine);
					}
					if (tempOptStrings.size() == natoms) optStrings = tempOptStrings;
				}
				//when the AMS calculation is successful, energy is always printed after the geometry, and therefore its fine to get the geometry here.

				failed = false;
				if (outLine.size() >= match3.size())
				{
					bool matching = true;
					for (int m = 0; m < match3.size(); m++)
					{
						if (outLine[m] != match3[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				if (outLine.size() >= match4.size())
				{
					bool matching = true;
					for (int m = 0; m < match4.size(); m++)
					{
						if (outLine[m] != match4[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				if (outLine.size() >= match5.size())
				{
					bool matching = true;
					for (int m = 0; m < match5.size(); m++)
					{
						if (outLine[m] != match5[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				/*if (failed)
				{
					return nullptr;
				}*/
			}


		}
		std::cout << "finished output: " << finishedOutput_energy << std::endl;

		if (timeOut && (!currentEnergy || !currentGeometry))
		{
			failed = true;
			//return nullptr;
		}
		if (failed && optStrings.size() != natoms)
		{
			//got no structure; complete failure
			//std::cout << "ADF optimization failure: opt strings size: " << optStrings.size() << std::endl;
			return nullptr;
		}

		double seconds_since_start = difftime(time(0), start);


		std::cout << "ADF time: " << seconds_since_start << std::endl;
		//std::cout << "RML opt strings size" << optStrings.size();
		structure* st = nullptr;



		std::vector<std::string> elements = {};
		std::vector<std::vector<atom>> set = {};

		for (int i = 0; i < optStrings.size(); i++)
		{
			std::cout << optStrings[i] << std::endl;

			std::string cline = optStrings[i];

			std::string element = "";
			if (cline[12] == ' ') element = cline.substr(13, 1);//it goes the other way ....
			else element = cline.substr(12, 2);
			//x
			double x;
			double y;
			double z;
			bool sign = false;
			if (cline[18] == '-') sign = true;
			if (sign) x = std::stod(cline.substr(18, 11));
			else x = std::stod(cline.substr(19, 10));
			if (cline[33] == '-') sign = true;
			else sign = false;
			if (sign) y = std::stod(cline.substr(33, 11));
			else y = std::stod(cline.substr(34, 10));
			if (cline[48] == '-') sign = true;
			else sign = false;
			if (sign) z = std::stod(cline.substr(48, 11));
			else z = std::stod(cline.substr(49, 10));



			//figure out which set this belongs to.
			int index = -1;
			for (int i = 0; i < elements.size(); i++)
			{
				if (elements[i] == element)
				{
					index = i;
				}
			}
			//if (buggo) std::cout << "making atom for elemental index " << index << " with coords " << x << " " << y << " " << z << std::endl;
			if (index == -1)
			{
				//if (buggo) std::cout << "not yet indexed" << std::endl;
				//atom type not yet indexed
				elements.push_back(element);
				std::vector < atom> alist = { atom(x,y,z) };
				set.push_back(alist);
			}
			else {
				//if (buggo) std::cout << "pushing back" << std::endl;
				set[index].push_back(atom(x, y, z));
			}
		}
		//if (buggo) std::cout << "this structure is not the last one starting over with next" << std::endl;

		if (st != nullptr) delete st;
		st = new structure(set, elements);


		std::cout << "found energy" << std::endl;

		std::string matchString;
		if (!timeOut) matchString = "Energy (hartree)";
		else {
			matchString = "current energy";
			energyString = currentEnergyString;
		}
		std::string num = "";
		for (int i = matchString.length(); i < energyString.length(); i++)
		{
			if (energyString[i] != ' ') num += energyString[i];
		}
		double energy = std::stod(num);
		st->energy = energy;
		// Close the file
		//if there was no convergence then nullptr will be returned
		//if (buggo) std::cout << "closing" << std::endl;
		aOutFile->close();
		delete aOutFile;


		//remove files


		return st;
		st->assignment = x->assignment;
		//return energyHartree;	
	}
	else if (calculator == "Gaussian")
	{
		std::string commands = "";
		commands += "cd " + dirname + "\n";
		commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "nohup g16 <" + dftFileName + iFT + " >" + dftFileName + oFT + " &\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::cout << "checking for Gaussian completion" << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timedOut = false;
		while (!getline(*aOutFile, checker) && !timedOut)
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);


			/*
			Timeout code
			*/
			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1)
			{
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: Skipping structure optimization due to error with computer" << std::endl;
					timedOut = true;
					//return nullptr;
				}
			}
			if (timeoutminutes != -1)
			{
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: Gaussian did not start" << std::endl;
					timedOut = true;
				}
			}
		}
		std::cout << "Gaussian started" << std::endl;

		bool finishedOutput = false;
		bool finishedEnergy = false;
		bool finishedStructure = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;


		while (!timedOut && (!finishedOutput || !finishedEnergy || !finishedStructure))
		{
			if (finishedOutput && !(finishedEnergy && finishedStructure))
			{
				std::cout << "finished output but did not finish finished output: " << finishedOutput << " , finished energy: " << finishedEnergy << " finished structure:" << finishedStructure << std::endl;
				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				std::system(ac);
				delete[] ac;
				//timedOut = true;
				return nullptr;

			}
			if (timeLeftMinutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "pkill g16\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;

				}

			}
			if (timeoutminutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeoutminutes * 60)
				{
					commands = "pkill g16\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					//return nullptr;
				}

			}

			optStrings = {};
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

			std::string match = " Job cpu time:";
			std::string matchString = "                          Input orientation:";
			std::string eMatch = " SCF Done:  E";
			std::string fail = " Error termination";
			while (getline(*aOutFile, outLine)) {
				if (outLine.size() > match.size())
				{
					bool matching = true;
					for (int m = 0; m < match.size(); m++)
					{
						if (outLine[m] != match[m]) matching = false;
					}
					if (matching) finishedOutput = true;//cannot be made false
				}
				if (outLine.size() > matchString.size())
				{
					bool matching = true;
					for (int m = 0; m < matchString.size(); m++)
					{
						if (outLine[m] != matchString[m]) matching = false;
					}
					if (matching) finishedStructure = true;//cannot be made false
				}
				if (outLine.size() > eMatch.size())
				{
					bool matching = true;
					for (int m = 0; m < eMatch.size(); m++)
					{
						if (outLine[m] != eMatch[m]) matching = false;
					}
					if (matching) finishedEnergy = true;//cannot be made false
				}
				if (outLine.size() > fail.size())
				{
					bool matching = true;
					for (int m = 0; m < fail.size(); m++)
					{
						if (outLine[m] != fail[m]) matching = false;
					}
					if (matching) {
						std::cout << "Optimization failed due to a problem in gaussian. error code: " << outLine << std::endl;
						//std::cout << "finished output but did not finish finished output: " << finishedOutput << " , finished energy: " << finishedEnergy << " finished structure:" << finishedStructure << std::endl;

						if (!timedOut && ( !finishedEnergy || !finishedStructure))return nullptr;
						else std::cout << "continuing as calculation was still successful" << std::endl;
					}
				}
			}





		}
		std::cout << "finished output: " << finishedOutput << std::endl;
		double seconds_since_start = difftime(time(0), start);

		if (!finishedEnergy || !finishedOutput)
		{
			return nullptr;//otherwise we at least can get an energy and geometry
		}

		aOutFile->close();
		delete aOutFile;
		std::ifstream* gOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::string cline;
		structure* st = nullptr;

		std::string matchString = "Input orientation:";
		std::string eMatch = " SCF Done:  E";
		std::string energyLine = "";

		std::cout << "reading file" << std::endl;

		while (getline(*gOutFile, cline)) {
			if (cline.length() > matchString.length())
			{

				bool matching = false;
				int matchIt = 0;
				for (int l = 0; l < cline.length() - matchString.length(); l++)
				{
					if (cline[l] == matchString[matchIt])
					{
						matchIt++;
						matching = true;
					}
					else if (matchIt < 18) {
						matchIt = 0;
						matching = false;
					}
					else {
						//the strings are already matched let the match go through
					}
				}
				//if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
				if (matching && matchIt >= matchString.length())
				{
					std::vector<std::vector<atom>> set;
					std::vector<std::string> elements;
					getline(*gOutFile, cline);
					getline(*gOutFile, cline);
					getline(*gOutFile, cline);
					getline(*gOutFile, cline);
					char space = ' ';
					while (getline(*gOutFile, cline) && cline[1] != '-')
					{
						std::string cnumber = "";
						int num = 0;
						std::string element;
						double x = 0;
						double y = 0;
						double z = 0;
						//std::cout << "cline: " << cline << std::endl;
						for (int i = 0; i < cline.length(); i++)
						{
							if (cline[i] == space)
							{
								//std::cout << "space detected, cnumber: " << cnumber << std::endl;
								if (cnumber.length() > 0)
								{
									num++;
									switch (num)
									{
									case 1:
										//center number, not sure what this is for. do nothing
										break;
									case 2:
										//atomic number
										element = atomicNumber(std::stoi(cnumber));
										//std::cout << "element: " << element << std::endl;
										break;
									case 3:
										//atomic type, not sure what this is. do nothing.
										break;
									case 4:
										//x coordinate
										x = std::stof(cnumber);
										//std::cout << "x: " << x << std::endl;
										break;
									case 5:
										//y
										y = std::stof(cnumber);
										//std::cout << "y: " << y << std::endl;
										break;
									case 6:
										//z
										z = std::stof(cnumber);
										//std::cout << "z: " << z << std::endl;
										break;
									default:
										break;

									}
									cnumber = "";
								}
							}
							else {
								cnumber += cline[i];
							}
						}
						//make sure z is added
						z = std::stof(cnumber);
						//std::cout << "z: " << z << std::endl;
						int index = -1;
						for (int i = 0; i < elements.size(); i++)
						{
							if (elements[i] == element)
							{
								index = i;
							}
						}
						if (index == -1)
						{
							elements.push_back(element);
							std::vector < atom> alist = { atom(x,y,z) };
							set.push_back(alist);
						}
						else {
							set[index].push_back(atom(x, y, z));
						}
					}
					if (st != nullptr) delete st;
					st = new structure(set, elements);


				}
			}
			if (cline.size() > eMatch.size())
			{
				bool matching = true;
				for (int m = 0; m < eMatch.size(); m++)
				{
					if (cline[m] != eMatch[m]) matching = false;
				}
				if (matching) energyLine = cline;
			}
			

		}
		gOutFile->close();
		delete gOutFile;

		//find the start of the energy
		int eStart = 0;
		int eFin = 0;
		bool equalPassed = false;
		for (int i = 0; i < energyLine.size(); i++)
		{
			if (eFin) {}
			else if (eStart) {
				if (energyLine[i] == ' ')
				{
					eFin = i;
				}
			}
			else if (equalPassed) {
				if (energyLine[i] == ' ') {}
				else eStart = i;
			}
			else if (energyLine[i] == '=') equalPassed = true;

		}
		//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
		double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));
		st->energy = energyHartree;

		std::cout << "Gaussian time: " << seconds_since_start << std::endl;
		st->assignment = x->assignment;
		return st;
	}
	return nullptr;
}
double energyDFT(std::string calculator, structure* x, std::string dftFileName,std::string pythonFile, double timeoutminutes, double timeLeftMinutes, std::string xyzfilename, int calculationNo)
{
	/*

	*/

	std::string dirname = "calc" + std::to_string(calculationNo);
	char* ac;
	std::string mkdir = "mkdir " + dirname + "\ncd " + dirname + "\n";
	ac = new char[mkdir.length() + 1];
	for (int i = 0; i < mkdir.length(); i++)
	{
		ac[i] = mkdir[i];
	}
	std::cout << "command: " << ac << std::endl;
	std::system(ac);
	delete[] ac;
	writeToXyz({ x }, dirname + "/" + xyzfilename + ".xyz");

	std::string iFT = "";
	std::string oFT = "";

	if (calculator == "ADF") {
		iFT = ".run"; oFT = ".out";
	}
	else if (calculator == "Gaussian")
	{
		iFT = ".com"; oFT = ".log";
	}
	//std::cout << "calcilator:"<<calculator << "|"<<iFT<<std::endl;
	std::ifstream dftFile(dftFileName + iFT);

	//std::cout << dirname + "/" + dftFileName + iFT;
	std::ofstream newDFTFile(dirname + "/" + dftFileName + iFT);

	if (calculator == "Gaussian")
	{
		//gaussian does not allow xyz input.
		std::string gline = "";
		while (getline(dftFile, gline))
		{
			if (gline.substr(0, std::string("GEOM").size()) == "GEOM")
			{
				std::string addline = "";
				for (int e = 0; e < x->elements.size(); e++)
				{

					for (int a = 0; a < x->set[e].size(); a++) {
						newDFTFile << x->elements[e] << " " << std::to_string(x->set[e][a].x) << " " << std::to_string(x->set[e][a].y) << " " << std::to_string(x->set[e][a].z) << std::endl;
					}

				}
			}
			else {
				newDFTFile << gline << std::endl;
			}
		}
		dftFile.close();
		newDFTFile.close();
	}
	else {


		newDFTFile << dftFile.rdbuf();
		dftFile.close();
		newDFTFile.close();
	}

	if (calculator == "P")
	{
		std::string commands = "";
		//commands += "cd " + dirname + "\n";
		//commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "nohup python " + pythonFile + " E " + dirname + " &\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + "e.txt");
		std::cout << "checking for completion" << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timedOut = false;
		
		bool finished = false;
		bool failure = false;
		double energyHartree = DBL_MAX;
		std::string outLine;
		while (!timedOut && (!finished))
		{
			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1)
			{

				
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					return DBL_MAX;

				}

			}
			if (timeoutminutes != -1)
			{
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					return DBL_MAX;
				}
			}
			
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + "e.txt");

			if (getline(*aOutFile, outLine)) {
				if (outLine == "F")
				{
					std::cout << "calculation failed." << std::endl;
					finished = true;
					failure = true;
				}
				else {
					finished = true;
					failure = false;
					std::cout << "energy out " << outLine << std::endl;
					energyHartree = std::stod(outLine);
				}
			}





		}

		double seconds_since_start = difftime(time(0), start);

		if (failure == true || finished == false)
		{
			return DBL_MAX;
		}
		else {
			return energyHartree;
		}
	}

	if (calculator == "W")
	{
		std::string commands = "";
		//commands += "cd " + dirname + "\n";
		//commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "Python " + pythonFile + " E " + dirname + "\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + "e.txt");
		std::cout << "checking for completion" << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timedOut = false;

		bool finished = false;
		bool failure = false;
		double energyHartree = DBL_MAX;
		std::string outLine;
		while (!timedOut && (!finished))
		{
			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1)
			{


				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					return DBL_MAX;

				}

			}
			if (timeoutminutes != -1)
			{
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					commands = "python " + pythonFile + " cancel\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					timedOut = true;
					return DBL_MAX;
				}
			}

			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + "e.txt");

			if (getline(*aOutFile, outLine)) {
				if (outLine == "F")
				{
					std::cout << "calculation failed." << std::endl;
					finished = true;
					failure = true;
				}
				else {
					finished = true;
					failure = false;
					std::cout << "energy out " << outLine << std::endl;
					energyHartree = std::stod(outLine);
				}
			}





		}

		double seconds_since_start = difftime(time(0), start);

		if (failure == true || finished == false)
		{
			return DBL_MAX;
		}
		else {
			return energyHartree;
		}
	}

	if (calculator == "ADF") {
		std::string commands = "";
		commands += "cd " + dirname + "\n";
		commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "nohup ./" + dftFileName + iFT + " > " + dftFileName + oFT + " &\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;



		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::cout << "checking for ADF completion" << std::endl;

		std::string checker;
		std::cout << "checking for adf completion" << std::endl;




		while (!getline(*aOutFile, checker))
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1) {
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: Skipping structure energy calculation due to error with computer" << std::endl;
					return DBL_MAX;
					//return nullptr;
				}
			}
			if (timeoutminutes != -1) {
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: ADF never started" << std::endl;
					return DBL_MAX;
					//return nullptr;
				}
			}
			if (seconds_since_start > 60)
			{
				std::cout << "WARNING: ADF never started" << std::endl;
				return DBL_MAX;
				//return nullptr;
			}
		}
		std::cout << "ADF started" << std::endl;

		//the file is writing complete?
		bool finishedOutput = false;
		std::string outLine;
		std::string energyString;
		while (!finishedOutput)
		{
			if (timeLeftMinutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					//kill the task so we dont waste cpu resource
					commands = "pkill ams\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					return DBL_MAX;

				}

			}
			if (timeoutminutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeoutminutes * 60)
				{
					commands = "pkill ams\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					std::cout << "WARNING: Returning DBL_MAX for structure energy calculation due to error with energy calculator" << std::endl;
					return DBL_MAX;
				}

			}


			//reopen file to get new updates
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

			//look for the energy line

			while (getline(*aOutFile, outLine)) {
				//come back here plz
				//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

				std::string match = "Energy (hartree)";
				std::string match4 = "BAD CORE INTEGRAL";
				std::string match5 = "Process received SIGTERM";
				if (outLine.size() >= match.size())
				{
					bool matching = true;
					for (int m = 0; m < match.size(); m++)
					{
						if (outLine[m] != match[m]) matching = false;
					}
					if (matching) {
						finishedOutput = true;//cannot be made false
						energyString = outLine;
					}
				}
				bool failed = false;
				if (outLine.size() >= match4.size())
				{
					bool matching = true;
					for (int m = 0; m < match4.size(); m++)
					{
						if (outLine[m] != match4[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				if (outLine.size() >= match5.size())
				{
					bool matching = true;
					for (int m = 0; m < match5.size(); m++)
					{
						if (outLine[m] != match5[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				if (failed)
				{
					return DBL_MAX;
				}
			}



		}
		std::cout << "finished output: " << finishedOutput << std::endl;


		double seconds_since_start = difftime(time(0), start);
		std::cout << "ADF time: " << seconds_since_start << std::endl;


		std::cout << "found energy" << std::endl;
		std::string matchString = "Energy (hartree)";
		std::string num = "";
		for (int i = matchString.length(); i < energyString.length(); i++)
		{
			if (energyString[i] != ' ') num += energyString[i];
		}
		double energy = std::stod(num);
		// Close the file
		//if there was no convergence then nullptr will be returned
		//if (buggo) std::cout << "closing" << std::endl;
		aOutFile->close();
		delete aOutFile;


		//remove files
		return energy;
	}
	else if (calculator == "Gaussian")
	{
		

		std::string commands = "";
		commands += "cd " + dirname + "\n";
		commands += "dos2unix " + dftFileName + iFT + "\n" + "chmod u+x " + dftFileName + iFT + "\n";
		commands += "nohup g16 <" + dftFileName + iFT + " >" + dftFileName + oFT + " &\n";
		std::cout << "printing input" << std::endl;


		ac = new char[commands.length() + 1];
		for (int i = 0; i < commands.length(); i++)
		{
			ac[i] = commands[i];
		}
		ac[commands.length()] = '\0';
		std::cout << "commands: " << ac << std::endl;

		std::cout << "writing command" << std::endl;
		std::cout << "commands: " << ac << std::endl;
		time_t start = time(0);
		std::system(ac);
		delete[] ac;



		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::cout << "checking for Gaussian completion" << std::endl;

		std::string checker;
		std::cout << "checking for Gaussian completion" << std::endl;




		while (!getline(*aOutFile, checker))
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

			auto seconds_since_start = difftime(time(0), start);
			if (timeLeftMinutes != -1) {
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: Skipping structure energy calculation due to error with computer" << std::endl;
					return DBL_MAX;
					//return nullptr;
				}
			}
			if (timeoutminutes != -1) {
				if (seconds_since_start > timeoutminutes * 60 - 5 * 60)
				{
					std::cout << "WARNING: Gaussian never started" << std::endl;
					return DBL_MAX;
					//return nullptr;
				}
			}
		}
		std::cout << "Gaussian started" << std::endl;

		
		bool finishedOutput = false;
		bool finishedEnergy = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;



		std::string energyLine = "";
		while (!finishedOutput || !finishedEnergy)
		{
			if (finishedOutput && !(finishedEnergy))
			{

				commands = "pkill g16\n";
				ac = new char[commands.length() + 1];
				for (int i = 0; i < commands.length(); i++)
				{
					ac[i] = commands[i];
				}
				ac[commands.length()] = '\0';
				std::cout << "commands: " << ac << std::endl;
				std::cout << "writing command" << std::endl;
				std::system(ac);
				delete[] ac;
				return DBL_MAX;

			}
			if (timeLeftMinutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeLeftMinutes * 60 - 5 * 60)
				{
					commands = "pkill g16\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					return DBL_MAX;

				}

			}
			if (timeoutminutes != -1)
			{

				auto seconds_since_start = difftime(time(0), start);
				if (seconds_since_start > timeoutminutes * 60)
				{
					commands = "pkill g16\n";
					ac = new char[commands.length() + 1];
					for (int i = 0; i < commands.length(); i++)
					{
						ac[i] = commands[i];
					}
					ac[commands.length()] = '\0';
					std::cout << "commands: " << ac << std::endl;
					std::cout << "writing command" << std::endl;
					time_t start = time(0);
					std::system(ac);
					delete[] ac;
					return DBL_MAX;
				}

			}

			optStrings = {};
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

			std::string match = " Job cpu time:";
			std::string eMatch = " SCF Done:  E";
			std::string fail = " Error termination";
			while (getline(*aOutFile, outLine)) {
				if (outLine.size() > match.size())
				{
					bool matching = true;
					for (int m = 0; m < match.size(); m++)
					{
						if (outLine[m] != match[m]) matching = false;
					}
					if (matching) finishedOutput = true;//cannot be made false
					else {
						//std::cout << match << "\n vs \n" << outLine << std::endl;
					}
				}
				if (outLine.size() > eMatch.size())
				{
					bool matching = true;
					for (int m = 0; m < eMatch.size(); m++)
					{
						if (outLine[m] != eMatch[m]) matching = false;
					}
					if (matching) {
						finishedEnergy = true;//cannot be made false
						energyLine = outLine;
					}
				}
				if (outLine.size() > fail.size())
				{
					bool matching = true;
					for (int m = 0; m < fail.size(); m++)
					{
						if (outLine[m] != fail[m]) matching = false;
					}
					if (matching) {
						std::cout << "Energy failed due to a problem in gaussian" << std::endl;
						return DBL_MAX;
					}
				}
			}


		}
		std::cout << "finished output: " << finishedOutput << std::endl;
		double seconds_since_start = difftime(time(0), start);


		aOutFile->close();
		delete aOutFile;

		//find the start of the energy
		int eStart = 0;
		int eFin = 0;
		bool equalPassed = false;
		for (int i = 0; i < energyLine.size(); i++)
		{
			if (eFin) {}
			else if (eStart) {
				if (energyLine[i] == ' ')
				{
					eFin = i;
				}
			}
			else if (equalPassed) {
				if (energyLine[i] == ' ') {}
				else eStart = i;
			}
			else if (energyLine[i] == '=') equalPassed = true;

		}
		//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
		double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));

		std::cout << "Gaussian time: " << seconds_since_start << std::endl;
		return energyHartree;
	}
	return DBL_MAX;
}

structure* readOptFile(std::string calculator, std::string dftFileName, int calculationNo,std::vector<int> composition, std::vector<std::string> elements)
{
	/*

	*/

	std::string dirname = "calc" + std::to_string(calculationNo);

	std::string oFT = "";

	if (calculator == "ADF") {
		oFT = ".out";
	}
	else if (calculator == "Gaussian")
	{
		oFT = ".log";
	}
	
	if (calculator == "ADF") {
		
		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::cout << "checking for ADF completion" << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		while (!getline(*aOutFile, checker))
		{
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

		}
		std::cout << "ADF started" << std::endl;

		//the file is writing complete?
		bool finishedOutput_energy = false;
		bool finishedOutput_geometry = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;


		bool currentGeometry = false;//checks if there is at least some geometry optimized from input
		bool currentEnergy = false;
		std::string currentEnergyString = "";

		int natoms = 0;
		for (int i = 0; i < composition.size(); i++) natoms += composition[i];
		//std::cout << "RML natoms: " << natoms << std::endl;
		bool failed = false;

			optStrings = {};
			//reopen file to get new updates
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + ".out");

			//look for the energy line

			while (getline(*aOutFile, outLine) && !failed) {

				//if getline is included in the while statement conditional then energy line will advance after the loop closes and point to the wrong line

				std::string match = "Energy (hartree)";
				std::string match2 = "Optimized geometry:";
				std::string match3 = "Geometry optimization failed";
				std::string match4 = "BAD CORE INTEGRAL";
				std::string match5 = "Process received SIGTERM";

				std::string someGeometry = "current energy";
				std::string someOptimization = "--------";
				std::string geometry = "Geometry";
				if (outLine.size() >= someGeometry.size())
				{
					bool matching = true;
					for (int m = 0; m < someGeometry.size(); m++)
					{
						if (outLine[m] != someGeometry[m]) matching = false;
					}
					if (matching) {
						currentEnergy = true;//cannot be made false
						currentEnergyString = outLine;
					}
				}
				bool intermediateGeomCheck = false;
				if (outLine.size() >= someOptimization.size())
				{
					bool matching = true;
					for (int m = 0; m < someOptimization.size(); m++)
					{
						if (outLine[m] != someOptimization[m]) matching = false;
					}
					if (matching) {
						intermediateGeomCheck = true;//cannot be made false
					}
				}
				if (intermediateGeomCheck)
				{
					bool getOpt = false;
					if (!getline(*aOutFile, outLine))
					{


						if (outLine.size() >= geometry.size())
						{
							bool matching = true;
							for (int m = 0; m < geometry.size(); m++)
							{
								if (outLine[m] != geometry[m]) matching = false;
							}
							if (matching) {
								getOpt = true;//cannot be made false
							}
						}
					}
					if (getOpt)
					{
						std::vector<std::string> tempOptStrings = {};

						bool stop = false;
						for (int i = 0; i < 4; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
						for (int e = 0; e < elements.size(); e++)  for (int a = 0; a < composition[e]; a++) {
							if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
							tempOptStrings.push_back(outLine);
						}
						if (tempOptStrings.size() == natoms) {
							optStrings = tempOptStrings;
							currentGeometry = true;
						}
					}
				}



				if (outLine.size() >= match.size())
				{
					bool matching = true;
					for (int m = 0; m < match.size(); m++)
					{
						if (outLine[m] != match[m]) matching = false;
					}
					if (matching) {
						finishedOutput_energy = true;//cannot be made false
						energyString = outLine;
					}
				}

				bool getopt = false;

				if (outLine.size() >= match2.size())
				{
					bool matching = true;
					for (int m = 0; m < match2.size(); m++)
					{
						if (outLine[m] != match2[m]) matching = false;
					}
					if (matching) {
						getopt = true;//cannot be made false
						finishedOutput_geometry = true;//geometry needs its own bool now because there are cases where energy calculates but geometry does not and crash line geometry optimization failed hasnt been printed yet
					}
				}
				if (getopt)
				{
					std::vector<std::string> tempOptStrings = {};
					bool stop = false;
					for (int i = 0; i < 7; i++) if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
					for (int e = 0; e < elements.size(); e++)  for (int a = 0; a < composition[e]; a++) {
						if (!stop) if (!getline(*aOutFile, outLine)) stop = true;
						tempOptStrings.push_back(outLine);
					}
					if (tempOptStrings.size() == natoms) optStrings = tempOptStrings;
				}
				//when the AMS calculation is successful, energy is always printed after the geometry, and therefore its fine to get the geometry here.

				failed = false;
				if (outLine.size() >= match3.size())
				{
					bool matching = true;
					for (int m = 0; m < match3.size(); m++)
					{
						if (outLine[m] != match3[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				if (outLine.size() >= match4.size())
				{
					bool matching = true;
					for (int m = 0; m < match4.size(); m++)
					{
						if (outLine[m] != match4[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				if (outLine.size() >= match5.size())
				{
					bool matching = true;
					for (int m = 0; m < match5.size(); m++)
					{
						if (outLine[m] != match5[m]) matching = false;
					}
					if (matching) failed = true;//cannot be made false
				}
				/*if (failed)
				{
					return nullptr;
				}*/
			}
	
		std::cout << "finished output: " << finishedOutput_energy << std::endl;

		if ( (!currentEnergy || !currentGeometry))
		{
			failed = true;
			//return nullptr;
		}
		if (failed && optStrings.size() != natoms)
		{
			//got no structure; complete failure
			//std::cout << "ADF optimization failure: opt strings size: " << optStrings.size() << std::endl;
			return nullptr;
		}

		structure* st = nullptr;



		std::vector<std::string> elements = {};
		std::vector<std::vector<atom>> set = {};

		for (int i = 0; i < optStrings.size(); i++)
		{
			std::cout << optStrings[i] << std::endl;

			std::string cline = optStrings[i];

			std::string element = "";
			if (cline[12] == ' ') element = cline.substr(13, 1);//it goes the other way ....
			else element = cline.substr(12, 2);
			//x
			double x;
			double y;
			double z;
			bool sign = false;
			if (cline[18] == '-') sign = true;
			if (sign) x = std::stod(cline.substr(18, 11));
			else x = std::stod(cline.substr(19, 10));
			if (cline[33] == '-') sign = true;
			else sign = false;
			if (sign) y = std::stod(cline.substr(33, 11));
			else y = std::stod(cline.substr(34, 10));
			if (cline[48] == '-') sign = true;
			else sign = false;
			if (sign) z = std::stod(cline.substr(48, 11));
			else z = std::stod(cline.substr(49, 10));



			//figure out which set this belongs to.
			int index = -1;
			for (int i = 0; i < elements.size(); i++)
			{
				if (elements[i] == element)
				{
					index = i;
				}
			}
			//if (buggo) std::cout << "making atom for elemental index " << index << " with coords " << x << " " << y << " " << z << std::endl;
			if (index == -1)
			{
				//if (buggo) std::cout << "not yet indexed" << std::endl;
				//atom type not yet indexed
				elements.push_back(element);
				std::vector < atom> alist = { atom(x,y,z) };
				set.push_back(alist);
			}
			else {
				//if (buggo) std::cout << "pushing back" << std::endl;
				set[index].push_back(atom(x, y, z));
			}
		}
		//if (buggo) std::cout << "this structure is not the last one starting over with next" << std::endl;

		if (st != nullptr) delete st;
		st = new structure(set, elements);


		std::cout << "found energy" << std::endl;

		std::string matchString = "Energy (hartree)";

		std::string num = "";
		for (int i = matchString.length(); i < energyString.length(); i++)
		{
			if (energyString[i] != ' ') num += energyString[i];
		}
		double energy = std::stod(num);
		st->energy = energy;
		// Close the file
		//if there was no convergence then nullptr will be returned
		//if (buggo) std::cout << "closing" << std::endl;
		aOutFile->close();
		delete aOutFile;


		//remove files


		return st;
		//return energyHartree;	
	}
	else if (calculator == "Gaussian")
	{
		

		std::ifstream* aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::cout << "checking for Gaussian completion" << std::endl;

		std::string checker;

		/*
		TODO THIS WHILE LOOP AND ALL OTHERS LIKE IT NEED TO BE DEBUGGED
		*/
		bool timedOut = false;
		
		
		bool finishedOutput = false;
		bool finishedEnergy = false;
		bool finishedStructure = false;
		std::vector<std::string> optStrings = {};
		std::string outLine;
		std::string energyString;


			if (finishedOutput && !(finishedEnergy && finishedStructure))
			{
				return nullptr;
			}
			optStrings = {};
			aOutFile->close();
			delete aOutFile;
			aOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);

			std::string match = " Job cpu time:";
			std::string matchString = "                          Input orientation:";
			std::string eMatch = " SCF Done:  E";
			std::string fail = " Error termination";
			while (getline(*aOutFile, outLine)) {
				if (outLine.size() > match.size())
				{
					bool matching = true;
					for (int m = 0; m < match.size(); m++)
					{
						if (outLine[m] != match[m]) matching = false;
					}
					if (matching) finishedOutput = true;//cannot be made false
				}
				if (outLine.size() > matchString.size())
				{
					bool matching = true;
					for (int m = 0; m < matchString.size(); m++)
					{
						if (outLine[m] != matchString[m]) matching = false;
					}
					if (matching) finishedStructure = true;//cannot be made false
				}
				if (outLine.size() > eMatch.size())
				{
					bool matching = true;
					for (int m = 0; m < eMatch.size(); m++)
					{
						if (outLine[m] != eMatch[m]) matching = false;
					}
					if (matching) finishedEnergy = true;//cannot be made false
				}
				if (outLine.size() > fail.size())
				{
					bool matching = true;
					for (int m = 0; m < fail.size(); m++)
					{
						if (outLine[m] != fail[m]) matching = false;
					}
					if (matching) {
						std::cout << "Optimization failed due to a problem in gaussian. error code: " << outLine << std::endl;
						//std::cout << "finished output but did not finish finished output: " << finishedOutput << " , finished energy: " << finishedEnergy << " finished structure:" << finishedStructure << std::endl;

						if (!timedOut && (!finishedEnergy || !finishedStructure))return nullptr;
						else std::cout << "continuing as calculation was still successful" << std::endl;
					}
				}
			}


		std::cout << "finished output: " << finishedOutput << std::endl;
		

		if (!finishedEnergy || !finishedOutput)
		{
			return nullptr;//otherwise we at least can get an energy and geometry
		}

		aOutFile->close();
		delete aOutFile;
		std::ifstream* gOutFile = new std::ifstream(dirname + "/" + dftFileName + oFT);
		std::string cline;
		structure* st = nullptr;

		matchString = "Input orientation:";
		eMatch = " SCF Done:  E";
		std::string energyLine = "";

		std::cout << "reading file" << std::endl;

		while (getline(*gOutFile, cline)) {
			if (cline.length() > matchString.length())
			{

				bool matching = false;
				int matchIt = 0;
				for (int l = 0; l < cline.length() - matchString.length(); l++)
				{
					if (cline[l] == matchString[matchIt])
					{
						matchIt++;
						matching = true;
					}
					else if (matchIt < 18) {
						matchIt = 0;
						matching = false;
					}
					else {
						//the strings are already matched let the match go through
					}
				}
				//if (matching) std::cout << "matching ends check length: " << matchIt << " " << matchString.length() << std::endl;
				if (matching && matchIt >= matchString.length())
				{
					std::vector<std::vector<atom>> set;
					std::vector<std::string> elements;
					getline(*gOutFile, cline);
					getline(*gOutFile, cline);
					getline(*gOutFile, cline);
					getline(*gOutFile, cline);
					char space = ' ';
					while (getline(*gOutFile, cline) && cline[1] != '-')
					{
						std::string cnumber = "";
						int num = 0;
						std::string element;
						double x = 0;
						double y = 0;
						double z = 0;
						//std::cout << "cline: " << cline << std::endl;
						for (int i = 0; i < cline.length(); i++)
						{
							if (cline[i] == space)
							{
								//std::cout << "space detected, cnumber: " << cnumber << std::endl;
								if (cnumber.length() > 0)
								{
									num++;
									switch (num)
									{
									case 1:
										//center number, not sure what this is for. do nothing
										break;
									case 2:
										//atomic number
										element = atomicNumber(std::stoi(cnumber));
										//std::cout << "element: " << element << std::endl;
										break;
									case 3:
										//atomic type, not sure what this is. do nothing.
										break;
									case 4:
										//x coordinate
										x = std::stof(cnumber);
										//std::cout << "x: " << x << std::endl;
										break;
									case 5:
										//y
										y = std::stof(cnumber);
										//std::cout << "y: " << y << std::endl;
										break;
									case 6:
										//z
										z = std::stof(cnumber);
										//std::cout << "z: " << z << std::endl;
										break;
									default:
										break;

									}
									cnumber = "";
								}
							}
							else {
								cnumber += cline[i];
							}
						}
						//make sure z is added
						z = std::stof(cnumber);
						//std::cout << "z: " << z << std::endl;
						int index = -1;
						for (int i = 0; i < elements.size(); i++)
						{
							if (elements[i] == element)
							{
								index = i;
							}
						}
						if (index == -1)
						{
							elements.push_back(element);
							std::vector < atom> alist = { atom(x,y,z) };
							set.push_back(alist);
						}
						else {
							set[index].push_back(atom(x, y, z));
						}
					}
					if (st != nullptr) delete st;
					st = new structure(set, elements);


				}
			}
			if (cline.size() > eMatch.size())
			{
				bool matching = true;
				for (int m = 0; m < eMatch.size(); m++)
				{
					if (cline[m] != eMatch[m]) matching = false;
				}
				if (matching) energyLine = cline;
			}


		}
		gOutFile->close();
		delete gOutFile;

		//find the start of the energy
		int eStart = 0;
		int eFin = 0;
		bool equalPassed = false;
		for (int i = 0; i < energyLine.size(); i++)
		{
			if (eFin) {}
			else if (eStart) {
				if (energyLine[i] == ' ')
				{
					eFin = i;
				}
			}
			else if (equalPassed) {
				if (energyLine[i] == ' ') {}
				else eStart = i;
			}
			else if (energyLine[i] == '=') equalPassed = true;

		}
		//if(DEBUG > 1) std::cout << energyLine.substr(eStart, eFin - eStart) << std::endl;
		double energyHartree = std::stod(energyLine.substr(eStart, eFin - eStart));
		st->energy = energyHartree;
		return st;
	}
	return nullptr;
}

