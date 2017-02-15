#pragma once 

struct BandParameters 
{
	double Eg;  // bandgap at T = 0, in eV
	double Ep;  // Kane energy, in eV
	double VBO; // valence band offset, in unit eV
	double KaneF;  // Kane parameter(F)
	double DeltaSO;  // Spin-Orbit Split Off in eV
	double dielectric0;        // static dielectric constant, Room Temperature, 300 K
	double dielectricinf;       // high frequency dielectric constant, Room Temperature
	double E_LO;  // LO phonon energy, in eV
	double gamma1, gamma2, gamma3;  // luttinger parameters
	// ***************** strain parameters below *****************
	double lc;  // lattice consant, in angstrom
	// hydrostatic & shear deformation potential
	double ac; // eV
	double av; // eV
	double b;  // eV
	double d;  // eV
	double C11; // GPa
	double C12; // GPa
	double C44; // GPa
};


struct Strain
{
	double exx, eyy, ezz;
	double cPe, vPe, Qe, Ehh, Elh, Eso, alphamix, betamix, Amix, Bmix, Cmix;
};