// declare all functions here

#pragma once

#include <iostream>
#include "f2c.h"
#include "clapack.h"
#include "physconst.h"
#include <fstream>
#include "structdef.h"


// convert bandparameters into SI units
void BandParaSIunits(BandParameters A, BandParameters &A_SI);

void SolveLinEq(double ** A, double *b, int Nd, double *X);

void ZSolveLinEq(doublecomplex ** A, doublecomplex *b, int Nd, doublecomplex *X);

double my_interp_eqdx(double x[], double y[], int N, double x0);

double my_interp_neqdx(double x[], double y[], int N, double x0);

int findlayer(double z, double * border, int Nlayer);

void GetNormalized(double WaveFunc[], int Nsize);

double factorial (int num);

double HermiteH(int n, double x);

double phiHarmonic(int n, double aH, double x);

void fsaveVqform(char fVqform_name[], double * qform_array, double * Vqform_array, int Nqform);
	
void freadVqform(char fVqform_name[], double * &qform_array, double * &Vqform_array, int &Nqform);

void SetBandParaGaAs(BandParameters &GaAs, double T_lattice);

void SetBandParaInAs(BandParameters &InAs, double T_lattice);

void SetBandParaGaInAs(BandParameters &GaInAs, BandParameters GaAs, BandParameters InAs, double T_lattice, double x_GaInAs);

void GetStrain(BandParameters Barrier, BandParameters Well, BandParameters Substrate, Strain &strB, Strain &strW);

double effmzE(BandParameters Material, Strain strM, double DelE);

double effmxyE(BandParameters Material, Strain strM, double DelE);

double effmzHH(BandParameters Material, Strain strM, double DelE);

double effmxyHH(BandParameters Material, Strain strM, double DelE);

double effmzLH(BandParameters Material, Strain strM, double DelE);

double effmxyLH(BandParameters Material, Strain strM, double DelE);

void getBandStruct(BandParameters Barrier, BandParameters Well, Strain strB, Strain strW,
	double * layer, bool * isbarrier, double * border, int Nlayer, double * z_grid, double * V_grid, double * z_grid_halfdzshift, double * V_grid_halfdzshift, int N_grid, 
	double (*effmzEdep)(BandParameters, Strain, double delE), double (*effmxyEdep)(BandParameters, Strain, double delE),
	int N_level, int &nEbound, double * eigE, double ** eigV, double * effmxy_ave);

double CoulombHolePP(double n2D, double epsr0, double effmR, double kappa, double Cq);

double n2DEF(double * k_array, double * Ek_array, int Nk, double mu, double T);

double n2DEF_dxdy(double * k_array, double * Ek_array, int Nk, double mu, double T);

double FermiLevel(double * k_array, double * Ek_array, int Nk, double n2D, double T);

double Polarizability_normal(double * k_array, double * Ek_array, int Nk, double n2D, double T, double q);

double Polarizability(double * k_array, double * Ek_array, int Nk, double n2D, double T, double q);

double PolTeq0(double effm, double mu, double q);

double PolT(double effm, double mu, double T, double q);

void fsaveVector(char fVname[], double * V, int N);