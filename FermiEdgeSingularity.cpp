#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <new>
#include <algorithm>
#include <stdlib.h>
#include "physconst.h"
#include "f2c.h"
#include "clapack.h"
#include "functions.h"


char StructureName[50] = "FES";

//double n2D = 0.5e12; // in cm^(-2)

double kmaxkF = 15.0;
//bool Vs_calculated = 0; // = 1 if calculated.
bool Vqform_calculated = 1;

double gammaBR = 2e-3*echarge;

double T_lattice = 5;  // lattice temperature, in K
double T_electron = 5;  // electronic temperature, this program is valid for low T only.

double Lw = 8e-9;  // well width, in meter

int Nlayer1p = 2;
// 1 period of structure
// layer thickness, in Angstrom
double layer1p[] = {150,80};
// isbarrier1p = 1 for barrier, = 0 for well
bool isbarrier1p[] = {1,0};
// doping density, in cm^(-3)
double ndope1p[] = {0,0};

double Eext = 0;  // external electric field, in kV/cm


double Te_eq_0_tr = 1; // if T_electron is smaller than this, then Te = 0 is assumed.


BandParameters Well, Barrier, Substrate;  // global, all quantities are in SI units in these three structs, and they will be directly used.


int main()
{
	using namespace std;
	
	// change the layer thickness and doping density into SI units
	for (int kk = 0; kk < Nlayer1p; kk++)
	{
		layer1p[kk] = layer1p[kk]*1e-10;
		ndope1p[kk] = ndope1p[kk]*1e6;
	}

	Eext = Eext*1e5;

	//n2D = n2D*1e4;

	BandParameters GaAs, InAs, GaInAs;

	SetBandParaGaAs(GaAs, T_lattice);

	SetBandParaInAs(InAs, T_lattice);

	double x_GaInAs = 0.20;  // Ga_(1-x)In_(x)As
	SetBandParaGaInAs(GaInAs, GaAs, InAs, T_lattice, x_GaInAs);
	
	Well = GaInAs;
	Barrier = GaAs;
	Substrate = GaAs;
	BandParaSIunits(GaInAs, Well); 
	BandParaSIunits(GaAs, Barrier);
	BandParaSIunits(GaAs, Substrate);

	Strain strB, strW;  // the strain for barrier and well

	GetStrain(Barrier, Well, Substrate, strB, strW);

	cout << effmxyE(Barrier, strB, 0) << endl;
	cout << effmxyE(Well, strW, 0) << endl;
		
	double Te;
	Te = T_electron;

	double epsr0 = Well.dielectric0;
	double dip = eGauss*hbar/Well.Eg*sqrt(Well.Ep/(2*m0));

	const bool twoperiods = 0;  // = 1 for calculating in two periods, = 0 for one period
	const double dz = 0.1e-10;  // step length, in m
	const int Ne_level = 15;  // number of levels to calculate in conduction band
	const int Nhh_level = 15;  // number of levels to calculate in heavy hole band

	// *********** specify the structure begin *********** 
	double Lperiod;
	Lperiod = 0;
	for (int kk = 0; kk < Nlayer1p; kk++)
		Lperiod += layer1p[kk];
	int Nperiod = static_cast<int>(floor(Lperiod/dz));
			
	double n2D1p;
	n2D1p = 0;
	for (int kk = 0; kk < Nlayer1p; kk++)
		n2D1p += ndope1p[kk]*layer1p[kk];

	const int Nleadingbarrier = 1;
	int Nlayer;
	Nlayer = (1+static_cast<int>(twoperiods))*Nlayer1p + Nleadingbarrier;
	
	double * layer = new double[Nlayer];
	bool * isbarrier = new bool[Nlayer];
	double * border = new double[Nlayer];

	if (twoperiods)
	{		
		for (int kk = 0; kk < Nlayer1p; kk++)
		{
			layer[kk] = layer1p[kk];
			layer[kk+Nlayer1p] = layer1p[kk];
			isbarrier[kk] = isbarrier1p[kk];
			isbarrier[kk+Nlayer1p] = isbarrier1p[kk];
		}
		for (int kk = 0; kk < Nleadingbarrier; kk++)
		{
			layer[2*Nlayer1p+kk] = layer1p[kk];
			isbarrier[2*Nlayer1p+kk] = isbarrier1p[kk];
		}
	}
	else
	{
		for (int kk = 0; kk < Nlayer1p; kk++)
		{
			layer[kk] = layer1p[kk];
			isbarrier[kk] = isbarrier1p[kk];
		}
		for (int kk = 0; kk < Nleadingbarrier; kk++)
		{
			layer[Nlayer1p+kk] = layer1p[kk];
			isbarrier[Nlayer1p+kk] = isbarrier1p[kk];
		}
	}		
	border[0] = layer[0];
	for (int kk = 1; kk < Nlayer; kk++)
		border[kk] = border[kk-1] + layer[kk];

	double L_structure = border[Nlayer-1];	
	// *********** specify the structure end *********** 

	int N_grid = static_cast<int>(floor(L_structure/dz));

	double * z_grid = new double [N_grid];
	double * z_grid_halfdzshift = new double[N_grid];
	z_grid[0] = 0;
	for (int kg = 1; kg < N_grid; kg++)
		z_grid[kg] = z_grid[kg-1] + dz;

	for (int kg = 0; kg < N_grid; kg++)
		z_grid_halfdzshift[kg] = z_grid[kg] + 0.5*dz;

	double * Ve_grid = new double [N_grid];
	double * Ve_grid_halfdzshift = new double[N_grid];
	for (int kg = 0; kg < N_grid; kg++)
	{
		int klayer;
		klayer = findlayer(z_grid[kg], border, Nlayer);
		if (isbarrier[klayer])
			Ve_grid[kg] = (Barrier.VBO + Barrier.Eg + strB.cPe) - echarge*Eext*z_grid[kg];  // the potential due to Eext is set to be 0 at z = 0 
		else
			Ve_grid[kg] = (Well.VBO + Well.Eg + strW.cPe) - echarge*Eext*z_grid[kg];
	}
	for (int kg = 0; kg < N_grid; kg++)
	{
		int klayer;
		klayer = findlayer(z_grid_halfdzshift[kg], border, Nlayer);
		if (isbarrier[klayer])
			Ve_grid_halfdzshift[kg] = (Barrier.VBO + Barrier.Eg + strB.cPe) - echarge*Eext*z_grid_halfdzshift[kg];
		else
			Ve_grid_halfdzshift[kg] = (Well.VBO + Well.Eg + strW.cPe) - echarge*Eext*z_grid_halfdzshift[kg];
	}

	double * Vhh_grid = new double [N_grid];
	double * Vhh_grid_halfdzshift = new double[N_grid];
	for (int kg = 0; kg < N_grid; kg++)
	{
		int klayer;
		klayer = findlayer(z_grid[kg], border, Nlayer);
		if (isbarrier[klayer])
			Vhh_grid[kg] = - ((Barrier.VBO + strB.Ehh) - echarge*Eext*z_grid[kg]);  // the potential due to Eext is set to be 0 at z = 0 
		else
			Vhh_grid[kg] = - ((Well.VBO + strW.Ehh) - echarge*Eext*z_grid[kg]);
	}
	for (int kg = 0; kg < N_grid; kg++)
	{
		int klayer;
		klayer = findlayer(z_grid_halfdzshift[kg], border, Nlayer);
		if (isbarrier[klayer])
			Vhh_grid_halfdzshift[kg] = - ((Barrier.VBO + strB.Ehh) - echarge*Eext*z_grid_halfdzshift[kg]);
		else
			Vhh_grid_halfdzshift[kg] = - ((Well.VBO + strW.Ehh) - echarge*Eext*z_grid_halfdzshift[kg]);
	}

	int nEe;
	double * eigEe, ** eigVe, * effmxye_ave;
	eigEe = new double [Ne_level];
	eigVe = new double * [Ne_level];
	effmxye_ave = new double [Ne_level];

	
	getBandStruct(Barrier, Well, strB, strW, layer, isbarrier, border, Nlayer, z_grid, Ve_grid, z_grid_halfdzshift, Ve_grid_halfdzshift, N_grid, 
	effmzE, effmxyE,
	Ne_level, nEe, eigEe, eigVe, effmxye_ave);

	cout << "nEe:\t" << nEe << endl;

	int nEhh;
	double * eigEhh, * eigEhhv, ** eigVhh, * effmxyhh_ave;
	eigEhh = new double [Nhh_level];
	eigEhhv = new double [Nhh_level];
	eigVhh = new double * [Nhh_level];
	effmxyhh_ave = new double [Nhh_level];
	getBandStruct(Barrier, Well, strB, strW, layer, isbarrier, border, Nlayer, z_grid, Vhh_grid, z_grid_halfdzshift, Vhh_grid_halfdzshift, N_grid, 
	effmzHH, effmxyHH,
	Nhh_level, nEhh, eigEhh, eigVhh, effmxyhh_ave);
	cout << "nEhh:\t" << nEhh << endl;
	for (int kn = 0; kn < nEhh; kn++)
		eigEhhv[kn] = - eigEhh[kn];

	double effmE, effmH, effmR;
	effmE = effmxye_ave[0];
	effmH = effmxyhh_ave[0];
	effmR = 1.0/(1.0/effmE+1.0/effmH);
	cout << "effmE, effmH:\t" << effmE/m0 << "\t" << effmH/m0 << "\t" << effmR/m0 << endl;

	double effmEz, effmHz;
	effmEz = effmzE(Well, strW, 0.0);
	effmHz = effmzHH(Well, strW, 0.0);
	cout << "Ee, Eh:\t" << pow(PI*hbar/Lw,2.0)/(2.0*effmEz)/echarge*1000.0 << "\t" << pow(PI*hbar/Lw,2.0)/(2.0*effmHz)/echarge*1000.0 << endl;


	double Egeff0;
	//Egeff0 = eigEe[0] - eigEhhv[0];
	Egeff0 = 1.338*echarge;  // got from the absorption fitting data
	cout << Egeff0/echarge << endl;	
	//getchar();

	//cout << "LO: " << hbar*sqrt(2*PI*pow(eGauss,2.0)*n2D*(1.0/Lw)/(epsr0*effmR))/echarge*1000 << endl;

	
	char fVqform_name[200];
	sprintf_s(fVqform_name, 200, "%s_Tl%.0f_Vqfrom.txt",StructureName,T_lattice);

	double * qform_array, * Vqform_array;
	int Nqform;

	// calculating the form factor in V(q)
	if (!Vqform_calculated)
	{
		cout << "calculating the form factor in V(q) ...\n"; 
		double qform_max = 2.05*5.0*sqrt(2.0*PI*5e12*1e4);
		Nqform = 5000;
		double dqform = qform_max/static_cast<double>(Nqform-1);
		qform_array = new double[Nqform];
		Vqform_array = new double[Nqform];
		for (int nq = 0; nq < Nqform; nq++)
		{
			double qq = static_cast<double>(nq)*dqform;
			double intze = 0;
			for (int kge = 0; kge < N_grid; kge++)
			{
				double ze, wavefe;
				ze = z_grid[kge];
				wavefe = eigVe[0][kge];
				double intzh = 0;
				for (int kgh = 0; kgh < N_grid; kgh++)
				{
					double zh, wavefh;
					zh = z_grid[kgh];
					wavefh = eigVhh[0][kgh];
					double addintzh = wavefe*wavefe*wavefh*wavefh*exp(-qq*fabs(ze-zh));
					intzh += addintzh;
				}
				double addintze = intzh;
				intze += addintze;
			}
			qform_array[nq] = qq;
			Vqform_array[nq] = intze;
			cout << nq << "\t" << Vqform_array[nq] << "\t";
		}
		cout << endl;
		fsaveVqform(fVqform_name, qform_array, Vqform_array, Nqform);
	}
	else
	{
		qform_array = NULL;
		Vqform_array = NULL;
		//cout << Vqform_array << endl;
		freadVqform(fVqform_name, qform_array, Vqform_array, Nqform);
		//cout << Vqform_array << endl;
	}	
	//for (int nq = 0; nq < Nq; nq++)
		//cout << q_array[nq]*Lw << "\t" << Vqform_array[nq] << endl;

	double n2Dmin = 1.0e11;
	double n2Dmax = 5.0e12;
	double dn2D = 1.0e11;
	int Nn2D = static_cast<int>((n2Dmax-n2Dmin)/dn2D)+1;
	for (int kn = 0; kn < Nn2D; kn++)
	{
		double n2D = n2Dmin + static_cast<double>(kn)*dn2D;
		cout << "\nn2D: " << n2D << endl;

		n2D *= 1.0e4;

	double kF = sqrt(2.0*PI*n2D);	

	cout << "kF\t" << kF << endl;

	double qmax;
	//qmax = 3.0*kmaxkF*kF;
	qmax = 3.0*sqrt(2.0*PI*5e12*1e4);
	int Nq = 10000;
	double dq = qmax/static_cast<double>(Nq-1);
	double * q_array = new double [Nq];
	double * PolEq_array = new double [Nq];
	double * PolHq_array = new double [Nq];
	double * Vsq_array = new double[Nq];
	double * VsPPq_array = new double[Nq];
	
	for (int nq = 0; nq < Nq; nq++)
	{
		double qq = static_cast<double>(nq)*dq;
		q_array[nq] = qq;
	}
	
	cout << "density of states (e, h):\t" << effmE/(PI*hbar*hbar) << "\t" << effmH/(PI*hbar*hbar) << endl;
	
	double Vqcoef = 2.0*PI*pow(eGauss,2.0)/epsr0;

	double muE0, muH0;  // chemical potential
	muE0 = kB*Te*log(exp(PI*hbar*hbar*n2D/(effmE*kB*Te))-1);
	muH0 = kB*Te*log(exp(PI*hbar*hbar*n2D/(effmH*kB*Te))-1);

		

	for (int nq = 0; nq < Nq; nq++)
	{
		double qq;
		qq = q_array[nq];
		PolEq_array[nq] = PolT(effmE, muE0, Te, qq);
		PolHq_array[nq] = PolT(effmH, muH0, Te, qq);
		//cout << q_array[nq]/kF << "\t" << PolEq_array[nq] << "\t"  << PolHq_array[nq] << endl;
	}

	/*

	double kappaE, kappaH, kappa;
	if (Te < Te_eq_0_tr)
	{
		kappaE = 2*effmE*pow(eGauss,2.0)/(epsr0*hbar*hbar);
		kappaH = 2*effmH*pow(eGauss,2.0)/(epsr0*hbar*hbar);
	}
	else
	{
		kappaE = 2*effmE*pow(eGauss,2.0)/(epsr0*hbar*hbar)*(1-exp(-PI*hbar*hbar*n2D/(effmE*kB*Te)));
		kappaH = 2*effmH*pow(eGauss,2.0)/(epsr0*hbar*hbar)*(1-exp(-PI*hbar*hbar*n2D/(effmH*kB*Te)));
	}
	kappa = kappaE + kappaH;

	// plasmon-pole approximation
	double Cq = 1.0;
	for (int nq = 0; nq < Nq; nq++)
	{
		double qq;
		qq = q_array[nq];
		double wplsqcoef, Cqcoef, Vqform;
		wplsqcoef = 2.0*PI*pow(eGauss,2.0)*n2D/(epsr0*effmR);
		Cqcoef = (Cq/4.0)*pow(hbar/(2*effmR),2.0)/wplsqcoef;
		Vqform = my_interp_eqdx(qform_array, Vqform_array, Nqform, qq);
		//Vqform = 1.0;
		VsPPq_array[nq] = Vqcoef*Vqform*(1.0/kappa+Cqcoef*pow(qq,2.0))/(1.0+qq/kappa+Cqcoef*pow(qq,3.0));		
	}	
	
	
	double DelECH_PP;  // Coulomb hole energy
	DelECH_PP = CoulombHolePP(n2D, epsr0, effmR, kappa, Cq);
	cout << "Coulomb hole (PP): " << DelECH_PP/echarge*1000 << endl;

	// calculate Coulomb hole
	double qCHmax = 5.0*kF;
	int NqCH = 10000;
	double dqCH = qCHmax/static_cast<double>(NqCH-1);
	double intqCH = 0;
	for (int nq = 0; nq < NqCH; nq++)
	{
		double qCH = static_cast<double>(nq)*dqCH;
		double Formq, PolEq, PolHq, Polq;
		Formq = my_interp_eqdx(qform_array, Vqform_array, Nqform, qCH);
		PolEq = my_interp_eqdx(q_array, PolEq_array, Nq, qCH);
		PolHq = my_interp_eqdx(q_array, PolHq_array, Nq, qCH);
		Polq = PolEq + PolHq;
		double addintqCH = Formq*Formq*Polq/(qCH-Vqcoef*Formq*Polq);
		if (nq == 0 || nq == NqCH-1)
			addintqCH *= 0.5;
		intqCH += addintqCH*dqCH;
	}
	double DelECH;
	DelECH = intqCH*Vqcoef*Vqcoef/(2.0*PI);
	cout << "Coulomb hole: " << DelECH/echarge*1000 << endl;
	*/


	double Egeff;
	//Egeff = Egeff0 + DelECH;
	Egeff = Egeff0;

	for (int nq = 0; nq < Nq; nq++)
	{
		double qq;

		qq = q_array[nq];
		double Vqform;
		Vqform = my_interp_eqdx(qform_array, Vqform_array, Nqform, qq);
		//Vqform = 1.0;
		double polE, polH;
		polE = PolEq_array[nq];
		polH = PolHq_array[nq];
		Vsq_array[nq] = Vqcoef*Vqform/(qq-Vqcoef*Vqform*(polE+polH));
	}	

	
	cout << "Calculating nonparabolic E(k) ...\n\n";
	double kkmax;
	kkmax = qmax/2.05;
	//kkmax = (kmaxkF+2.0)*kF;
	//kkmax = 15.0*kF;
	int Nkk = 5000;
	double dkk = kkmax/static_cast<double>(Nkk-1);
	double * kk_array = new double [Nkk];
	double * Eekk_array = new double [Nkk];
	double * Ehkk_array = new double [Nkk];	

	for (int nkk = 0; nkk < Nkk; nkk++)
	{
		double kk;
		kk = static_cast<double>(nkk)*dkk;
		kk_array[nkk] = kk;
	}
	
	for (int nkk = 0; nkk < Nkk; nkk++)
	{
		double kk = kk_array[nkk];

		if (kk==0)
			Eekk_array[nkk] = 0.0;
		else
		{
			double effm0, E1, E2, effmE1;
			{
				double E_temp = eigEe[0];
				double effmxy_ave_inverse = 0;
				for (int kg = 0; kg < N_grid; kg++)
				{
					int klayer;
					double add_effmxy_ave_inverse;
					klayer = findlayer(z_grid[kg], border, Nlayer);
					if (isbarrier[klayer])
						add_effmxy_ave_inverse = pow(eigVe[0][kg],2.0)/effmxyE(Barrier, strB, E_temp - Ve_grid[kg]);
					else
						add_effmxy_ave_inverse = pow(eigVe[0][kg],2.0)/effmxyE(Well, strW, E_temp - Ve_grid[kg]);
					if (add_effmxy_ave_inverse < 0)
					{
						cout << add_effmxy_ave_inverse << endl;
						cerr << "Negative in-plane effective mass occured!" << endl;
						exit(1);
					}
					effmxy_ave_inverse += add_effmxy_ave_inverse;
				}
				effm0 = 1.0/effmxy_ave_inverse;
			}
			E1 = hbar*hbar*kk*kk/(2.0*effm0);
			do
			{
				E2 = E1;
				{
					double E_temp = E1 + eigEe[0];
					double effmxy_ave_inverse = 0;
					for (int kg = 0; kg < N_grid; kg++)
					{
						int klayer;
						double add_effmxy_ave_inverse;
						klayer = findlayer(z_grid[kg], border, Nlayer);
						if (isbarrier[klayer])
							add_effmxy_ave_inverse = pow(eigVe[0][kg],2.0)/effmxyE(Barrier, strB, E_temp - Ve_grid[kg]);
						else
							add_effmxy_ave_inverse = pow(eigVe[0][kg],2.0)/effmxyE(Well, strW, E_temp - Ve_grid[kg]);
						if (add_effmxy_ave_inverse < 0)
						{
							cout << E_temp - Ve_grid[kg] << endl;
							cout << effmxyE(Barrier, strB, E_temp - Ve_grid[kg])/m0 << endl;
							cout << effmxyE(Well, strW, E_temp - Ve_grid[kg])/m0 << endl;
							cout << add_effmxy_ave_inverse/pow(eigVe[0][kg],2.0) << endl;
							cerr << "Negative in-plane effective mass occured!" << endl;
							exit(1);
						}
						effmxy_ave_inverse += add_effmxy_ave_inverse;
					}
					effmE1 = 1.0/effmxy_ave_inverse;
				}
				E1 = hbar*hbar*kk*kk/(2.0*effmE1);
			}while(fabs((E2-E1)/E1)>1e-5);
			Eekk_array[nkk] = E1;
		}		
	}

	//cout << "bb\t" << Eekk_array[Nkk-1]/(hbar*hbar*pow(kk_array[Nkk-1],2.0)/(2.0*effmxyE(Well, strW, 0.0))) << endl;
	//cout << kk_array[Nkk-1]/kF << endl;
	//cout << Eekk_array[Nkk-1]/echarge*1000 << endl;

	for (int nkk = 0; nkk < Nkk; nkk++)
	{
		double kk = kk_array[nkk];
		Ehkk_array[nkk] = hbar*hbar*kk*kk/(2.0*effmH);
	}
	
	
	
/*	for (int nkk = 1; nkk < Nkk; nkk++)
	{
		double kk;
		kk = static_cast<double>(nkk)*dkk;
		double aaa = Eekk_array[nkk];
		double bbb = hbar*hbar*kk*kk/(2.0*effmE);
		printf("%10.6e\t",((bbb-aaa)/bbb));
	}
	cout << endl;
*/	
	cout << "Calculating exchange self energy ...\n\n";
	double * Eexkk_array = new double[Nkk];
	for (int nkk = 0; nkk < Nkk; nkk++)
	{
		double kk;
		kk = kk_array[nkk];
		int Nkp = 200;
		double dkp = kF/static_cast<double>(Nkp-1);  // assuming T = 0, so it's only valid for low T
		double intkp = 0.0;
		for (int nkp = 0; nkp < Nkp; nkp++)
		{
			double addintkp;
			double kp = static_cast<double>(nkp)*dkp;
			int Ntheta = 50;
			double dtheta = PI/static_cast<double>(Ntheta-1);
			double inttheta = 0.0;
			for (int nth = 0; nth < Ntheta; nth++)
			{
				double addintth;
				double theta = static_cast<double>(nth)*dtheta;
				double k_kp = sqrt(kk*kk+kp*kp-2.0*kk*kp*cos(theta));
				addintth = my_interp_eqdx(q_array,Vsq_array,Nq,k_kp)*dtheta;
				if (nth == 0 || nth == Ntheta-1)
					addintth *= 0.5;
				inttheta += addintth;
			}
			addintkp = 2.0*inttheta*kp*dkp;  // times 2.0 since we only integrate from 0 to Pi
			if (nkp == 0 || nkp == Nkp-1)
				addintkp *= 0.5;
			intkp += addintkp;
		}
		Eexkk_array[nkk] = intkp/pow(2.0*PI,2.0);
	}

	double EexkF = my_interp_eqdx(kk_array,Eexkk_array,Nkk,kF);
	cout << "&&&&&&&&&&&&\t" << Eexkk_array[0] << "\t" << EexkF <<"\t" << Eexkk_array[Nkk-1] << endl;
		
	double * EeRkk_array = new double[Nkk];
	double * EhRkk_array = new double[Nkk];

	for (int nkk = 0; nkk < Nkk; nkk++)
	{
		double kk = kk_array[nkk];
		EeRkk_array[nkk] = Eekk_array[nkk] - Eexkk_array[nkk];
		EhRkk_array[nkk] = Ehkk_array[nkk] - Eexkk_array[nkk];
	}

	
	double muE, muH;
	//muE = my_interp_eqdx(kk_array,EeRkk_array,Nkk,kF);
	//cout << hbar*hbar*kF*kF/(2.0*effmE)/echarge*1000 << endl;
	//cout << my_interp_eqdx(kk_array,Eekk_array,Nkk,kF)/echarge*1000 << endl;
	//cout << muE/echarge*1000 << endl;
	//muH = my_interp_eqdx(kk_array,EhRkk_array,Nkk,kF);
	//cout << "============\t" << muE << "\t" << muH << endl;

	muE = FermiLevel(kk_array, EeRkk_array, Nkk, n2D, Te);
	muH = FermiLevel(kk_array, EhRkk_array, Nkk, n2D, Te);
	
	//scout << "============\t" << muE << "\t" << muH << endl;
	//cout << "kB*T:\t" << kB*Te/echarge*1000 << endl;
	//cout << (muE-EeRkk_array[0])/echarge*1000 << endl;
	//cout << (muH-EhRkk_array[0])/echarge*1000 << endl;


	double kmin, kmax;
	kmin = 0.0*kF;
	//kmax = kmaxkF*kF;
	kmax = kkmax;
	int Nk = 1000;
	double dk = (kmax-kmin)/static_cast<double>(Nk-1);
	double * k_array = new double [Nk];
	for (int nk = 0; nk < Nk; nk++)
		k_array[nk] =  kmin + static_cast<double>(nk)*dk;


	double ** Vskkp = new double * [Nk];
	for (int nk = 0; nk < Nk; nk++)
		Vskkp[nk] = new double [Nk];
	
	cout << "Calculating Vskkp[nk1][nk2] ... \n\n";
	for (int nk1 = 0; nk1 < Nk; nk1++)
	{
		for (int nk2 = 0; nk2 < Nk; nk2++)
		{
			double k1, k2;
			k1 = k_array[nk1];
			k2 = k_array[nk2];
			int Ntheta = 200;
			double dtheta = PI/static_cast<double>(Ntheta-1);
			double inttheta = 0;
			for (int nth = 0; nth < Ntheta; nth++)
			{
				double theta = static_cast<double>(nth)*dtheta;
				double qq = sqrt(k1*k1+k2*k2-2*k1*k2*cos(theta));
				if (qq<q_array[0] || qq>q_array[Nq-1])
				{
					cout << "qq out of range\n";
					exit(1);
				}
				double Vsq = my_interp_eqdx(q_array, Vsq_array, Nq, qq);
				
				double addint = Vsq*dtheta;
				if (nth == 0 || nth == Ntheta-1)
					addint *= 0.5;
				inttheta += addint;
			}
			Vskkp[nk1][nk2] = 2.0*inttheta;  // 2.0 since we only integrate from 0 to Pi
		}
	}

	
	double * fFDE = new double [Nk];
	double * fFDH = new double [Nk];
	for (int nk = 0; nk < Nk; nk++)
	{
		double kk, Eek, Ehk;
		kk = k_array[nk];
		Eek = my_interp_eqdx(kk_array,EeRkk_array,Nkk,kk);
		Ehk = my_interp_eqdx(kk_array,EhRkk_array,Nkk,kk);
		if (Te < Te_eq_0_tr)
		{
			if (kk <= kF)
			{
				fFDE[nk] = 1.0;
				fFDH[nk] = 1.0;
			}
			else
			{
				fFDE[nk] = 0.0;
				fFDH[nk] = 0.0;
			}
		}
		else
		{
			fFDE[nk] = 1.0/(exp((Eek-muE)/(kB*Te))+1.0);
			fFDH[nk] = 1.0/(exp((Ehk-muH)/(kB*Te))+1.0);
		}
	}


	double * dipk_array = new double [Nk];
	for (int nk = 0; nk < Nk; nk++)
	{
		double kk, Eek, Ehk, Ee0, Eh0;
		kk = k_array[nk];
		Eek = my_interp_eqdx(kk_array,EeRkk_array,Nkk,kk);
		Ehk = my_interp_eqdx(kk_array,EhRkk_array,Nkk,kk);
		Ee0 = my_interp_eqdx(kk_array,EeRkk_array,Nkk,0.0);
		Eh0 = my_interp_eqdx(kk_array,EhRkk_array,Nkk,0.0);
		dipk_array[nk] = dip*(Egeff+Ee0+Eh0)/(Egeff+Eek+Ehk);
	}

	//cout << "dipk:\t" << dipk_array[0] << "\t" << dipk_array[Nk-1] << endl;

	
	double Ephmin, Ephmax;
	Ephmin = Egeff+EeRkk_array[0]+EhRkk_array[0] - 0.2*(muE-EeRkk_array[0]+muH-EhRkk_array[0]);
	Ephmax = Egeff+EeRkk_array[0]+EhRkk_array[0] + 2.0*(muE-EeRkk_array[0]+muH-EhRkk_array[0]);
	//Ephmin = 1.26*echarge;
	//Ephmax = 1.4*echarge;


	cout << "Calculating susceptibility ...\n\n";
	int Nph = 200;
	double dEph = (Ephmax - Ephmin)/static_cast<double>(Nph-1);
	double * Eph = new double [Nph];
	double * GainPH = new double [Nph];
	double * GainPH0 = new double [Nph];
	for (int np = 0; np < Nph; np++)
		Eph[np] = Ephmin + static_cast<double>(np)*dEph;

	for (int np = 0; np < Nph; np++)
	{
		cout << np << "\t";
		double Ephoton = Eph[np];

		double * gammak_array = new double [Nk];
		double Ealpha = 100e-3*echarge;
		for (int nk = 0; nk < Nk; nk++)
		{
			double kk, Eek, Ehk;
			kk = k_array[nk];
			Eek = my_interp_eqdx(kk_array,Eekk_array,Nkk,kk);
			Ehk = my_interp_eqdx(kk_array,Ehkk_array,Nkk,kk);
			gammak_array[nk] = gammaBR;
		}

		doublecomplex * chi = new doublecomplex [Nk];
		doublecomplex * chi0 = new doublecomplex [Nk];
		for (int nk = 0; nk < Nk; nk++)
		{
			double kk, Eek, Ehk, gammak;
			kk = k_array[nk];
			gammak = gammak_array[nk];
			Eek = my_interp_eqdx(kk_array,EeRkk_array,Nkk,kk);
			Ehk = my_interp_eqdx(kk_array,EhRkk_array,Nkk,kk);
			chi0[nk].r = dipk_array[nk]*(fFDE[nk]+fFDH[nk]-1.0)*(Ephoton-Egeff-Eek-Ehk)/(pow(Ephoton-Egeff-Eek-Ehk,2.0)+pow(gammak,2.0));
			chi0[nk].i = - dipk_array[nk]*(fFDE[nk]+fFDH[nk]-1.0)*gammak/(pow(Ephoton-Egeff-Eek-Ehk,2.0)+pow(gammak,2.0));
		}		

		doublecomplex ** A, * X, * b;
		A = new doublecomplex * [Nk];
		for (int nk = 0; nk < Nk; nk++)
			A[nk] = new doublecomplex [Nk];
		X = new doublecomplex [Nk];
		b = new doublecomplex [Nk];
		for (int nk = 0; nk < Nk; nk++)
			for (int nkp = 0; nkp < Nk; nkp++)
			{
				A[nk][nkp].r = -1.0/pow(2.0*PI,2.0)*k_array[nkp]*dk*Vskkp[nk][nkp]/dipk_array[nk]*chi0[nk].r;
				A[nk][nkp].i = -1.0/pow(2.0*PI,2.0)*k_array[nkp]*dk*Vskkp[nk][nkp]/dipk_array[nk]*chi0[nk].i;
				if (nkp == 0 || nkp == Nk-1)
				{
					A[nk][nkp].r /= 2.0;
					A[nk][nkp].i /= 2.0;
				}
				if (nkp == nk)
				{
					A[nk][nkp].r += 1.0;
				}
			}
		for (int nk = 0; nk < Nk; nk++)
		{
			b[nk].r = chi0[nk].r;
			b[nk].i = chi0[nk].i;
		}
		ZSolveLinEq(A, b, Nk, X);

		for (int nk = 0; nk < Nk; nk++)
			chi[nk] = X[nk];
		
		for (int nk = 0; nk < Nk; nk++)
		{
			//cout << chi0[nk].i << "\t" << chi[nk].i << endl;
		}

		for (int nk = 0; nk < Nk; nk++)
			delete[] A[nk];
		delete [] A;
		delete[] X;
		delete[] b;

		doublecomplex chiw;
		doublecomplex intchiw;
		intchiw.r = 0;
		intchiw.i = 0;
		for (int nk = 0; nk < Nk; nk++)
		{
			doublecomplex addint;
			addint.r = dipk_array[nk]*k_array[nk]*dk*chi[nk].r;
			addint.i = dipk_array[nk]*k_array[nk]*dk*chi[nk].i;
			if (nk == 0 || nk == Nk-1)
			{
				addint.r *= 0.5;
				addint.i *= 0.5;
			}
			intchiw.r += addint.r;
			intchiw.i += addint.i;
		}
		chiw.r = 2.0*intchiw.r/(2.0*PI);  // the first 2.0 is for spin degeneracy
		chiw.i = 2.0*intchiw.i/(2.0*PI);  // the first 2.0 is for spin degeneracy

		doublecomplex chiw0;
		doublecomplex intchiw0;
		intchiw0.r = 0;
		intchiw0.i = 0;
		for (int nk = 0; nk < Nk; nk++)
		{
			doublecomplex addint;
			addint.r = dipk_array[nk]*k_array[nk]*dk*chi0[nk].r;
			addint.i = dipk_array[nk]*k_array[nk]*dk*chi0[nk].i;
			if (nk == 0 || nk == Nk-1)
			{
				addint.r *= 0.5;
				addint.i *= 0.5;
			}
			intchiw0.r += addint.r;
			intchiw0.i += addint.i;
		}
		chiw0.r = 2.0*intchiw0.r/(2.0*PI);  // the first 2.0 is for spin degeneracy
		chiw0.i = 2.0*intchiw0.i/(2.0*PI);  // the first 2.0 is for spin degeneracy
		
		delete [] chi;
		delete[] chi0;
		delete[] gammak_array;
		
		double gainw, gainw0;
		gainw = - 4.0*PI*(Ephoton/hbar)/(sqrt(Well.dielectricinf)*lightspeed)*chiw.i;
		GainPH[np] = gainw;
		
		gainw0 = - 4.0*PI*(Ephoton/hbar)/(sqrt(Well.dielectricinf)*lightspeed)*chiw0.i;
		GainPH0[np] = gainw0;
	}
	

	for (int nk = 0; nk < Nk; nk++)
		delete [] Vskkp[nk];
	delete [] Vskkp;
	delete[] fFDE;
	delete[] fFDH;
	delete[] k_array;

/*	double maxgain = 0;
	int nmaxgain = 0;
	for (int np = 0; np < Nph; np++)
		if (GainPH[np]>maxgain)
		{
			maxgain = GainPH[np];
			nmaxgain = np;
		}

	cout << maxgain << "\t" << nmaxgain << endl;
*/

	const int Nsavef = 200;
	char fgain_name[Nsavef];
	char fgain0_name[Nsavef];
	char prefix[100];
	sprintf_s(prefix, 100, "%s_n2D%3.2e_Tl%.0f_Te%.0f_gamma%.1f_kmax%.1f_",StructureName,n2D*1e-4,T_lattice,T_electron,gammaBR/echarge*1000,kmax/kF);
	sprintf_s(fgain_name, Nsavef, "%sGainPH.txt", prefix);
	sprintf_s(fgain0_name, Nsavef, "%sGainPH0.txt", prefix);

	cout << GainPH[Nph/2] << "\t" << GainPH0[Nph/2] << endl;
	//getchar();

	ofstream gain_out(fgain_name);
	if (!gain_out.is_open())
	{
		cerr << "Cannot open fgain file for output" << endl;
		exit(1);
	}
	for (int kph = 0; kph < Nph; kph++)
		gain_out << Eph[kph]/echarge << "\t";
	gain_out << endl << endl;
	for (int kph = 0; kph < Nph; kph++)
		gain_out << GainPH[kph] << "\t";
	gain_out.close();

	ofstream gain0_out(fgain0_name);
	if (!gain0_out.is_open())
	{
		cerr << "Cannot open fgain0 file for output" << endl;
		exit(1);
	}
	for (int kph = 0; kph < Nph; kph++)
		gain0_out << Eph[kph]/echarge << "\t";
	gain0_out << endl << endl;
	for (int kph = 0; kph < Nph; kph++)
		gain0_out << GainPH0[kph] << "\t";
	gain0_out.close();


	delete [] dipk_array;

	delete [] q_array;
	delete[] Vsq_array;
	delete[] kk_array;
	delete[] Eekk_array;
	delete[] Ehkk_array;
	delete[] Eexkk_array;
	delete[] EeRkk_array;
	delete[] EhRkk_array;

	}
	
	delete[] qform_array;
	delete[] Vqform_array;
	delete [] z_grid;
	delete[] Ve_grid;
	delete[] Ve_grid_halfdzshift;
	delete[] Vhh_grid;
	delete[] Vhh_grid_halfdzshift;

	delete [] layer;
	delete[] isbarrier;
	delete[] border;

	

	for (int n = 0; n < nEe; n++)
		delete [] eigVe[n];
	delete [] eigEe;
	delete[] eigVe;
	delete[] effmxye_ave;

	for (int n = 0; n < nEhh; n++)
		delete [] eigVhh[n];
	delete [] eigEhh;
	delete[] eigEhhv;
	delete[] eigVhh;
	delete[] effmxyhh_ave;
	
	
	return 0;
}

