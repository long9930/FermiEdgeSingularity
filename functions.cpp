// functions.cpp, put all functions in this file

#include "stdafx.h"
#include "functions.h"


// convert bandparameters into SI units
void BandParaSIunits(BandParameters A, BandParameters &A_SI)
{
	A_SI.Eg = A.Eg * echarge;
	A_SI.Ep = A.Ep * echarge;
	A_SI.VBO = A.VBO * echarge;
	A_SI.KaneF = A.KaneF;
	A_SI.DeltaSO = A.DeltaSO * echarge;
	A_SI.dielectric0 = A.dielectric0;
	A_SI.dielectricinf = A.dielectricinf;
	A_SI.E_LO = A.E_LO * echarge;
	A_SI.lc = A.lc * 1.0e-10;
	A_SI.ac = A.ac * echarge;
	A_SI.av = A.av * echarge;
	A_SI.b = A.b * echarge;
	A_SI.d = A.d * echarge;
	A_SI.C11 = A.C11 * 1.0e9;
	A_SI.C12 = A.C12 * 1.0e9;
	A_SI.C44 = A.C44 * 1.0e9;
}


void SolveLinEq(double ** A, double *b, int Nd, double *X)
{
	// solve A*X = b, A is Nd by Nd matrix, b and X are column vectors with Nd elements.
	doublereal * A_;
	A_ = new doublereal[Nd*Nd];
	int kk = 0;
	for (int col = 0; col < Nd; col++)
		for (int row = 0; row < Nd; row++)
		{
			A_[kk] = A[row][col];
			kk++;
		}
	integer N_ = Nd;
	integer nrhs_ = 1;
	integer lda_ = Nd;
	integer * ipiv_;
	ipiv_ = new integer [N_];
	doublereal *b_;
	b_ = new doublereal [N_];
	for (int k = 0; k < N_; k++)
		b_[k] = b[k];
	integer ldb_ = N_;
	integer info_;
	dgesv_(&N_, &nrhs_, A_, &lda_, ipiv_, b_, &ldb_, &info_);
	if(info_ != 0) // failed 
	{
		printf("Error in solving linear equations.\n");
		fprintf(stderr, "dsyev_ fails %d\n", info_);
		exit(1);
	}
	else
	{
		for (int k = 0; k < Nd; k++)
			X[k] = b_[k];
	}
	delete [] A_;
	delete [] b_;
	delete [] ipiv_;
}

void ZSolveLinEq(doublecomplex ** A, doublecomplex *b, int Nd, doublecomplex *X)
{
	// solve A*X = b, A is Nd by Nd matrix, b and X are column vectors with Nd elements.
	doublecomplex * A_;
	A_ = new doublecomplex[Nd*Nd];
	int kk = 0;
	for (int col = 0; col < Nd; col++)
		for (int row = 0; row < Nd; row++)
		{
			A_[kk] = A[row][col];
			kk++;
		}
	integer N_ = Nd;
	integer nrhs_ = 1;
	integer lda_ = Nd;
	integer * ipiv_;
	ipiv_ = new integer [N_];
	doublecomplex *b_;
	b_ = new doublecomplex [N_];
	for (int k = 0; k < N_; k++)
		b_[k] = b[k];
	integer ldb_ = N_;
	integer info_;
	zgesv_(&N_, &nrhs_, A_, &lda_, ipiv_, b_, &ldb_, &info_);
	if(info_ != 0) // failed 
	{
		printf("Error in solving linear equations.\n");
		fprintf(stderr, "dsyev_ fails %d\n", info_);
		exit(1);
	}
	else
	{
		for (int k = 0; k < Nd; k++)
			X[k] = b_[k];
	}

	delete []A_;
	delete [] b_;
	delete [] ipiv_;
}


double my_interp_eqdx(double x[], double y[], int N, double x0)
{
	// linear interpolation
	using namespace std;
	if (x0 < x[0] || x0 > x[N-1])
	{
		cout << "The input value is out of range of the known function1!" << endl;
		exit(1);
	}
	double dx, y0;
	dx = x[1]-x[0];
	int n0 = static_cast<int>(floor((x0-x[0])/dx));
	if (n0 == N-1)
		y0 = y[N-1];
	else
		y0 = y[n0] + (y[n0+1]-y[n0])/dx*(x0-x[n0]);
	return y0;
}

double my_interp_neqdx(double x[], double y[], int N, double x0)
{
	// linear interpolation
	using namespace std;
	if (x0 < x[0] || x0 >= x[N-1])
	{
		cout << x[0] << "\t" << x[N-1] << "\t" << x0 << endl;
		cout << "The input value is out of range of the known function2!" << endl;
		exit(1);
	}
	int n = 0;
	while (x0 >= x[n])
		n++;	
	double dx, y0;
	dx = x[n]-x[n-1];
	y0 = y[n-1] + (y[n]-y[n-1])/dx*(x0-x[n-1]);
	return y0;
}


int findlayer(double z, double * border, int Nlayer)
{
	using namespace std;
	if (z < 0 || z > border[Nlayer-1])
	{
		cout << "The point z is not in the considered region! \n";
		exit(1);
	}
	int klayer = -1;
	for (int kb = 0; kb < Nlayer; kb++)
	{
		if (z < border[kb])
//		if (z <= border[kb])
		{
			klayer = kb;
			break;
		}		
	}
	if (z == border[Nlayer-1])
		klayer = Nlayer-1;
	if (klayer < 0)
	{
		cout << "Error in defining effmz!" << endl;
		exit(1);
	}
	return klayer;
}


void GetNormalized(double WaveFunc[], int Nsize)
{
	double psi_sq_sum = 0, psi_norm;
	for (int kg = 0; kg < Nsize; kg++)
		psi_sq_sum += pow(WaveFunc[kg],2);
	psi_norm = sqrt(psi_sq_sum);
	for (int kg = 0; kg < Nsize; kg++)
		WaveFunc[kg] = WaveFunc[kg]/psi_norm;
}


double factorial (int num)
{
	using namespace std;
	if (num < 0)
	{
		cerr << "no negative number when calculating factorials!" << endl;
		exit(1);
	}
	double result = 1;
	for (int i=1; i<=num; ++i)
		result *= static_cast<double>(i);
	return result;
}


double HermiteH(int n, double x)
{
	int halfn = static_cast<int>(floor(static_cast<double>(n)/2.0));
	double y = 0;
	for (int m = 0; m <= halfn; m++)
	{
		double addy;
		addy = 1.0/(factorial(m)*factorial(n-2*m))*pow(2.0*x,static_cast<double>(n-2*m));
		if (m%2 != 0)
			addy *= -1.0;
		y += addy;
	}
	y *= factorial(n);
	return y;
}


double phiHarmonic(int n, double aH, double x)
{
	double norm, y;
	norm = 1.0/sqrt(sqrt(PI)*aH*pow(2.0,static_cast<double>(n))*factorial(n));
	y = exp(-x*x/(2.0*aH*aH))*HermiteH(n,x/aH);
	y = y*norm;
	return y;
}


void fsaveVqform(char fVqform_name[], double * qform_array, double * Vqform_array, int Nqform)
{
	using namespace std;
	ofstream fout(fVqform_name);
	if (!fout.is_open())
	{
		cerr << "Cannot open fVqform file for output!" << endl;
		exit(1);
	}
	for (int nq = 0; nq < Nqform; nq++)
		fout << qform_array[nq] << "\t";
	fout << "\n\n";
	for (int nq = 0; nq < Nqform; nq++)
		fout << Vqform_array[nq] << "\t";
	fout << "\n\n";		
	fout.close();
}


void freadVqform(char fVqform_name[], double * &qform_array, double * &Vqform_array, int &Nqform)
{
	using namespace std;
	ifstream fin(fVqform_name);
	if (!fin.is_open())
	{
		cerr << "Cannot open fVqform file for input!" << endl;
		exit(1);
	}
	double * readtemp = new double [20000];
	int kt = 0;
	while (fin.good())
	{
		fin >> readtemp[kt];
		//cout << kt << "\t" << readtemp[kt] << endl;
		kt++;
	}
	kt--;
	//cout << "kt\t" << kt << endl;
	if (kt%2 != 0)
	{
		cerr << "error in reading Vqform!" << endl;
		exit(1);
	}
	Nqform = kt/2;
	qform_array = new double [Nqform];
	Vqform_array = new double [Nqform];
	for (int nq = 0; nq < Nqform; nq++)
		qform_array[nq] = readtemp[nq];
	for (int nq = 0; nq < Nqform; nq++)
		Vqform_array[nq] = readtemp[nq+Nqform];
	fin.close();
	delete [] readtemp;
}




void SetBandParaGaAs(BandParameters &GaAs, double T_lattice)
// band parameters for GaAs
{
	double Varshnia, Varshnib;
	Varshnia = 0.5405e-3;
	Varshnib = 204;
	GaAs.Eg = 1.519 - Varshnia * pow(T_lattice,2) / (T_lattice + Varshnib);
	GaAs.Ep = 28.8;
	GaAs.VBO = -0.8;
	GaAs.KaneF = -1.94;
	GaAs.DeltaSO = 0.341;
	GaAs.dielectric0 = 12.9;
	GaAs.dielectricinf = 10.89;
	GaAs.E_LO = 35e-3;
	GaAs.gamma1 = 6.98;
	GaAs.gamma2 = 2.06;
	GaAs.gamma3 = 2.93;
	GaAs.lc = 5.65325 + 3.88e-5*(T_lattice-300);
	GaAs.ac = -7.17;
	GaAs.av = -1.16;
	GaAs.b = -2.0;
	GaAs.d = -4.8;
	GaAs.C11 = 1221;
	GaAs.C12 = 566;
	GaAs.C44 = 600;
}


void SetBandParaInAs(BandParameters &InAs, double T_lattice)
// band parameters for InAs
{
	double Varshnia, Varshnib;
	Varshnia = 0.276e-3;
	Varshnib = 93;
	InAs.Eg = 0.417 - Varshnia * pow(T_lattice,2) / (T_lattice + Varshnib);
	InAs.Ep = 21.5;
	InAs.VBO = -0.59;
	InAs.KaneF = -2.90;
	InAs.DeltaSO = 0.39;
	InAs.dielectric0 = 15.15;
	InAs.dielectricinf = 12.3;
	InAs.E_LO = 30e-3;
	InAs.gamma1 = 20.0;
	InAs.gamma2 = 8.5;
	InAs.gamma3 = 9.2;
	InAs.lc = 6.0583 + 2.74e-5*(T_lattice-300);
	InAs.ac = -5.08;
	InAs.av = -1.00;
	InAs.b = -1.8;
	InAs.d = -3.6;
	InAs.C11 = 832.9;
	InAs.C12 = 452.6;
	InAs.C44 = 395.9;
}


void SetBandParaGaInAs(BandParameters &GaInAs, BandParameters GaAs, BandParameters InAs, double T_lattice, double x_GaInAs)
// band parameters for Ga_(1-x)In_(x)As
{
	double x = x_GaInAs;
	GaInAs.Eg = (1-x)*GaAs.Eg + x*InAs.Eg - x*(1-x)*0.477;
	GaInAs.Ep = (1-x)*GaAs.Ep + x*InAs.Ep - x*(1-x)*(-1.48);
	GaInAs.VBO = (1-x)*GaAs.VBO + x*InAs.VBO - x*(1-x)*(-0.38);
	GaInAs.KaneF = (1-x)*GaAs.KaneF + x*InAs.KaneF - x*(1-x)*(1.77);
	GaInAs.DeltaSO = (1-x)*GaAs.DeltaSO + x*InAs.DeltaSO - x*(1-x)*(0.15);
	GaInAs.dielectric0 = 15.1 - 2.87*(1-x) + 0.67*pow(1-x,2);
	GaInAs.dielectricinf = 12.3 - 1.4*(1-x);
	GaInAs.E_LO = (1-x)*GaAs.E_LO + x*InAs.E_LO;
	GaInAs.gamma1 = (1-x)*GaAs.gamma1 + x*InAs.gamma1;
	GaInAs.gamma2 = (1-x)*GaAs.gamma2 + x*InAs.gamma2;
	GaInAs.gamma3 = (1-x)*GaAs.gamma3 + x*InAs.gamma3;
	GaInAs.lc = (1-x)*GaAs.lc + x*InAs.lc;
	GaInAs.ac = (1-x)*GaAs.ac + x*InAs.ac - x*(1-x)*2.61;
	GaInAs.av = (1-x)*GaAs.av + x*InAs.av;
	GaInAs.b = (1-x)*GaAs.b + x*InAs.b;
	GaInAs.d = (1-x)*GaAs.d + x*InAs.d;
	GaInAs.C11 = (1-x)*GaAs.C11 + x*InAs.C11;
	GaInAs.C12 = (1-x)*GaAs.C12 + x*InAs.C12;
	GaInAs.C44 = (1-x)*GaAs.C44 + x*InAs.C44;
}


void GetStrain(BandParameters Barrier, BandParameters Well, BandParameters Substrate, Strain &strB, Strain &strW)
{
	strW.exx = (Substrate.lc-Well.lc)/Well.lc;
	strW.eyy = strW.exx;
	strW.ezz = -(2.0*Well.C12/Well.C11)*strW.exx;
	strW.cPe = Well.ac*(strW.exx+strW.eyy+strW.ezz);
	strW.vPe = Well.av*(strW.exx+strW.eyy+strW.ezz);
	strW.Qe = Well.b/2.0*(2*strW.ezz-strW.exx-strW.eyy);
	strW.Ehh = -strW.vPe-strW.Qe;
	strW.Elh = -strW.vPe+1.0/2.0*(strW.Qe-Well.DeltaSO+sqrt(pow((Well.DeltaSO),2)+2*Well.DeltaSO*strW.Qe+9*pow(strW.Qe,2)));
	strW.Eso = -strW.vPe+1.0/2.0*(strW.Qe-Well.DeltaSO-sqrt(pow((Well.DeltaSO),2)+2*Well.DeltaSO*strW.Qe+9*pow(strW.Qe,2)));
	if (fabs(strW.exx) < 1e-5)
	{
		strW.alphamix = 1;
		strW.betamix = 0;
	}
	else
	{
		strW.Amix = Well.DeltaSO+strW.Qe;
		strW.Bmix = sqrt(pow(Well.DeltaSO,2)+2*Well.DeltaSO*strW.Qe+9*pow(strW.Qe,2));
		strW.Cmix = sqrt(2*strW.Bmix*(strW.Bmix-strW.Amix));
		strW.alphamix = 2*sqrt(2.0)*fabs(strW.Qe)/strW.Cmix;
		strW.betamix = (strW.Amix-strW.Bmix)*fabs(strW.Qe)/(strW.Cmix*strW.Qe);
	}

	strB.exx = (Substrate.lc-Barrier.lc)/Barrier.lc;
	strB.eyy = strB.exx;
	strB.ezz = -(2.0*Barrier.C12/Barrier.C11)*strB.exx;
	strB.cPe = Barrier.ac*(strB.exx+strB.eyy+strB.ezz);
	strB.vPe = Barrier.av*(strB.exx+strB.eyy+strB.ezz);
	strB.Qe = Barrier.b/2.0*(2*strB.ezz-strB.exx-strB.eyy);
	strB.Ehh = -strB.vPe-strB.Qe;
	strB.Elh = -strB.vPe+1.0/2.0*(strB.Qe-Barrier.DeltaSO+sqrt(pow((Barrier.DeltaSO),2)+2*Barrier.DeltaSO*strB.Qe+9*pow(strB.Qe,2)));
	strB.Eso = -strB.vPe+1.0/2.0*(strB.Qe-Barrier.DeltaSO-sqrt(pow((Barrier.DeltaSO),2)+2*Barrier.DeltaSO*strB.Qe+9*pow(strB.Qe,2)));
	if (fabs(strB.exx) < 1e-5)
	{
		strB.alphamix = 1;
		strB.betamix = 0;
	}
	else
	{
		strB.Amix = Barrier.DeltaSO+strB.Qe;
		strB.Bmix = sqrt(pow(Barrier.DeltaSO,2)+2*Barrier.DeltaSO*strB.Qe+9*pow(strB.Qe,2));
		strB.Cmix = sqrt(2*strB.Bmix*(strB.Bmix-strB.Amix));
		strB.alphamix = 2*sqrt(2.0)*fabs(strB.Qe)/strB.Cmix;
		strB.betamix = (strB.Amix-strB.Bmix)*fabs(strB.Qe)/(strB.Cmix*strB.Qe);
	}
}


double effmzE(BandParameters Material, Strain strM, double DelE)
{
	double KaneF, Ep, Eg, alphamix, betamix, cPe, Elh, Eso;
	KaneF = Material.KaneF;
	Ep = Material.Ep;
	Eg = Material.Eg;
	alphamix = strM.alphamix;
	betamix = strM.betamix;
	cPe = strM.cPe;
	Elh = strM.Elh;
	Eso = strM.Eso;
	double effmz_inverse;
	effmz_inverse =  1+2*KaneF+Ep/3*(pow((sqrt(2.0)*alphamix-betamix),2)/(Eg+cPe+DelE-Elh) + pow((sqrt(2.0)*betamix+alphamix),2)/(Eg+cPe+DelE-Eso));
	return m0/effmz_inverse;
}


double effmxyE(BandParameters Material, Strain strM, double DelE)
{
	double KaneF, Ep, Eg, alphamix, betamix, cPe, Elh, Ehh, Eso;
	KaneF = Material.KaneF;
	Ep = Material.Ep;
	Eg = Material.Eg;
	alphamix = strM.alphamix;
	betamix = strM.betamix;
	cPe = strM.cPe;
	Elh = strM.Elh;
	Ehh = strM.Ehh;
	Eso = strM.Eso;
	double effmxy_inverse;
	effmxy_inverse = 1+2*KaneF+Ep/6*( 3/(Eg+cPe+DelE-Ehh) + pow((alphamix-sqrt(2.0)*betamix),2.0)/(Eg+cPe+DelE-Elh)
		+ pow((sqrt(2.0)*alphamix+betamix),2.0)/(Eg+cPe+DelE-Eso));
	return m0/effmxy_inverse;
}


double effmzHH(BandParameters Material, Strain strM, double DelE)
{
	return m0 / (Material.gamma1 - 2*Material.gamma2);
}


double effmxyHH(BandParameters Material, Strain strM, double DelE)
{
	return m0 / (Material.gamma1 + Material.gamma2);
}


double effmzLH(BandParameters Material, Strain strM, double DelE)
{
	return m0 / (Material.gamma1 + 2*Material.gamma2);
}


double effmxyLH(BandParameters Material, Strain strM, double DelE)
{
	return m0 / (Material.gamma1 - Material.gamma2);
}


void getBandStruct(BandParameters Barrier, BandParameters Well, Strain strB, Strain strW,
	double * layer, bool * isbarrier, double * border, int Nlayer, double * z_grid, double * V_grid, double * z_grid_halfdzshift, double * V_grid_halfdzshift, int N_grid, 
	double (*effmzEdep)(BandParameters, Strain, double delE), double (*effmxyEdep)(BandParameters, Strain, double delE),
	int N_level, int &nEbound, double * eigE, double ** eigV, double * effmxy_ave)
{
	using namespace std;

	const double dE = 0.1e-3*echarge;  // step for increasing energy, it's valid unless two levels residing in this range.

	double dz = z_grid[1] - z_grid[0];

	double maxV, minV;
	maxV = V_grid[0];
	minV = V_grid[0];
	for (int kg = 0; kg < N_grid; kg++)
	{
		if (V_grid[kg] < minV)
			minV = V_grid[kg];
		if (V_grid[kg] > maxV)
			maxV = V_grid[kg];
	}
	
	int nE = 0;
	while (nE < N_level)
	{
		cout << nE << endl;
		double E1, E2;
		double * psi1 = new double[N_grid], * psi2 = new double[N_grid];
		double * effmz1 = new double[N_grid], * effmz2 = new double[N_grid];
		double * effmz_halfdzshift1 = new double[N_grid], * effmz_halfdzshift2 = new double[N_grid];
		double psi_rboundary1, psi_rboundary2;
		double E_temp, psi_rboundary_temp;
		double * psi_temp = new double[N_grid], * effmz_temp = new double[N_grid], * effmz_halfdzshift_temp = new double[N_grid];
		double psi_rboundary_limit = 1e-10;
		double convergence_limitE = 1e-15*maxV;
		double psi_i0 = 0, psi_i1 = 1;
		
		if (nE == 0)
			E1 = minV + 0.01*dE;
		else
			E1 = eigE[nE-1] + 0.01*dE;
		
		// find the right boundary value of the wave function for the initial E1
		E1 = E1 - dE;
		bool negaeffm = 0;
		do
		{
		E1 += dE;
			negaeffm = 0;
			for (int kg = 0; kg < N_grid; kg++)
			{
				int klayer;
				klayer = findlayer(z_grid[kg], border, Nlayer);
				if (isbarrier[klayer])
					effmz1[kg] = effmzEdep(Barrier, strB, E1-V_grid[kg]);
				else
					effmz1[kg] = effmzEdep(Well, strW, E1-V_grid[kg]);
				if (effmz1[kg] < 0)
				{
					negaeffm = 1;
					break;
				}
				int klayer_halfdzshift;
				klayer_halfdzshift = findlayer(z_grid_halfdzshift[kg], border, Nlayer);
				if (isbarrier[klayer_halfdzshift])
					effmz_halfdzshift1[kg] = effmzEdep(Barrier, strB, E1-V_grid_halfdzshift[kg]);
				else
					effmz_halfdzshift1[kg] = effmzEdep(Well, strW, E1-V_grid_halfdzshift[kg]);
			}
			if (negaeffm == 0)
			{
				psi1[0] = psi_i0;
				psi1[1] = psi_i1;
				for (int kg = 1; kg < N_grid-1; kg++)
				{
					double m1, m2, Vz;
					m1 = effmz_halfdzshift1[kg-1];
					m2 = effmz_halfdzshift1[kg];
					// m1 = (effmz1[kg-1] + effmz1[kg])/2;  // used in Harrison's book.
					// m2 = (effmz1[kg+1] + effmz1[kg])/2;
					Vz = V_grid[kg];
					psi1[kg+1] = m2 * ( (2*pow(dz/hbar,2)*(Vz-E1) + 1/m1 + 1/m2)*psi1[kg] - psi1[kg-1]/m1 );
				}
				GetNormalized(psi1, N_grid);
				psi_rboundary1 = psi1[N_grid-1];
			}
		} while (negaeffm == 1);
				
		// try to locate the position of the solution
		while(1)
		{			
			E2 = E1 + dE;
			for (int kg = 0; kg < N_grid; kg++)
			{
				int klayer;
				klayer = findlayer(z_grid[kg], border, Nlayer);
				if (isbarrier[klayer])
					effmz2[kg] = effmzEdep(Barrier, strB, E2-V_grid[kg]);
				else
					effmz2[kg] = effmzEdep(Well, strW, E2-V_grid[kg]);
				int klayer_halfdzshift;
				klayer_halfdzshift = findlayer(z_grid_halfdzshift[kg], border, Nlayer);
				if (isbarrier[klayer_halfdzshift])
					effmz_halfdzshift2[kg] = effmzEdep(Barrier, strB, E2-V_grid_halfdzshift[kg]);
				else
					effmz_halfdzshift2[kg] = effmzEdep(Well, strW, E2-V_grid_halfdzshift[kg]);
			}
			psi2[0] = psi_i0;
			psi2[1] = psi_i1;
			for (int kg = 1; kg < N_grid-1; kg++)
			{
				double m1, m2, Vz;
				m1 = effmz_halfdzshift2[kg-1];
				m2 = effmz_halfdzshift2[kg];				
				Vz = V_grid[kg];
				psi2[kg+1] = m2 * ( (2*pow(dz/hbar,2)*(Vz-E2) + 1/m1 + 1/m2)*psi2[kg] - psi2[kg-1]/m1 );
			}
			GetNormalized(psi2, N_grid);
			psi_rboundary2 = psi2[N_grid-1];
			if (psi_rboundary1*psi_rboundary2 < 0)
				break;
			else
			{
				E1 = E2;
				psi_rboundary1 = psi_rboundary2;
			}				
		} 	

		// try to find the eigenvalue 
		while (fabs(E2-E1) > convergence_limitE )
		{
			E_temp = (E1 + E2)/2;
			for (int kg = 0; kg < N_grid; kg++)
			{
				int klayer;
				klayer = findlayer(z_grid[kg], border, Nlayer);
				if (isbarrier[klayer])
					effmz_temp[kg] = effmzEdep(Barrier, strB, E_temp-V_grid[kg]);						
				else
					effmz_temp[kg] = effmzEdep(Well, strW, E_temp-V_grid[kg]);
				int klayer_halfdzshift;
				klayer_halfdzshift = findlayer(z_grid_halfdzshift[kg], border, Nlayer);
				if (isbarrier[klayer_halfdzshift])
					effmz_halfdzshift_temp[kg] = effmzEdep(Barrier, strB, E_temp-V_grid_halfdzshift[kg]);
				else
					effmz_halfdzshift_temp[kg] = effmzEdep(Well, strW, E_temp-V_grid_halfdzshift[kg]);
			}
			psi_temp[0] = psi_i0;
			psi_temp[1] = psi_i1;
			for (int kg = 1; kg < N_grid-1; kg++)
			{
				double m1, m2, Vz;
				m1 = effmz_halfdzshift_temp[kg-1];
				m2 = effmz_halfdzshift_temp[kg];
				Vz = V_grid[kg];
				psi_temp[kg+1] = m2 * ( (2*pow(dz/hbar,2)*(Vz-E_temp) + 1/m1 + 1/m2)*psi_temp[kg] - psi_temp[kg-1]/m1 );
			}
			GetNormalized(psi_temp, N_grid);
			psi_rboundary_temp = psi_temp[N_grid-1];
					
			if (psi_rboundary_temp*psi_rboundary1 > 0)
			{
				E1 = E_temp;
				psi_rboundary1 = psi_rboundary_temp;
			}
			else
			{
				E2 = E_temp;
				psi_rboundary2 = psi_rboundary_temp;
			}
		}

		if (E1 > maxV)
		{
			cout << "eigen energy has gone beyond bounded limit. finished!" << endl;
			break;
		}

		eigE[nE] = E_temp;

		double effmxy_ave_inverse = 0;
		for (int kg = 0; kg < N_grid; kg++)
		{
			int klayer;
			double add_effmxy_ave_inverse;
			klayer = findlayer(z_grid[kg], border, Nlayer);
			if (isbarrier[klayer])
				add_effmxy_ave_inverse = pow(psi_temp[kg],2.0)/effmxyEdep(Barrier, strB, E_temp-V_grid[kg]);
			else
				add_effmxy_ave_inverse = pow(psi_temp[kg],2.0)/effmxyEdep(Well, strW, E_temp-V_grid[kg]);
			if (add_effmxy_ave_inverse < 0)
			{
				cerr << "Negative in-plane effective mass occured!" << endl;
				exit(1);
			}
			effmxy_ave_inverse += add_effmxy_ave_inverse;
		}
		effmxy_ave[nE] = 1.0/effmxy_ave_inverse;
		
		eigV[nE] = new double [N_grid];

		for (int kg = 0; kg < N_grid; kg++)
			eigV[nE][kg] = psi_temp[kg];

		nE++;
		delete[] psi1;
		delete[] psi2;
		delete[] effmz1;
		delete[] effmz2;
		delete[] effmz_halfdzshift1;
		delete[] effmz_halfdzshift2;
		delete[] psi_temp;
		delete[] effmz_temp;
		delete[] effmz_halfdzshift_temp;
	}
	
	cout << "The eigen energies are (meV):" << endl;
	for (int kE = 0; kE < nE; kE++)
		cout << eigE[kE]/echarge*1000 << endl;

	nEbound = nE;
}

	
	
double CoulombHolePP(double n2D, double epsr0, double effmR, double kappa, double Cq)
{
	using namespace std;
	double qmax;
	qmax = pow(100.0*32.0*PI*effmR*eGauss*eGauss*n2D/(Cq*hbar*hbar*epsr0),1.0/3.0);
	cout << qmax << "\t" << qmax/kappa << "\t" << (Cq/4.0)*pow(hbar/(2*effmR),2.0)/(2*PI*pow(eGauss,2.0)*n2D/(epsr0*effmR))*pow(qmax,3.0) << endl;
	//cout << pow(1000/((Cq/4.0)*pow(hbar/(2*effmR),2.0)/(2*PI*pow(eGauss,2.0)*n2D/(epsr0*effmR))),1.0/3.0) << endl;
	int Nq = 5000;
	double dq = qmax/static_cast<double>(Nq-1);
	double intq = 0;
	for (int nq = 0; nq < Nq; nq++)
	{
		double qq = static_cast<double>(nq)*dq;
		double addintq;
		addintq = - 1.0/(1.0 + qq/kappa + (Cq/4.0)*pow(hbar/(2*effmR),2.0)/(2*PI*pow(eGauss,2.0)*n2D/(epsr0*effmR))*pow(qq,3.0));
		if (nq==0 || nq==Nq-1)
			addintq *= 0.5;
		intq += addintq*dq;
	}
	return intq*eGauss*eGauss/epsr0;
}


double n2DEF(double * k_array, double * Ek_array, int Nk, double mu, double T)
{
	double kmax;
	if (mu+10.0*kB*T > Ek_array[Nk-1])
	{
		std::cerr << "Fermi level plus 10*kB*T is out of range!\n";
		exit(1);
	}
	kmax = my_interp_neqdx(Ek_array, k_array, Nk, mu+10.0*kB*T);
	int Nkint = 2000;
	double dk = kmax/static_cast<double>(Nkint-1);
	double intk = 0;
	for (int nk = 0; nk < Nkint; nk++)
	{
		double kk, Ekk, fFDkk, addintk;
		kk = static_cast<double>(nk)*dk;
		Ekk = my_interp_eqdx(k_array, Ek_array, Nk, kk);
		fFDkk = 1.0/(exp((Ekk-mu)/(kB*T))+1.0);
		addintk = kk*fFDkk;
		if (nk == 0 || nk == Nkint-1)
			addintk *= 0.5;
		intk += addintk*dk;
	}
	return intk/PI;
}


double n2DEF_dxdy(double * k_array, double * Ek_array, int Nk, double mu, double T)
{
	double KMAX;
	if (mu+10.0*kB*T > Ek_array[Nk-1])
	{
		std::cerr << "Fermi level plus 10*kB*T is out of range!\n";
		exit(1);
	}
	KMAX = my_interp_neqdx(Ek_array, k_array, Nk, mu+10.0*kB*T);
	int Nkint = 1000;
	double dk = KMAX/static_cast<double>(Nkint-1);
	double kxmin, kxmax, kymin, kymax;
	kxmin = - KMAX;
	kxmax = KMAX;
	kymin = - KMAX;
	kymax = KMAX;
	double dkx, dky;
	dkx = (kxmax-kxmin)/static_cast<double>(Nkint-1);
	dky = (kymax-kymin)/static_cast<double>(Nkint-1);
	double intkx = 0;
	for (int nkx = 0; nkx < Nkint; nkx++)
	{
		double kx = kxmin + static_cast<double>(nkx)*dkx;
		double intky = 0;
		for (int nky = 0; nky < Nkint; nky++)
		{
			double ky = kymin + static_cast<double>(nky)*dky;
			double kk = sqrt(kx*kx + ky*ky);
			double Ekk, fFDkk;
			Ekk = my_interp_eqdx(k_array, Ek_array, Nk, kk);
			fFDkk = 1.0/(exp((Ekk-mu)/(kB*T))+1.0);
			double addintky;
			addintky = fFDkk;
			if (nky == 0 || nky == Nkint-1)
				addintky *= 0.5;
			intky += addintky*dky;
		}
		double addintkx;
		addintkx = intky;
		if (nkx == 0 || nkx == Nkint-1)
			addintkx *= 0.5;
		intkx += addintkx*dkx;
	}
	return 2.0/pow(2.0*PI,2.0)*intkx;
}


double FermiLevel(double * k_array, double * Ek_array, int Nk, double n2D, double T)
{
	using namespace std;
	double mu1, mu2, mutemp, n2D1, n2D2, n2Dtemp;
	mu1 = Ek_array[0]-5.0*kB*T;
	mu2 = Ek_array[Nk-1]-15.0*kB*T;
	n2D1 = n2DEF(k_array, Ek_array, Nk, mu1, T);
	n2D2 = n2DEF(k_array, Ek_array, Nk, mu2, T);
	if ((n2D1-n2D)*(n2D2-n2D)>0)
	{
		cerr << "Error in finding Fermi level at begining!\n";
		exit(1);
	}
	do
	{
		mutemp = (mu1+mu2)/2.0;
		n2Dtemp = n2DEF(k_array, Ek_array, Nk, mutemp, T);
		if ((n2Dtemp-n2D)*(n2D1-n2D)>0)
		{
			mu1 = mutemp;
			n2D1 = n2Dtemp;
		}
		else
		{
			mu2 = mutemp;
			n2D2 = n2Dtemp;
		}
		if ((n2D1-n2D)*(n2D2-n2D)>0)
		{
			cerr << "Error in finding Fermi level!\n";
			exit(1);
		}
	}while(fabs((mu2-mu1)/mu1) > 1e-5);
	
	return mu1;
}


double Polarizability_normal(double * k_array, double * Ek_array, int Nk, double n2D, double T, double q)
{
	using namespace std;
	double muF;
	muF = FermiLevel(k_array, Ek_array, Nk, n2D, T);
	//cout << "muF\t" << muF << endl;
	double qx, qy;
	qx = q;
	qy = 0;
	double KMAX;
	KMAX = my_interp_neqdx(Ek_array, k_array, Nk, muF+15.0*kB*T) + q;
	//double EKMAX = my_interp_eqdx(k_array, Ek_array, Nk, KMAX);
	//cout << 1.0/(exp((EKMAX-muF)/(kB*T))+1.0) << endl;
	//getchar();
	double kxmin, kxmax, kymin, kymax;
	kxmin = - KMAX;
	kxmax = KMAX;
	kymin = - KMAX;
	kymax = KMAX;
	int Nkint = 1000;
	double dkx, dky;
	dkx = (kxmax-kxmin)/static_cast<double>(Nkint-1);
	dky = (kymax-kymin)/static_cast<double>(Nkint-1);
	double intkx = 0;
	for (int nkx = 0; nkx < Nkint; nkx++)
	{
		double kx = kxmin + static_cast<double>(nkx)*dkx;
		double intky = 0;
		for (int nky = 0; nky < Nkint; nky++)
		{
			double ky = kymin + static_cast<double>(nky)*dky;
			double kk = sqrt(kx*kx + ky*ky);
			double kqx, kqy, kq;
			kqx = kx + qx;
			kqy = ky + qy;
			kq = sqrt(kqx*kqx + kqy*kqy);
			//cout << (kq-kk) << endl;
			//getchar();
			double Ekk, Ekq, fFDkk, fFDkq;
			Ekk = my_interp_eqdx(k_array, Ek_array, Nk, kk);
			Ekq = my_interp_eqdx(k_array, Ek_array, Nk, kq);
			fFDkk = 1.0/(exp((Ekk-muF)/(kB*T))+1.0);
			fFDkq = 1.0/(exp((Ekq-muF)/(kB*T))+1.0);
			double addintky;
			if (fabs(fFDkk-fFDkq) < 1e-6)
			{
				//cout << fFDkk << endl;
				//getchar();
				double f_1;
				f_1 = - 1.0/(exp(-(Ekk-muF)/(kB*T))+1.0);
				addintky = fFDkk*f_1/(kB*T);
				//addintky = fFDkk*(fFDkk-1.0)/(kB*T);
			}
			else
				addintky = (fFDkq-fFDkk)/(Ekq-Ekk);
			if (nky == 0 || nky == Nkint-1)
				addintky *= 0.5;
			intky += addintky*dky;
		}
		double addintkx;
		addintkx = intky;
		if (nkx == 0 || nkx == Nkint-1)
			addintkx *= 0.5;
		intkx += addintkx*dkx;
	}
	double Pol = 2.0*intkx/pow(2.0*PI,2.0);
	return Pol;
}


double Polarizability(double * k_array, double * Ek_array, int Nk, double n2D, double T, double q)
{
	using namespace std;
	if ( q < 0)
	{
		cout << "q cannot be negative in calculating polarizability!\n";
		exit(1);
	}
	double muF;
	muF = FermiLevel(k_array, Ek_array, Nk, n2D, T);
	//cout << "muF\t" << muF << endl;
	double qx, qy;
	qx = q;
	qy = 0;  // qy must be set to be 0 in this function
	double KMAX;
	KMAX = my_interp_neqdx(Ek_array, k_array, Nk, muF+15.0*kB*T);
	double kxmin, kxmax, kymin, kymax;
	kxmin = - KMAX;
	kxmax = KMAX;
	kymin = - KMAX;
	kymax = KMAX;
	double qkxmin, qkxmax, qkymin, qkymax;
	qkxmin = -qx - KMAX;
	qkxmax = min(kxmin, -qx + KMAX);
	qkymin = -qy - KMAX;
	qkymax = -qy + KMAX;
	int Nkint = 400;
	double dkx, dky;
	dkx = (kxmax-kxmin)/static_cast<double>(Nkint-1);
	dky = (kymax-kymin)/static_cast<double>(Nkint-1);
	double dqkx, dqky;
	dqkx = (qkxmax-qkxmin)/static_cast<double>(Nkint-1);
	dqky = (qkymax-qkymin)/static_cast<double>(Nkint-1);
	double intkx = 0;
	for (int nkx = 0; nkx < Nkint; nkx++)
	{
		double kx = kxmin + static_cast<double>(nkx)*dkx;
		double intky = 0;
		for (int nky = 0; nky < Nkint; nky++)
		{
			double ky = kymin + static_cast<double>(nky)*dky;
			double kk = sqrt(kx*kx + ky*ky);
			double kqx, kqy, kq;
			kqx = kx + qx;
			kqy = ky + qy;
			kq = sqrt(kqx*kqx + kqy*kqy);
			double Ekk, Ekq, fFDkk, fFDkq;
			Ekk = my_interp_eqdx(k_array, Ek_array, Nk, kk);
			Ekq = my_interp_eqdx(k_array, Ek_array, Nk, kq);
			fFDkk = 1.0/(exp((Ekk-muF)/(kB*T))+1.0);
			fFDkq = 1.0/(exp((Ekq-muF)/(kB*T))+1.0);
			double addintky;
			if (fabs(fFDkk-fFDkq) < 1e-6)
			{
				double f_1;
				f_1 = - 1.0/(exp(-(Ekk-muF)/(kB*T))+1.0);
				addintky = fFDkk*f_1/(kB*T);
				//addintky = fFDkk*(fFDkk-1.0)/(kB*T);
			}
			else
				addintky = (fFDkq-fFDkk)/(Ekq-Ekk);
			if (nky == 0 || nky == Nkint-1)
				addintky *= 0.5;
			intky += addintky*dky;
		}
		double addintkx;
		addintkx = intky;
		if (nkx == 0 || nkx == Nkint-1)
			addintkx *= 0.5;
		intkx += addintkx*dkx;
	}

	double intqkx = 0;
	for (int nqkx = 0; nqkx < Nkint; nqkx++)
	{
		double qkx = qkxmin + static_cast<double>(nqkx)*dqkx;
		double intqky = 0;
		for (int nqky = 0; nqky < Nkint; nqky++)
		{
			double qky = qkymin + static_cast<double>(nqky)*dqky;
			//if (qkx >= kxmin && qkx <= kxmax && qky >= kymin && qky <= kymax)
				//continue;
			double kk = sqrt(qkx*qkx + qky*qky);
			double kqx, kqy, kq;
			kqx = qkx + qx;
			kqy = qky + qy;
			kq = sqrt(kqx*kqx + kqy*kqy);
			double Ekk, Ekq, fFDkk, fFDkq;
			Ekk = my_interp_eqdx(k_array, Ek_array, Nk, kk);
			Ekq = my_interp_eqdx(k_array, Ek_array, Nk, kq);
			fFDkk = 1.0/(exp((Ekk-muF)/(kB*T))+1.0);
			fFDkq = 1.0/(exp((Ekq-muF)/(kB*T))+1.0);
			double addintqky;
			if (fabs(fFDkk-fFDkq) < 1e-6)
			{
				double f_1;
				f_1 = - 1.0/(exp(-(Ekk-muF)/(kB*T))+1.0);
				addintqky = fFDkk*f_1/(kB*T);
				//addintqky = fFDkk*(fFDkk-1.0)/(kB*T);
			}
			else
				addintqky = (fFDkq-fFDkk)/(Ekq-Ekk);
			if (nqky == 0 || nqky == Nkint-1)
				addintqky *= 0.5;
			intqky += addintqky*dqky;
		}
		double addintqkx;
		addintqkx = intqky;
		if (nqkx == 0 || nqkx == Nkint-1)
			addintqkx *= 0.5;
		intqkx += addintqkx*dqkx;
	}

	double Pol = 2.0/pow(2.0*PI,2.0)*(intkx + intqkx);
	return Pol;
}



double PolTeq0(double effm, double mu, double q)
{
	using namespace std;
	if ( q < 0)
	{
		cout << "q cannot be negative in calculating PolTeq0()!\n";
		exit(1);
	}
	if (mu < 0)
	{
		cout << "Error! Fermi level cannot be negative when T = 0!\n";
		exit(1);
	}
	if (mu==0)
		return 0.0;
	else
	{
		double xsq, pol;
		xsq = hbar*hbar*q*q/(8.0*effm*mu);
		if (xsq <= 1.0)
			pol = effm/(PI*hbar*hbar);
		else
			pol = effm/(PI*hbar*hbar)*(1.0 - sqrt(1.0-1.0/xsq));
		return -pol;
	}
}


double PolT(double effm, double mu, double T, double q)
{
	using namespace std;
	if ( q < 0)
	{
		cout << "q cannot be negative in calculating PolT()!\n";
		exit(1);
	}
	if (T==0)
	{
		cout << "Error!, function PolT() can only be used to calculate the finite temperature case!\n";
		exit(1);
	}
	double Dmupmin, Dmupmax, mupmin, mupmax, dmup;
	int Nmup = 5000;
	Dmupmin = - 10.0*(2.0*kB*T);
	Dmupmax = 10.0*(2.0*kB*T);
	mupmin = mu + Dmupmin;
	if (mupmin < 0)
		mupmin = 0;
	mupmax = mu + Dmupmax;
	dmup = (mupmax - mupmin)/static_cast<double>(Nmup-1);
	double intmup = 0;
	for (int nmup = 0; nmup < Nmup; nmup++)
	{
		double Dmup, mup, pol0mup, coshfactor;
		mup = mupmin + static_cast<double>(nmup)*dmup;
		Dmup = mup - mu;
		pol0mup = PolTeq0(effm, mup, q);
		coshfactor = cosh(Dmup/(2.0*kB*T));
		double addintmup;
		addintmup = pol0mup/pow(coshfactor,2.0);
		if (nmup == 0 || nmup == Nmup-1)
			addintmup *= 0.5;
		intmup += addintmup*dmup;
	}
	intmup /= 4.0*kB*T;
	return intmup;
}


void fsaveVector(char fVname[], double * V, int N)
{
	using namespace std;
	ofstream fout(fVname);
	if (!fout.is_open())
	{
		cerr << "Cannot open fVname file for output!" << endl;
		exit(1);
	}
	for (int n = 0; n < N; n++)
		fout << V[n] << "\t";
	fout << "\n\n";
	fout.close();
}