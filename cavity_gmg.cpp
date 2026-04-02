//#define USE_S3M

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "analytics.h"
#include "grid.h"
#include "mat2.h"
#include "discr.h"
#include "gmg_stokes.h"
#include "vtk.h"



//doi: 10.1016/j.jnnfm.2004.12.003

double mu_s = 1.0; // solvent viscosity
double mu_p = mu_s; // polymer viscosity
double We = 1.0;

int glev = 6; // n1 = 2^glev+1
double t = 0, T = 1.0;
double dt = 1.0, dt0 = 0.0;
double cfl = 0.5; // Kurganov-Tadmor scheme is stable for cfl <= 0.5
int nframe = 64;

int test = 0; // 0 -- regularized cavity, 1,2 -- analytical tests for Stokes problem alone
bool divgrad = false; // Delta_h = 2 div_h D_h (true), standard five-point Laplace (false)


double cfl_time_step(double t, double dt, int level)
{
	double maxu = max_exact(test, 1);
	double hx = pint::hx(level), hy = pint::hy(level);
	double h = std::min(hx,hy);
	double dt_cfl = cfl/(maxu+1.0e-12) * h;
	return dt_cfl;
}

void update_time()
{
	dt0 = dt;
	t += dt;
	int level = 0; // TODO
	dt = std::min(dt, cfl_time_step(t, dt, level));
}


void setup_exact(dof_vector<double>& UVP)
{
	int level = 0;
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(1,i,j,level);
		UVP(0,pdof) = exact(test, pdof.grid,pdof.x(),pdof.y());
	}
	for(int i = 1; i < n1; ++i)
		for(int j = 0; j < n2; ++j)
	{
		pint pdof(2,i,j,level);
		UVP(1,pdof) = exact(test, pdof.grid,pdof.x(),pdof.y());
	}
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		UVP(2,pdof) = exact(test, pdof.grid,pdof.x(),pdof.y());
	}
}

int main2()
{
	pint::glev = glev;
	test_gmg1();
	test_gmg2();
	return 0;
}

int main(int argc, char ** argv)
{
	Solver::Initialize(&argc, &argv, "database.xml");
	if(argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " glev [T=" << T
			<< "] [nframe=" << nframe
			<< "] [We=" << We
			<< "]\n";
		std::cout << "NOTE for parameter levels: n1 = 2^(glev)+1 is amount of nodes in X and Y direction!\n";
		return 0;
	}
	
	glev = atoi(argv[1]);
	if(glev < 2) glev = 2;
	pint::glev = glev;
	if(argc > 2)
		T = atof(argv[2]);
	if(argc > 3)
		nframe = atoi(argv[3]);
	if(argc > 4)
		We = atof(argv[4]);

	//test = 1;
	
	int nlevels = level0+3;
	int level = level0; // TODO
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	dt = cfl * std::min(hx,hy);
	std::cout << "n1 " << pint::nx(1,level) << " n2 " << pint::ny(2,level) << " test " << test << " dt " << dt << std::endl;
	std::cout << "nframe " << nframe << " T " << T << " We " << We << std::endl;
	
	const int grids1[3] = {1,2,3};	// Uh(1), Vh(2), Ph(3)
	const int ncomps1[3] = {1,1,1};
	const int grids2[1] = {3};	// Psi(3), Psi=(Psi_xx, Psi_yy, Psi_xy)
	const int ncomps2[1] = {3};

	gmg_layout lay1(3,grids1,ncomps1,nlevels), lay2(1,grids2,ncomps2,1), lay3(3,grids1,ncomps1,1);
	dof_vector<double> UVP0(&lay1), RHS0(&lay3);
	dof_vector<symmat2> PSI(&lay2), PSI0(&lay2), PSI1(&lay2), TAU(&lay2), TAU0(&lay2), TAU1(&lay2);

	int nstokes = lay1.size(level);//n1*(n2-1)+n2*(n1-1)+(n1-1)*(n2-1); // Uh,Vh,Ph
	int ntau = lay2.size(0);//3*(n1-1)*(n2-1);
	std::cout << "nstokes " << nstokes << " ntau " << ntau << std::endl;
	// psi = (psi_xx (Ph), psi_yy (Ph), psi_xy (Ph))
	// psi -- psi(n), psi0 -- psi(n-1), psi1 -- psi(n+1)
	// tau -- tau(n+1), tau0 -- tau(n-1), tau1 -- tau(n)

	bool success = false, finish = false;
	double res = 0.0, atol = 1.0e-4;
	int iter = 0, maxiters = 10000;

	multigrid_params common;
	common.nlevels = nlevels; common.maxiters = 30;
	common.atol = 1.0e-6; common.rtol = 1.0e-12;
	common.smooth_iters = 4; common.schedule = 2;
	//common.integral_constraint = false;
	gmg_stokes gmg(common, lay1);
	dof_vector<double> &UVP = gmg.UVP, &RHS = gmg.RHS;
	if(test) setup_exact(UVP);
	save_vtk("tmp0.vtk", PSI, TAU, UVP);
	UVP.zero(0);
	do
	{
		RHS.zero(level);
		if(test == 0) fill_rhs(TAU, RHS, mu_p); // RHS <- RHS + div(TAU)
		success = gmg.solve();
		if(test)
		{
			finish = true;
			for(int i = 0; i < n1; ++i)
				for(int j = 1; j < n2; ++j)
			{
				pint pdof(1,i,j,level);
				UVP0(0,pdof) = fabs(UVP(0,pdof)-exact(test,1,pdof.x(),pdof.y()));
			}
			for(int j = 0; j < n2; ++j)
				for(int i = 1; i < n1; ++i)
			{
				pint pdof(2,i,j,level);
				UVP0(1,pdof) = fabs(UVP(1,pdof)-exact(test,2,pdof.x(),pdof.y()));
			}
			for(int i = 1; i < n1; ++i)
				for(int j = 1; j < n2; ++j)
			{
				pint pdof(3,i,j,level);
				UVP0(2,pdof) = fabs(UVP(2,pdof)-exact(test,3,pdof.x(),pdof.y()));
			}
			save_vtk("tmp1.vtk", PSI, TAU, UVP0, level);
			save_vtk("tmp2.vtk", PSI, TAU, UVP);
			//save_vtk("sol"+std::to_string(iter)+".vtk", PSI, TAU, UVP);
			//++iter;
		}
		else if(success)
		{
			update_psi(UVP, UVP0, PSI, PSI0, PSI1);
			UVP0.assign(UVP, level);
			update_tau(PSI, TAU1, TAU);
			res = residual(TAU0, TAU);
			TAU0.assign(TAU1, level);
			std::cout << "t " << t << "/" << T << " iter " << iter << " " << " resid " << res << std::endl;
			//std::cin.ignore();
			if(iter % nframe == 0)
			{
				int frame = iter/nframe;
				save_vtk("out" + std::to_string(frame) + ".vtk", PSI, TAU, UVP);
			}
			++iter;
			finish = res < atol && t > T;
			if(iter >= maxiters || !success)
			{
				finish = false;
				break;
			}
			update_time();
		}
		else break;
	}
	while(!finish);

	Solver::Finalize();
	return 0;
}
