#include <inmost.h>
#include "grid.h"
#include "gmg_poisson.h"

using namespace INMOST;

double mu_s = 1.0;
int glev = 6;
int test = 2;
double We = 1.0;
double dt = 1.0, dt0 = 1.0;
double t = 0.0;

int main2(int argc, char ** argv)
{
	if(argc > 1)
		glev = atoi(argv[1]);
	pint::glev = glev;
	int nx[3], ny[3];
	for(int k = 0; k < 3; ++k)
	{
		nx[k] = pint::nx(1,k);
		ny[k] = pint::ny(2,k);
		std::cout << "level " << k << " nx " << nx[k] << " ny " << ny[k]  << std::endl;
		int level = k;
		int n1 = nx[level], n2 = ny[level];
		double hx = pint::hx(level), hy = pint::hy(level);
		double hx0 = pint::hx(0), hy0 = pint::hy(0);
		for(int i = 0; i < n1; ++i)
			for(int j = 1; j < n2; ++j)
		{
			pint pdof(1,i,j,level);
			assert(pdof.inside());
			if(i == 0) assert(fabs(pdof.x()-0.5*hx0) < 1.0e-6*hx);
			if(i == n1-1) assert(fabs(pdof.x()-(1-0.5*hx0)) < 1.0e-6*hx);
		}
		for(int j = 0; j < n2; ++j)
			for(int i = 1; i < n1; ++i)
		{
			pint pdof(2,i,j,level);
			assert(pdof.inside());
			if(j == 0) assert(fabs(pdof.y()-0.5*hy0) < 1.0e-6*hy);
			if(j == n2-1) assert(fabs(pdof.y()-(1-0.5*hy0)) < 1.0e-6*hy);
		}
		for(int i = 1; i < n1; ++i)
			for(int j = 1; j < n2; ++j)
		{
			pint pdof(3,i,j,level);
			assert(pdof.inside());
			if(i == 1) assert(fabs(pdof.x()-0.5*(hx+hx0)) < 1.0e-6*hx);
			if(i == n1-1) assert(fabs(pdof.x()-(1-0.5*(hx+hx0))) < 1.0e-6*hx);
			if(j == 1) assert(fabs(pdof.y()-0.5*(hy+hy0)) < 1.0e-6*hy);
			if(j == n2-1) assert(fabs(pdof.y()-(1-0.5*(hy+hy0))) < 1.0e-6*hy);
		}
	}

	return 0;
}

int main(int argc, char ** argv)
{
	Solver::Initialize(&argc, &argv,"database.xml");

	if(argc > 1)
		glev = atoi(argv[1]);
	pint::glev = glev;
	int nlevels = level0+4;
	int level = level0;

	test = 2;
	
	const int grids[1] = {3};
	const int ncomps[1] = {1};
	gmg_layout lay1(1,grids,ncomps, nlevels);
	
	const int grids2[1] = {3};	// Psi(3), Psi=(Psi_xx, Psi_yy, Psi_xy)
	const int ncomps2[1] = {3};
	gmg_layout lay2(1,grids2,ncomps2,1);

	multigrid_params common;
	common.nlevels = nlevels; common.maxiters = 30; common.atol = 1.0e-8;
	common.smooth_iters = 8; common.schedule = 2;
	common.integral_constraint = false;
	gmg_poisson gmg(common, lay1);
	dof_vector<double> &UVP = gmg.UVP, &RHS = gmg.RHS;
	dof_vector<double> UVP0(&lay1);

	bool success = gmg.solve();
	
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		UVP0(0,pdof) = fabs(UVP(0,pdof)-exact(test,3,pdof.x(),pdof.y()));
	}

	//save_vtk_poisson("tmp1.vtk", UVP);
	save_vtk_poisson("pois1.vtk", UVP0, level);
	
	Solver::Finalize();
	return 0;
}
