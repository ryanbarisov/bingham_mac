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

int main(int argc, char ** argv)
{
	Solver::Initialize(&argc, &argv,"database.xml");

	if(argc > 1)
		glev = atoi(argv[1]);
	pint::glev = glev;
	int nlevels = level0+1;
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
	common.smooth_iters = 4; common.schedule = 2;
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
