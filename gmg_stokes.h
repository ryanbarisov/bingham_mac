#ifndef _GMG_STOKES_H
#define _GMG_STOKES_H

#include "gmg.h"
#include "discr.h"
#include "vtk.h"

extern double mu_s;
extern bool divgrad;

// solve system with bordered matrix
// N=3:
//
// | D0       U0 | | SOL0 | = | RHS0 |
// |    D1    U1 | | SOL1 | = | RHS1 |
// |       D2 U2 | | SOL2 | = | RHS2 |
// | L0 L1 L2    | | SOL3 | = | RHS3 |
template<unsigned N>
bool solve_bordered(const double D[N-1], const double U[N-1], const double L[N-1], double RHS[N], double SOL[N])
{
	double numer = -RHS[N-1], denom = 0.0;
	bool fail = true;
	for(int i = 0; i < N-1; ++i)
	{
		assert(fabs(D[i]) > 1.0e-12);
		numer += L[i]/D[i] * RHS[i];
		denom += L[i]/D[i] * U[i];
		if(fabs(U[i]) > 1.0e-12) fail = false;
	}
	if(fail)
	{
		for(int i = 0; i < N-1; ++i)
		{
			assert(fabs(D[i]) > 1.0e-12);
			SOL[i] = RHS[i] / D[i];
		}
		SOL[N-1] = 0.0;
		return true;
	}
	assert(fabs(denom) > 1.0e-12);
	SOL[N-1] = numer/denom;
	for(int i = 0; i < N-1; ++i)
		SOL[i] = (RHS[i] - U[i] * SOL[N-1]) / D[i];
	return false;
}

void vanka_smoother(dof_vector<double>& UVP, const dof_vector<double>& RHS0, int i, int j, int level)
{
	double D[4] = {0,0,0,0}, U[4] = {0,0,0,0}, L[4], RHS[5] = {0,0,0,0,0}, SOL[5];
	pint dofs1[16], dofs2[16], dofs3[16], dof, puvp[5];
	puvp[0] = pint(1,i-1,j,level); puvp[1] = pint(1,i,j,level); puvp[2] = pint(2,i,j-1,level); puvp[3] = pint(2,i,j,level); puvp[4] = pint(3,i,j,level);

	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	double scale = hx*hy;
	int ndofs;

	for(int q = 0; q < 5; ++q)
	{
		assert(puvp[q].inside());
		RHS[q] = RHS0(puvp[q].grid-1, puvp[q])*scale;
		if(level == level0 && q != 4)
		       RHS[q] += rhs_exact(test,puvp[q].grid,puvp[q].x(),puvp[q].y()) * scale;
	}

	laplace(puvp[0],dofs1,dofs2, -mu_s);
	for(int q = 0; q < 5; ++q)
		if( dofs1[q].inside() && dofs1[q] == puvp[0])
			D[0] += dofs1[q].coef*scale;
		else if(dofs1[q].inside())
			RHS[0] -= dofs1[q].coef*scale*UVP(0,dofs1[q]);
		else if(level == level0)
			RHS[0] -= dofs1[q].coef*scale*exact(test,1,dofs1[q].x(),dofs1[q].y());
	ndofs = grad(puvp[0],dofs3);
	if(test != 0 || (puvp[0].i != 0 && puvp[0].i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
		for(int q = 0; q < ndofs; ++q)
			if( dofs3[q].inside() && dofs3[q] == puvp[4] )
				U[0] += dofs3[q].coef*scale;
			else if(dofs3[q].inside())
				RHS[0] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);


	laplace(puvp[1],dofs1,dofs2, -mu_s);
	for(int q = 0; q < 5; ++q)
		if( dofs1[q].inside() && dofs1[q] == puvp[1])
			D[1] += dofs1[q].coef*scale;
		else if(dofs1[q].inside())
			RHS[1] -= dofs1[q].coef*scale*UVP(0,dofs1[q]);
		else if(level == level0)
			RHS[1] -= dofs1[q].coef*scale*exact(test,1,dofs1[q].x(),dofs1[q].y());
	ndofs = grad(puvp[1],dofs3);
	if(test != 0 || (puvp[1].i != 0 && puvp[1].i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
		for(int q = 0; q < ndofs; ++q)
			if( dofs3[q].inside() && dofs3[q] == puvp[4] )
				U[1] += dofs3[q].coef*scale;
			else if(dofs3[q].inside())
				RHS[1] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);

	laplace(puvp[2],dofs1,dofs2, -mu_s);
	for(int q = 0; q < 5; ++q)
		if( dofs2[q].inside() && dofs2[q] == puvp[2])
			D[2] += dofs2[q].coef*scale;
		else if(dofs2[q].inside())
			RHS[2] -= dofs2[q].coef*scale*UVP(1,dofs2[q]);
		else if(level == level0)
			RHS[2] -= dofs2[q].coef*scale*exact(test,2,dofs2[q].x(),dofs2[q].y());
	ndofs = grad(puvp[2],dofs3);
	if(test != 0 || (puvp[2].j != 0 && puvp[2].j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
		for(int q = 0; q < ndofs; ++q)
			if( dofs3[q].inside() && dofs3[q] == puvp[4] )
				U[2] += dofs3[q].coef*scale;
			else if(dofs3[q].inside())
				RHS[2] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);
	
	laplace(puvp[3],dofs1,dofs2, -mu_s);
	for(int q = 0; q < 5; ++q)
		if( dofs2[q].inside() && dofs2[q] == puvp[3])
			D[3] += dofs2[q].coef*scale;
		else if(dofs2[q].inside())
			RHS[3] -= dofs2[q].coef*scale*UVP(1,dofs2[q]);
		else if(level == level0)
			RHS[3] -= dofs2[q].coef*scale*exact(test,2,dofs2[q].x(),dofs2[q].y());
	ndofs = grad(puvp[3],dofs3);
	if(test != 0 || (puvp[3].j != 0 && puvp[3].j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
		for(int q = 0; q < ndofs; ++q)
			if( dofs3[q].inside() && dofs3[q] == puvp[4] )
				U[3] += dofs3[q].coef*scale;
			else if(dofs3[q].inside())
				RHS[3] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);

	L[0] = -scale/hx;
	L[1] =  scale/hx;
	L[2] = -scale/hy;
	L[3] =  scale/hy;

	double alpha_uv = 1.0;
	for(int q = 0; q < 4; ++q) D[q] /= alpha_uv;

	bool fail = solve_bordered<5>(D,U,L,RHS,SOL);

	for(int q = 0; q < 5; ++q)
		UVP(puvp[q].grid-1,puvp[q]) = SOL[q];
}
void vanka_smoother(int iters, dof_vector<double>& UVP, const dof_vector<double>& RHS0, int level)
{
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static bool fwd = 0;
	for(int iter = 0; iter < iters; ++iter)
	{
		if(iter % 2 == fwd)
		//if(fwd)
		for(int i = 1; i < n1; ++i)
			for(int j = 1; j < n2; ++j)
				vanka_smoother(UVP, RHS0, i, j, level);
		else
		for(int i = n1-1; i >= 1; --i)
			for(int j = n2-1; j >= 1; --j)
				vanka_smoother(UVP, RHS0, i, j, level);
	}
	fwd = !fwd;
}

// distributive Gauss-Seidel https://www.math.uci.edu/~chenlong/226/MACcode.pdf
void gs_smoother(dof_vector<double>& UVP, const dof_vector<double>& RHS0, int level, dof_vector<double>& EP, dof_vector<double>& EP0, int iter = 0)
{
	pint dofs1[16], dofs2[16], dofs3[16];
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	int ndofs;
	double RHS, D;

	int iters = 1;
	bool swp = iter%2==0;

	// EP <- f - B^Tp^k
	EP.zero(level);
	int iU0 = 0, iU1 = n1, diU = 1; if(swp) {iU0 = n1-1; iU1 = -1; diU = -1;}
	int jU0 = 1, jU1 = n2, djU = 1; if(swp) {jU0 = n2-1; jU1 =  0; djU = -1;}
	int jV0 = 0, jV1 = n2, djV = 1; if(swp) {jV0 = n2-1; jV1 = -1; djV = -1;}
	int iV0 = 1, iV1 = n1, diV = 1; if(swp) {iV0 = n1-1; iV1 =  0; diV = -1;}
	int iP0 = 1, iP1 = n1, diP = 1; if(swp) {iP0 = n1-1; iP1 =  0; diP = -1;}
	int jP0 = 1, jP1 = n2, djP = 1; if(swp) {jP0 = n2-1; jP1 =  0; djP = -1;}

	//for(int i = 0; i < n1; ++i)
	//	for(int j = 1; j < n2; ++j)
	for(int i = iU0; i != iU1; i+=diU)
		for(int j = jU0; j != jU1; j+=djU)
	{
		pint pdof(1,i,j,level);
		RHS = RHS0(0,pdof);
		RHS += level == level0 ? rhs_exact(test,pdof.grid,pdof.x(),pdof.y()) : 0.0;
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (pdof.i != 0 && pdof.i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q)
				if(dofs3[q].inside())
					RHS -= dofs3[q].coef*UVP(2,dofs3[q]);
		EP(0,pdof) = RHS;
	}
	//for(int j = 0; j < n2; ++j)
	//	for(int i = 1; i < n1; ++i)
	for(int j = jV0; j != jV1; j+=djV)
		for(int i = iV0; i != iV1; i+=diV)
	{
		pint pdof(2,i,j,level);
		RHS = RHS0(1,pdof);
		RHS += level == level0 ? rhs_exact(test,pdof.grid,pdof.x(),pdof.y()) : 0.0;
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (pdof.j != 0 && pdof.j != n2-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q)
				if(dofs3[q].inside())
					RHS -= dofs3[q].coef*UVP(2,dofs3[q]);
		EP(1,pdof) = RHS;
	}
	// solve A (UV) = EP = r_u using Gauss-Seidel inner iterations
	// here A (UV) = -mu_s \Delta UV
	for(int it = 0; it < iters; ++it)
	{
		//for(int i = 0; i < n1; ++i)
		//	for(int j = 1; j < n2; ++j)
		for(int i = iU0; i != iU1; i+=diU)
			for(int j = jU0; j != jU1; j+=djU)
		{
			pint pdof(1,i,j,level);
			laplace(pdof,dofs1,dofs2,-mu_s);
			RHS = EP(0,pdof);
			for(int q = 0; q < 5; ++q)
			{
				if(dofs1[q] == pdof)
					D = dofs1[q].coef;
				else if(dofs1[q].inside())
					RHS -= dofs1[q].coef*UVP(0,dofs1[q]);
				else if(level == level0)
					RHS -= dofs1[q].coef*exact(test,dofs1[q].grid,dofs1[q].x(),dofs1[q].y());
			}
			UVP(0,pdof) = RHS/D;
		}
		//for(int j = 0; j < n2; ++j)
		//	for(int i = 1; i < n1; ++i)
		for(int j = jV0; j != jV1; j+=djV)
			for(int i = iV0; i != iV1; i+=diV)
		{
			pint pdof(2,i,j,level);
			laplace(pdof,dofs1,dofs2,-mu_s);
			RHS = EP(1,pdof);
			for(int q = 0; q < 5; ++q)
			{
				if(dofs2[q] == pdof)
					D = dofs2[q].coef;
				else if(dofs2[q].inside())
					RHS -= dofs2[q].coef*UVP(1,dofs2[q]);
				else if(level == level0)
					RHS -= dofs2[q].coef*exact(test,dofs2[q].grid,dofs2[q].x(),dofs2[q].y());
			}
			UVP(1,pdof) = RHS/D;
		}
	}
	// Delta e_p = g - div(UV)
	for(int it = 0; it < iters; ++it)
	{
		//for(int i = 1; i < n1; ++i)
		//	for(int j = 1; j < n2; ++j)
		for(int i = iP0; i != iP1; i+=diP)
			for(int j = jP0; j != jP1; j+=djP)
		{
			pint pdof(3,i,j,level);
			RHS = RHS0(2,pdof);
			ndofs = div(pdof, dofs2);
			for(int q = 0; q < ndofs; ++q)
				if(dofs2[q].inside())
					RHS -= dofs2[q].coef*UVP(dofs2[q].grid-1,dofs2[q]);
				else if(level == level0)
					RHS -= dofs2[q].coef*exact(test,dofs2[q].grid,dofs2[q].x(),dofs2[q].y());
			laplace(pdof,dofs1,dofs2, 1.0);// div(grad EP) = 1.0 * Delta(EP)
			for(int q = 0; q < 5; ++q)
			{
				if(dofs1[q] == pdof)
					D = dofs1[q].coef;
				else if(dofs1[q].inside())
					RHS -= dofs1[q].coef*EP(2,dofs1[q]);
			}
			EP(2,pdof) = RHS/D;
		}
	}
	// update U,V, and P
	// U += grad_x(EP)
	//for(int i = 0; i < n1; ++i)
	//	for(int j = 1; j < n2; ++j)
	for(int i = iU0; i != iU1; i+=diU)
		for(int j = jU0; j != jU1; j+=djU)
	{
		pint pdof(1,i,j,level);
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (i != 0 && i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q)
				if(dofs3[q].inside() )
					UVP(0,pdof) += dofs3[q].coef*EP(2,dofs3[q]);
	}
	// V += grad_y(EP)
	//for(int j = 0; j < n2; ++j)
	//	for(int i = 1; i < n1; ++i)
	for(int j = jV0; j != jV1; j+=djV)
		for(int i = iV0; i != iV1; i+=diV)
	{
		pint pdof(2,i,j,level);
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (j != 0 && j != n2-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dy != 0
			for(int q = 0; q < ndofs; ++q)
				if( dofs3[q].inside() )
					UVP(1,pdof) += dofs3[q].coef*EP(2,dofs3[q]);
	}
	// P -= Delta(EP) or P += div(UV) - g
	double Pmean = 0.0;
	//for(int i = 1; i < n1; ++i)
	//	for(int j = 1; j < n2; ++j)
	for(int i = iP0; i != iP1; i+=diP)
		for(int j = jP0; j != jP1; j+=djP)
	{
		pint pdof(3,i,j,level);
		RHS = -RHS0(2,pdof);
		ndofs = div(pdof,dofs1);
		for(int q = 0; q < ndofs; ++q)
			if(dofs1[q].inside())
				RHS += dofs1[q].coef*UVP(dofs1[q].grid-1, dofs1[q]);
			else if(level == level0)
				RHS += dofs1[q].coef*exact(test,dofs1[q].grid,dofs1[q].x(),dofs1[q].y());
		UVP(2,pdof) += RHS;
		/*RHS = 0.0;
		laplace(pdof,dofs1,dofs2,1.0);
		for(int q = 0; q < 5; ++q) if(dofs1[q].inside())
			RHS += dofs1[q].coef*EP(2,dofs1[q]);
		UVP(2,pdof) -= RHS;*/
		Pmean += UVP(2,pdof);
	}
	//for(int i = 1; i < n1; ++i)
	//	for(int j = 1; j < n2; ++j)
	for(int i = iP0; i != iP1; i+=diP)
		for(int j = jP0; j != jP1; j+=djP)
	{
		pint pdof(3,i,j,level);
		UVP(2,pdof) -= Pmean/((n1-1)*(n2-1));
	}
}
void gs_smoother(int iters, dof_vector<double>& UVP, const dof_vector<double>& RHS0, int level, dof_vector<double>& EP, dof_vector<double>& EP0)
{
	static bool fwd = true;
	for(int k = 0; k < iters; ++k)
		gs_smoother(UVP,RHS0,level,EP,EP0,fwd);
	fwd = !fwd;
}

// fill matrix corresponding to Stokes equations for U,V,P
// -mu \Delta u + \grad p = f
// div u = 0
// RHS <- exact() is needed only on the finest level
void fill_stokes(dof_vector<double>& RHS, multigrid_level* mgl, bool resid = false)
{
	int level = mgl->level;
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	double scale = pint::hx(level0)*pint::hy(level0);//hx*hy;
	const gmg_layout& lay1 = *(RHS.lay);
	dof_vector<double>& UVP = mgl->ptr->UVP;

	multigrid_coarse_level* mgcl = dynamic_cast<multigrid_coarse_level*>(mgl);
	bool filled = mgcl ? mgcl->filled : true;
	Sparse::Matrix* A = mgcl ? mgcl->A : NULL;

	pint dofs1[5], dofs2[5], dofs3[2];
	int ndofs;
	// equations for Uh
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(1,i,j,level);
		int k = lay1.dof(pdof, 0);
		laplace(pdof,dofs1,dofs2, -mu_s);
		if(level == level0) RHS(0,pdof) = rhs_exact(test, 1,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() )
			{
				if(mgcl && !filled)
					(*A)[k][lay1.dof(dofs1[q],0)] += dofs1[q].coef * scale;
				else if(resid)
					RHS(0,pdof) -= dofs1[q].coef*scale * UVP(0,dofs1[q]);
			}
			else if(level == level0) //if(test == 0)
				RHS(0,pdof) -= dofs1[q].coef * scale * exact(test, 1,dofs1[q].x(),dofs1[q].y());
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (i != 0 && i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
			{
				if(mgcl && !filled)
					(*A)[k][lay1.dof(dofs3[q],2)] += dofs3[q].coef * scale;
				else if(resid)
					RHS(0,pdof) -= dofs3[q].coef*scale * UVP(2,dofs3[q]);
			}
	}
	// equations for Vh
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		pint pdof(2,i,j,level);
		int k = lay1.dof(pdof, 1);
		laplace(pdof,dofs1,dofs2, -mu_s);
		if(level == level0) RHS(1,pdof) = rhs_exact(test, 2,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs2[q].inside() )
			{
				if(mgcl && !filled)
					(*A)[k][lay1.dof(dofs2[q],1)] += dofs2[q].coef * scale;
				else if(resid)
					RHS(1,pdof) -= dofs2[q].coef*scale * UVP(1,dofs2[q]);
			}
			else if(level == level0) //if(test == 0)
				RHS(1,pdof) -= dofs2[q].coef * scale * exact(test, 2,dofs2[q].x(),dofs2[q].y());
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (j != 0 && j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
			{
				if(mgcl && !filled)
					(*A)[k][lay1.dof(dofs3[q],2)] += dofs3[q].coef * scale;
				else if(resid)
					RHS(1,pdof) -= dofs3[q].coef*scale * UVP(2,dofs3[q]);
			}
	}
	// div u = 0, equations for Ph
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		int k = lay1.dof(pdof, 2);
		ndofs = div(pdof,dofs1);
		if(mgcl && !filled && mgcl->ptr->common.integral_constraint)
		{
			int nstokes = lay1.size(level);
			(*A)[k][nstokes] = scale;
			(*A)[nstokes][k] = scale;
		}
		for(int q = 0; q < ndofs; ++q)
			if( dofs1[q].inside() )
			{
				if(mgcl && !filled)
					(*A)[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
				else if(resid)
					RHS(2,pdof) -= dofs1[q].coef*scale * UVP(dofs1[q].grid-1,dofs1[q]);
			}
			else if(level == level0) //if(test == 0)
				RHS(2,pdof) -= dofs1[q].coef * scale * exact(test, dofs1[q].grid,dofs1[q].x(),dofs1[q].y());
	}
}


struct gmg_stokes : multigrid
{
	dof_vector<double> EP, EP0;
	gmg_stokes(const multigrid_params& common, const gmg_layout& lay)
		: multigrid(common, lay), EP(&lay), EP0(&lay) {}
	void fill_residual(dof_vector<double>& R, multigrid_level* mgl)
	{
		multigrid_coarse_level* mgcl = dynamic_cast<multigrid_coarse_level*>(mgl);
		if(!mgcl)
		{
			fill_stokes(R, mgl, true);
			return;
		}
		if(!mgcl->filled)
		{
			mgcl->setup();
			fill_stokes(R, mgcl, true);
			mgcl->S->SetMatrix(*(mgcl->A));
			mgcl->A->Save("A.mtx");
			mgcl->filled = true;
		}
		else
			fill_stokes(R, mgcl, true);
	}
	void smoothing(multigrid_level* mgl)
	{
		vanka_smoother(common.smooth_iters, UVP, RHS, mgl->level);
		// gauss-seidel smoother doesn't work for now:
		//gs_smoother(common.smooth_iters, UVP, RHS, mgl->level,EP,EP0);
	}
	void restriction(dof_vector<double>& U, multigrid_level* mgl);
	void prolongation(multigrid_level* mgl);
};


void gmg_stokes::restriction(dof_vector<double>& U, multigrid_level* mgl)
{
	int level = mgl->level;
	int n1 = pint::nx(1,level+1), n2 = pint::ny(2,level+1);
	static const double ucoef[2] = {1./8,1./4};
	// cycle over u-CV of coarse cells
	for(int iC = 0; iC < n1; ++iC)
		for(int jC = 1; jC < n2; ++jC)
	{
		int iF1 = 2*iC, iF0 = 2*iC-1, iF2 = 2*iC+1, jF0 = 2*jC-1, jF1 = 2*jC;
		pint pC(1,iC,jC,level+1), pF10(1,iF1,jF0,level), pF11(1,iF1,jF1,level);
		pint pF00(1,iF0,jF0,level), pF01(1,iF0,jF1,level);
		pint pF20(1,iF2,jF0,level), pF21(1,iF2,jF1,level);
		double uF10 = 0, uF11 = 0, uF00 = 0, uF01 = 0, uF20 = 0, uF21 = 0;
		uF10 = U(0,pF10); uF11 = U(0,pF11);
		if(pF00.inside()) uF00 = U(0,pF00);// else if(level==level0) uF00 = exact(test,1,pF00.x(),pF00.y());
		if(pF01.inside()) uF01 = U(0,pF01);// else if(level==level0) uF01 = exact(test,1,pF01.x(),pF01.y());
		if(pF20.inside()) uF20 = U(0,pF20);// else if(level==level0) uF20 = exact(test,1,pF20.x(),pF20.y());
		if(pF21.inside()) uF21 = U(0,pF21);// else if(level==level0) uF21 = exact(test,1,pF21.x(),pF21.y());
		// zero Neumann near boundary
		if(iC == 0)	{uF00 = uF10; uF01 = uF11;}
		if(jC == 0)	{uF00 = uF01; uF10 = uF11; uF20 = uF21;}
		if(iC == n1)	{uF20 = uF10; uF21 = uF11;}
		if(jC == n2-1)	{uF01 = uF00; uF11 = uF10; uF21 = uF20;}
		//U(0,pC) = 0.5*(uF10+uF11);
		U(0,pC) = ucoef[0]*(uF00+uF01+uF20+uF21) + ucoef[1]*(uF10+uF11);
	}
	// cycle over v-CV of coarse cells
	for(int jC = 0; jC < n2; ++jC)
		for(int iC = 1; iC < n1; ++iC)
	{
		int jF1 = 2*jC, jF0 = 2*jC-1, jF2 = 2*jC+1, iF0 = 2*iC-1, iF1 = 2*iC;
		pint pC(2,iC,jC,level+1), pF01(2,iF0,jF1,level), pF11(2,iF1,jF1,level);
		pint pF00(2,iF0,jF0,level), pF10(2,iF1,jF0,level);
		pint pF02(2,iF0,jF2,level), pF12(2,iF1,jF2,level);
		double vF00 = 0, vF10 = 0, vF01 = 0, vF11 = 0, vF02 = 0, vF12 = 0;
		vF01 = U(1,pF01); vF11 = U(1,pF11);
		if(pF00.inside()) vF00 = U(1,pF00);// else if(level==level0) vF00 = exact(test,2,pF00.x(),pF00.y());
		if(pF10.inside()) vF10 = U(1,pF10);// else if(level==level0) vF10 = exact(test,2,pF10.x(),pF10.y());
		if(pF02.inside()) vF02 = U(1,pF02);// else if(level==level0) vF02 = exact(test,2,pF02.x(),pF02.y());
		if(pF12.inside()) vF12 = U(1,pF12);// else if(level==level0) vF12 = exact(test,2,pF12.x(),pF12.y());
		// zero Neumann near boundary
		if(jC == 0)	{vF00 = vF01; vF10 = vF11;}
		if(iC == 0)	{vF00 = vF10; vF01 = vF11; vF02 = vF12;}
		if(jC == n2)	{vF02 = vF01; vF12 = vF11;}
		if(iC == n1-1)	{vF10 = vF00; vF11 = vF01; vF12 = vF02;}
		//U(1,pC) = 0.5*(U(1,pF01)+U(1,pF11));
		U(1,pC) = ucoef[0]*(vF00+vF10+vF02+vF12) + ucoef[1]*(vF01+vF11);
	}
	// cycle over p-CV of coarse cells
	for(int iC = 1; iC < n1; ++iC)
		for(int jC = 1; jC < n2; ++jC)
	{
		int iF1 = 2*iC-1, iF2 = 2*iC, jF1 = 2*jC-1, jF2 = 2*jC;
		pint pC(3,iC,jC,level+1), pF11(3,iF1,jF1,level), pF12(3,iF1,jF2,level), pF21(3,iF2,jF1,level), pF22(3,iF2,jF2,level);
		U(2,pC) = 0.25*(U(2,pF11)+U(2,pF12)+U(2,pF21)+U(2,pF22));
	}
}
void gmg_stokes::prolongation(multigrid_level* mgl)
{
	dof_vector<double>& U = UVP;
	int level = mgl->level;
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static const double ucoef[4] = {3./8, 3./8, 1./8, 1./8};
	static const double ucoef2[4] = {9./16, 3./32, 3./16, 1./32};//10.1016/j.jcp.2022.111500
	static const double pcoef[4] = {9./16, 3./16, 3./16, 1./16};
	// cycle over u-CV of coarse cells
	for(int iC = 0; iC < n1+1; ++iC)
		for(int jC = 0; jC < n2; ++jC)
	{
		pint dofC11 = pint(1,iC,jC,level), dofC12 = pint(1,iC,jC+1,level), dofC21 = pint(1,iC+1,jC,level), dofC22 = pint(1,iC+1,jC+1,level), dofC01 = pint(1,iC-1,jC,level), dofC02 = pint(1,iC-1,jC+1,level);
		double uC01 = 0, uC02 = 0, uC11 = 0, uC12 = 0, uC21 = 0, uC22 = 0;
		if(dofC01.inside()) uC01 = U(0,dofC01);
		if(dofC02.inside()) uC02 = U(0,dofC02);
		if(dofC11.inside()) uC11 = U(0,dofC11);
		if(dofC12.inside()) uC12 = U(0,dofC12);
		if(dofC21.inside()) uC21 = U(0,dofC21);
		if(dofC22.inside()) uC22 = U(0,dofC22);
		// use du/dn = 0 near boundary
		if(iC == 0)	{uC01 = uC11; uC02 = uC12;}
		if(jC == 0)	{uC01 = uC02; uC11 = uC12; uC21 = uC22;}
		if(iC == n1)	{uC21 = uC11; uC22 = uC12;}
		if(jC == n2-1)	{uC02 = uC01; uC12 = uC11; uC22 = uC21;}
		int iF = 2*iC, jF = 2*jC;
		pint dofF11 = pint(1,iF,jF,level-1), dofF12 = pint(1,iF,jF+1,level-1), dofF21 = pint(1,iF+1,jF,level-1), dofF22 = pint(1,iF+1,jF+1,level-1);
		double uF11, uF12, uF21, uF22;
		//uF11 = 2*ucoef[0]*uC11 + 2*ucoef[2]*uC12;
		//uF12 = 2*ucoef[0]*uC12 + 2*ucoef[2]*uC11;
		uF11 = ucoef2[0]*uC11 + ucoef2[1]*(uC21+uC01) + ucoef2[2]*uC12 + ucoef2[3]*(uC22+uC02);
		uF12 = ucoef2[0]*uC12 + ucoef2[1]*(uC22+uC02) + ucoef2[2]*uC11 + ucoef2[3]*(uC21+uC01);
		uF21 = ucoef[0]*uC11 + ucoef[1]*uC21 + ucoef[2]*uC12 + ucoef[3]*uC22;
		uF22 = ucoef[0]*uC12 + ucoef[1]*uC22 + ucoef[2]*uC11 + ucoef[3]*uC21;
		if(dofF11.inside()) U(0,dofF11) += uF11;
		if(dofF12.inside()) U(0,dofF12) += uF12;
		if(dofF21.inside()) U(0,dofF21) += uF21;
		if(dofF22.inside()) U(0,dofF22) += uF22;
	}
	// cycle over v-CV of coarse cells
	for(int iC = 0; iC < n1; ++iC)
		for(int jC = 0; jC < n2+1; ++jC)
	{
		pint dofC11 = pint(2,iC,jC,level), dofC12 = pint(2,iC,jC+1,level), dofC21 = pint(2,iC+1,jC,level), dofC22 = pint(2,iC+1,jC+1,level), dofC10 = pint(2,iC,jC-1,level), dofC20 = pint(2,iC+1,jC-1,level);
		double vC10 = 0, vC20 = 0, vC11 = 0, vC21 = 0, vC12 = 0, vC22 = 0;
		if(dofC10.inside()) vC10 = U(1,dofC10);
		if(dofC20.inside()) vC20 = U(1,dofC20);
		if(dofC11.inside()) vC11 = U(1,dofC11);
		if(dofC21.inside()) vC21 = U(1,dofC21);
		if(dofC12.inside()) vC12 = U(1,dofC12);
		if(dofC22.inside()) vC22 = U(1,dofC22);
		// use du/dn = 0 near boundary
		if(iC == 0)	{vC10 = vC20; vC11 = vC21; vC12 = vC22;}
		if(jC == 0)	{vC10 = vC11; vC20 = vC21;}
		if(iC == n1-1)	{vC20 = vC10; vC21 = vC11; vC22 = vC12;}
		if(jC == n2)	{vC12 = vC11; vC22 = vC21;}
		int iF = 2*iC, jF = 2*jC;
		pint dofF11 = pint(2,iF,jF,level-1), dofF12 = pint(2,iF,jF+1,level-1), dofF21 = pint(2,iF+1,jF,level-1), dofF22 = pint(2,iF+1,jF+1,level-1);
		double vF11, vF12, vF21, vF22;
		//vF11 = 2*ucoef[0]*vC11 + 2*ucoef[2]*vC21;
		//vF21 = 2*ucoef[0]*vC21 + 2*ucoef[2]*vC11;
		vF11 = ucoef2[0]*vC11 + ucoef2[1]*(vC12+vC10) + ucoef2[2]*vC21 + ucoef2[3]*(vC22+vC20);
		vF21 = ucoef2[0]*vC21 + ucoef2[1]*(vC22+vC20) + ucoef2[2]*vC11 + ucoef2[3]*(vC12+vC10);
		vF12 = ucoef[0]*vC11 + ucoef[1]*vC12 + ucoef[2]*vC21 + ucoef[3]*vC22;
		vF22 = ucoef[0]*vC21 + ucoef[1]*vC22 + ucoef[2]*vC11 + ucoef[3]*vC12;
		if(dofF11.inside()) U(1,dofF11) += vF11;
		if(dofF12.inside()) U(1,dofF12) += vF12;
		if(dofF21.inside()) U(1,dofF21) += vF21;
		if(dofF22.inside()) U(1,dofF22) += vF22;
	}
	// cycle over p-CV of coarse cells
	for(int iC = 0; iC < n1; ++iC)
		for(int jC = 0; jC < n2; ++jC)
	{
		pint dofC00 = pint(3,iC,jC,level), dofC01 = pint(3,iC,jC+1,level), dofC10 = pint(3,iC+1,jC,level), dofC11 = pint(3,iC+1,jC+1,level);
		double pC00 = 0, pC01 = 0, pC10 = 0, pC11 = 0;
		if(dofC00.inside()) pC00 = U(2,dofC00);
		if(dofC01.inside()) pC01 = U(2,dofC01);
		if(dofC10.inside()) pC10 = U(2,dofC10);
		if(dofC11.inside()) pC11 = U(2,dofC11);
		// use dp/dn = 0 near boundary
		if(iC == 0)	{pC00 = pC10; pC01 = pC11;}
		if(jC == 0)	{pC00 = pC01; pC10 = pC11;}
		if(iC == n1-1)	{pC10 = pC00; pC11 = pC01;}
		if(jC == n2-1)	{pC01 = pC00; pC11 = pC10;}
		int iF = 2*iC, jF = 2*jC;
		pint dofF00 = pint(3,iF,jF,level-1), dofF01 = pint(3,iF,jF+1,level-1), dofF10 = pint(3,iF+1,jF,level-1), dofF11 = pint(3,iF+1,jF+1,level-1);
		double pF00, pF01, pF10, pF11;
		pF00 = pcoef[0]*pC00 + pcoef[1]*pC01 + pcoef[2]*pC10 + pcoef[3]*pC11;
		pF01 = pcoef[0]*pC01 + pcoef[1]*pC00 + pcoef[2]*pC11 + pcoef[3]*pC10;
		pF10 = pcoef[0]*pC10 + pcoef[1]*pC00 + pcoef[2]*pC11 + pcoef[3]*pC01;
		pF11 = pcoef[0]*pC11 + pcoef[1]*pC01 + pcoef[2]*pC10 + pcoef[3]*pC00;
		if(dofF00.inside()) U(2,dofF00) += pF00;
		if(dofF01.inside()) U(2,dofF01) += pF01;
		if(dofF10.inside()) U(2,dofF10) += pF10;
		if(dofF11.inside()) U(2,dofF11) += pF11;
	}
}


static double func(int grid, double x, double y)
{
	if(grid == 1)
		return x+y+2;
	else if(grid == 2)
		return x-y-1;
	else if(grid == 3)
		return x;
	else return 0.0;
}
int test_gmg1()
{
	//n1 = pow(2, glev)+1;//atoi(argv[1]);
	//n2 = n1;//argv[2];
	int nlevels = 3;

	const int grids1[3] = {1,2,3};	// Uh(1), Vh(2), Ph(3)
	const int ncomps1[3] = {1,1,1};
	const int grids2[1] = {3};	// Psi(3), Psi=(Psi_xx, Psi_yy, Psi_xy)
	const int ncomps2[1] = {3};

	gmg_layout lay1(3,grids1,ncomps1,nlevels), lay2(1,grids2,ncomps2,1), lay3(3,grids1,ncomps1,1);
	dof_vector<double> UVP(&lay1), UVP0(&lay1), RHS(&lay1);
	dof_vector<symmat2> PSI(&lay2), TAU(&lay2);

	int level = 0;
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(1,i,j,level);
		UVP(0,pdof) = func(pdof.grid,pdof.x(),pdof.y());
	}
	for(int i = 1; i < n1; ++i)
		for(int j = 0; j < n2; ++j)
	{
		pint pdof(2,i,j,level);
		UVP(1,pdof) = func(pdof.grid,pdof.x(),pdof.y());
	}
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		UVP(2,pdof) = func(pdof.grid,pdof.x(),pdof.y());
	}
	save_vtk("tmp2.vtk", PSI, TAU, UVP);
	UVP0.assign(UVP, 0);

	multigrid_params common; common.nlevels = nlevels;
	gmg_stokes gmg(common, lay1);
	gmg.restriction(UVP, gmg.head);
	save_vtk("tmp3.vtk", PSI, TAU, UVP, 1);
	gmg.restriction(UVP, gmg.head->next);
	save_vtk("tmp4.vtk", PSI, TAU, UVP, 2);
	UVP.zero(2);
	gmg.prolongation(gmg.head->next->next);
	save_vtk("tmp5.vtk", PSI, TAU, UVP, 1);
	UVP.zero(1);
	gmg.prolongation(gmg.head->next);
	save_vtk("tmp6.vtk", PSI, TAU, UVP);

	for(int q = 0; q < UVP.size(); ++q)
		if(fabs(UVP.U[q] - UVP0.U[q]) > 1.0e-12)
			std::cerr << "component " << q
				<< " difference " << fabs(UVP.U[q] - UVP0.U[q])
				<< std::endl;
	return 0;
}

int test_gmg2()
{
	double D1[4] = {1,2,3,4};
	double U1[4] = {5,6,7,8};
	//double L1[4] = {9,10,11,12};
	//double RHS1[5] = {-26,-26,-44,-24,26};
	double L1[4] = {-9,9,-1,1};
	double RHS1[5] = {-26,-26,-44,-24,34};
	double SOL1[5] = {-1,2,-3,4,-5}, SOL2[5];
	solve_bordered<5>(D1,U1,L1,RHS1,SOL2);
	for(int k = 0; k < 5; ++k)
		if(fabs(SOL1[k]-SOL2[k]) > 1.0e-8)
			std::cout << "FAIL: solve_bordered k " << k
				<< " expected " << SOL1[k]
				<< " computed " << SOL2[k]
				<< std::endl;
	
	return 0;
}

#endif // _GMG_STOKES_H
