#ifndef _GMG_POISSON_H
#define _GMG_POISSON_H

#include "gmg.h"

extern double mu_s;
extern int test;

void fill_poisson(dof_vector<double>& RHS, multigrid_level* mgl, bool resid)
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
	// div grad u = 0, equations for Ph
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		int k = lay1.dof(pdof, 0);
		laplace(pdof,dofs1,dofs2,-mu_s);
		if(level == level0) RHS(0,pdof) = rhs_exact(test,3,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() )
			{
				if(mgcl && !filled)
					(*A)[k][lay1.dof(dofs1[q],0)] += dofs1[q].coef * scale;
				else if(resid)
					RHS(0,pdof) -= dofs1[q].coef*scale * UVP(0,dofs1[q]);
			}
			else if(level == level0) //if(test == 0)
				RHS(0,pdof) -= dofs1[q].coef * scale * exact(test,3,dofs1[q].x(),dofs1[q].y());
	}
}

void gs_smoother(int i, int j, dof_vector<double>& UVP, dof_vector<double>& RHS, int level)
{
	pint dofs1[5], dofs2[5];
	pint pdof(3,i,j,level);
	laplace(pdof, dofs1, dofs2, -mu_s);
	double val = RHS(0,pdof);
	if(level == level0) val += rhs_exact(test,pdof.grid,pdof.x(),pdof.y());
	for(int q = 0; q < 4; ++q)
		if(dofs1[q].inside())
			val -= dofs1[q].coef*UVP(0,dofs1[q]);
		else if(level == level0)
			val -= dofs1[q].coef*exact(test,dofs1[q].grid,dofs1[q].x(),dofs1[q].y());
	assert(pdof == dofs1[4]);
	UVP(0,pdof) = val / dofs1[4].coef;
}
void gs_smoother(int iters, dof_vector<double>& UVP, dof_vector<double>& RHS, int level)
{
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static bool fwd = 0;
	for(int k = 0; k < iters; ++k)
	{
		if(k % 2 == fwd)
		{
		for(int i = 1; i < n1; ++i)
			for(int j = 1; j < n2; ++j)
				gs_smoother(i,j,UVP,RHS,level);
		}
		else
		for(int i = n1-1; i >= 1; --i)
			for(int j = n2-1; j >= 1; --j)
				gs_smoother(i,j,UVP,RHS,level);
	}
	fwd = !fwd;
}

struct gmg_poisson : multigrid
{
	gmg_poisson(const multigrid_params& common, const gmg_layout& lay)
		: multigrid(common, lay) {}
	void fill_residual(dof_vector<double>& R, multigrid_level* mgl)
	{
		multigrid_coarse_level* mgcl = dynamic_cast<multigrid_coarse_level*>(mgl);
		if(!mgcl)
		{
			fill_poisson(R, mgl, true);
			return;
		}
		if(!mgcl->filled)
		{
			mgcl->setup();
			fill_poisson(R, mgcl, true);
			mgcl->S->SetMatrix(*(mgcl->A));
			mgcl->A->Save("A.mtx");
			mgcl->filled = true;
		}
		else
			fill_poisson(R, mgcl, true);
	}
	void smoothing(multigrid_level* mgl)
	{
		gs_smoother(common.smooth_iters, UVP, RHS, mgl->level);
	}
	void restriction(dof_vector<double>& U, multigrid_level* mgl);
	void prolongation(multigrid_level* mgl);
};

void gmg_poisson::restriction(dof_vector<double>& U, multigrid_level* mgl)
{
	int level = mgl->level;
	int n1 = pint::nx(1,level+1), n2 = pint::ny(2,level+1);
	// cycle over CV of coarse cells
	for(int iC = 1; iC < n1; ++iC)
		for(int jC = 1; jC < n2; ++jC)
	{
		int iF1 = 2*iC-1, iF2 = 2*iC, jF1 = 2*jC-1, jF2 = 2*jC;
		pint pC(3,iC,jC,level+1), pF11(3,iF1,jF1,level), pF12(3,iF1,jF2,level), pF21(3,iF2,jF1,level), pF22(3,iF2,jF2,level);
		U(0,pC) = 0.25*(U(0,pF11)+U(0,pF12)+U(0,pF21)+U(0,pF22));
	}
}

void gmg_poisson::prolongation(multigrid_level* mgl)
{
	dof_vector<double>& U = UVP;
	int level = mgl->level;
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static const double pcoef[4] = {9./16, 3./16, 3./16, 1./16};
	// cycle over CV of coarse cells
	for(int iC = 0; iC < n1; ++iC)
		for(int jC = 0; jC < n2; ++jC)
	{
		pint dofC00 = pint(3,iC,jC,level), dofC01 = pint(3,iC,jC+1,level), dofC10 = pint(3,iC+1,jC,level), dofC11 = pint(3,iC+1,jC+1,level);
		double pC00, pC01, pC10, pC11;
		if(dofC00.inside()) pC00 = U(0,dofC00);
		if(dofC01.inside()) pC01 = U(0,dofC01);
		if(dofC10.inside()) pC10 = U(0,dofC10);
		if(dofC11.inside()) pC11 = U(0,dofC11);
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
		if(dofF00.inside()) U(0,dofF00) += pF00;
		if(dofF01.inside()) U(0,dofF01) += pF01;
		if(dofF10.inside()) U(0,dofF10) += pF10;
		if(dofF11.inside()) U(0,dofF11) += pF11;
	}
}

#endif // _GMG_POISSON_H
