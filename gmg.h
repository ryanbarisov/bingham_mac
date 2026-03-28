#ifndef _GMG_H
#define _GMG_H

#include "vtk.h"

#if defined(USE_S3M)
#include "s3m/csrmatrix.h"
#include "s3m/bpcg.h"
#include "s3m/pcg.h"
#include "s3m/bicgstab.h"
#include "s3m/amg_ruge_stuben.h"
#include "s3m/gauss_seidel.h"
#include "s3m/ilduc.h"
#include "s3m/vanka.h"
#else
#include <inmost.h>
using namespace INMOST;
#endif


extern double mu_s;
extern bool divgrad;

void save_vtk(std::string filename, const dof_vector<double>& x, const gmg_layout& lay2, int level = 0, const dof_vector<symmat2>* psi = NULL, const dof_vector<symmat2>* tau = NULL)
{
	if(psi == NULL || tau == NULL)
	{
		dof_vector<symmat2> psi(&lay2);
		save_vtk(filename, psi, psi, x, level);
	}
	else save_vtk(filename, *psi, *tau, x, level);
}

void fill_vector(const dof_vector<double>& in, Sparse::Vector& out, int level = 0)
{
	//std::copy(in.U.begin(), in.U.end(), out.Begin());
	std::copy(in.U.begin()+in.offset(level), in.U.begin()+in.offset(level)+in.size(level), out.Begin());
}
void fill_vector(const Sparse::Vector& in, dof_vector<double>& out, int level = 0)
{
	//std::copy(in.Begin(), in.End(), out.U.begin());
	//std::copy(in.Begin(), in.End(), out.U.begin()+out.offset(level));
	for(int k = 0; k < out.size(level); ++k)
		out.U[out.offset(level)+k] = in[k];
	//std::copy(in.Begin(), in.Begin()+out.offset(level), out.U.begin()+out.offset(level));
}
void fill_vector(const dof_vector<double>& in, std::vector<double>& out, int level = 0)
{
	//std::copy(in.U.begin(), in.U.end(), out.begin());
	std::copy(in.U.begin()+in.offset(level), in.U.begin()+in.offset(level)+in.size(level), out.begin());
}
void fill_vector(const std::vector<double>& in, dof_vector<double>& out, int level = 0)
{
	//std::copy(in.begin(), in.end(), out.U.begin());
	std::copy(in.begin(), in.end(), out.U.begin()+out.offset(level));
}

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
	double scale = 1.0;//hx*hy;//pint::hx(0)*pint::hy(0);//1.0;//hx*hy;
	int ndofs;

	for(int q = 0; q < 5; ++q) if(puvp[q].inside())
		RHS[q] = RHS0(puvp[q].grid-1, puvp[q])*scale;

	bool divfree = false;
	if(divfree && i > 1)
	{
		D[0] = 1.0; U[0] = 0.0; RHS[0] = UVP(0,puvp[0]);
	}
	else
	{
		laplace(puvp[0],dofs1,dofs2, -mu_s);
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() && dofs1[q] == puvp[0])
				D[0] += dofs1[q].coef*scale;
			else if(dofs1[q].inside())
				RHS[0] -= dofs1[q].coef*scale*UVP(0,dofs1[q]);
			else if(level == 0)
				RHS[0] -= dofs1[q].coef*scale*exact(test,1,dofs1[q].x(),dofs1[q].y());
		if(divgrad) for(int q = 0; q < 4; ++q)
			if(dofs2[q].inside())
				RHS[0] -= dofs2[q].coef*scale*UVP(1,dofs2[q]);
			else if(level == 0)
				RHS[0] -= dofs2[q].coef*scale*exact(test,2,dofs2[q].x(),dofs2[q].y());
		ndofs = grad(puvp[0],dofs3);
		if(test != 0 || (puvp[0].i != 0 && puvp[0].i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q)
				if( dofs3[q].inside() && dofs3[q] == puvp[4] )
					U[0] += dofs3[q].coef*scale;
				else if(dofs3[q].inside())
					RHS[0] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);
	}

	laplace(puvp[1],dofs1,dofs2, -mu_s);
	for(int q = 0; q < 5; ++q)
		if( dofs1[q].inside() && dofs1[q] == puvp[1])
			D[1] += dofs1[q].coef*scale;
		else if(dofs1[q].inside())
			RHS[1] -= dofs1[q].coef*scale*UVP(0,dofs1[q]);
		else if(level == 0)
			RHS[1] -= dofs1[q].coef*scale*exact(test,1,dofs1[q].x(),dofs1[q].y());
	if(divgrad) for(int q = 0; q < 4; ++q)
		if(dofs2[q].inside())
			RHS[1] -= dofs2[q].coef*scale*UVP(1,dofs2[q]);
		else if(level == 0)
			RHS[1] -= dofs2[q].coef*scale*exact(test,2,dofs2[q].x(),dofs2[q].y());
	ndofs = grad(puvp[1],dofs3);
	if(test != 0 || (puvp[1].i != 0 && puvp[1].i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
		for(int q = 0; q < ndofs; ++q)
			if( dofs3[q].inside() && dofs3[q] == puvp[4] )
				U[1] += dofs3[q].coef*scale;
			else if(dofs3[q].inside())
				RHS[1] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);

	if(divfree && j > 1)
	{
		D[2] = 1.0; U[2] = 0.0; RHS[2] = UVP(1,puvp[2]);
	}
	else
	{
		laplace(puvp[2],dofs1,dofs2, -mu_s);
		for(int q = 0; q < 5; ++q)
			if( dofs2[q].inside() && dofs2[q] == puvp[2])
				D[2] += dofs2[q].coef*scale;
			else if(dofs2[q].inside())
				RHS[2] -= dofs2[q].coef*scale*UVP(1,dofs2[q]);
			else if(level == 0)
				RHS[2] -= dofs2[q].coef*scale*exact(test,2,dofs2[q].x(),dofs2[q].y());
		if(divgrad) for(int q = 0; q < 4; ++q)
			if(dofs1[q].inside())
				RHS[2] -= dofs1[q].coef*scale*UVP(0,dofs1[q]);
			else if(level == 0)
				RHS[2] -= dofs1[q].coef*scale*exact(test,1,dofs1[q].x(),dofs1[q].y());
		ndofs = grad(puvp[2],dofs3);
		if(test != 0 || (puvp[2].j != 0 && puvp[2].j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
			for(int q = 0; q < ndofs; ++q)
				if( dofs3[q].inside() && dofs3[q] == puvp[4] )
					U[2] += dofs3[q].coef*scale;
				else if(dofs3[q].inside())
					RHS[2] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);
	}
	laplace(puvp[3],dofs1,dofs2, -mu_s);
	for(int q = 0; q < 5; ++q)
		if( dofs2[q].inside() && dofs2[q] == puvp[3])
			D[3] += dofs2[q].coef*scale;
		else if(dofs2[q].inside())
			RHS[3] -= dofs2[q].coef*scale*UVP(1,dofs2[q]);
		else if(level == 0)
			RHS[3] -= dofs2[q].coef*scale*exact(test,2,dofs2[q].x(),dofs2[q].y());
	if(divgrad) for(int q = 0; q < 4; ++q)
		if(dofs1[q].inside())
			RHS[3] -= dofs1[q].coef*scale*UVP(0,dofs1[q]);
		else if(level == 0)
			RHS[3] -= dofs1[q].coef*scale*exact(test,1,dofs1[q].x(),dofs1[q].y());
	ndofs = grad(puvp[3],dofs3);
	if(test != 0 || (puvp[3].j != 0 && puvp[3].j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
		for(int q = 0; q < ndofs; ++q)
			if( dofs3[q].inside() && dofs3[q] == puvp[4] )
				U[3] += dofs3[q].coef*scale;
			else if(dofs3[q].inside())
				RHS[3] -= dofs3[q].coef*scale*UVP(2,dofs3[q]);

	L[0] = -1.0/hx;
	L[1] =  1.0/hx;
	L[2] = -1.0/hy;
	L[3] =  1.0/hy;

	double alpha_uv = 0.99;
	for(int q = 0; q < 4; ++q) D[q] /= alpha_uv;

	bool fail = solve_bordered<5>(D,U,L,RHS,SOL);
	UVP(0,puvp[0]) = SOL[0];
	UVP(0,puvp[1]) = SOL[1];
	UVP(1,puvp[2]) = SOL[2];
	UVP(1,puvp[3]) = SOL[3];
	/*if(!fail)*/ UVP(2,puvp[4]) = SOL[4];
}

void vanka_smoother(int iters, dof_vector<double>& UVP, const dof_vector<double>& RHS0, int level)
{
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static bool fwd = 0;
	for(int iter = 0; iter < iters; ++iter)
	{
		if(iter % 2 == !fwd)
		//if(fwd)
		for(int j = 1; j < n2; ++j)
			for(int i = 1; i < n1; ++i)
				vanka_smoother(UVP, RHS0, i, j, level);
		else
		for(int j = n2-1; j >= 1; --j)
			for(int i = n1-1; i >= 1; --i)
				vanka_smoother(UVP, RHS0, i, j, level);
	}
	fwd = !fwd;
}

typedef struct multigrid_params
{
	int nlevels = 2;
	int schedule = 1; // 1 -- V, 2 -- W
	int presmooth_iters = 2, postsmooth_iters = 2;
	int maxiters = 200;
	double rtol = 1.0e-6, atol = 1.0e-8;

	dof_vector<double> *UVP = NULL, *RHS0 = NULL;
	const gmg_layout* lay2;

	bool integral_constraint = true;

	void setup() {}
} multigrid_params;

typedef struct multigrid_level
{
	int level = 0;
	multigrid_params common;
	struct multigrid_level* next = NULL;
	
	multigrid_level(int level, const multigrid_params& common);

	multigrid_level(const multigrid_level& other)
		: level(other.level), common(other.common)
	{
		if(other.next != NULL) next = new multigrid_level(*(other.next));
	}
	multigrid_level& operator=(const multigrid_level& other)
	{
		level = other.level;
		common = other.common;
		if(other.next != NULL) next = new multigrid_level(*(other.next));
		return *this;
	}
	virtual ~multigrid_level() {if(next != NULL) delete next;}

	void smoothing(int iters)
	{
		vanka_smoother(iters, *(common.UVP), *(common.RHS0), level);
	}
	void restriction(dof_vector<double>& U);
	void prolongation(dof_vector<double>& U);
	void update();

	virtual void solve();
} multigrid_level;

struct multigrid_coarse_level : multigrid_level
{
	bool filled = false;
#if defined(USE_S3M)
	CSRMatrix* A;
	std::vector<double> x, b;
	//BramblePasciakCG< AMGRugeStuben< GaussSeidel, BICGSTAB<ILDUC> > > S;
	BICGSTAB< Vanka > S;
#else
	Sparse::Matrix *A;
	Sparse::Vector *b, *x;
	Solver S;
#endif
	
	multigrid_coarse_level(int level, multigrid_params common) : multigrid_level(level, common)
#if !defined(USE_S3M)
	, S(Solver::INNER_ILU2, "test")
#endif
	{
#if defined(USE_S3M)
		S.GetParameters().Load("params.txt");
		S.GetParameters().Set<int>("B_block_end", lay1.offset(2));//Uh+Vh
#endif
	}
	~multigrid_coarse_level()
	{
		if(filled)
		{
			delete A;
			delete b;
			delete x;
		}
	}
	void setup()
	{
		if(!filled)
		{
			int nstokes = common.UVP->size(level);
			if(common.integral_constraint) nstokes += 1; // integral constraint on pressure
#if defined(USE_S3M)
			A = new CSRMatrix();
			A->ReserveSize(nstokes);
			x = new std::vector<double>(nstokes);
			b = new std::vector<double>(nstokes);
#else
			A = new Sparse::Matrix("stokes", 0, nstokes);
			b = new Sparse::Vector("stokes", 0, nstokes);
			x = new Sparse::Vector("stokes", 0, nstokes);
#endif
		}
	}
	virtual void solve();
};

static multigrid_level* create_level(int level, const multigrid_params& common)
{
	if(level == common.nlevels - 1)
		return new multigrid_coarse_level(level, common);
	else if(level < common.nlevels - 1)
		return new multigrid_level(level, common);
	else return NULL;
}

typedef struct multigrid
{
	multigrid_params common;
	multigrid_level* head;
	dof_vector<double> res;

	multigrid(multigrid_params common) : common(common), res(common.UVP->lay)
	{
		head = create_level(0, common);
	}
	void step() { head->solve(); }
	double residual();
	bool solve()
	{
		bool success = false;
		double resid = 0.0, resid0 = 0.0;
		int iter = 0;
		// solve in a loop, using parameters from common
		do
		{
			step();
			resid = residual();
			if(iter == 0) resid0 = resid;
			if(resid/resid0 < common.rtol || resid < common.atol) success = true;
			else if(iter > common.maxiters) success = false;
			++iter;
			std::cout << "gmg iter " << iter << " resid " << resid << " resid/resid0 " << (resid/resid0) << std::endl;
			//success = true;
		} while(!success);
		std::cout << "gmg: solved? " << success << " iters " << iter << std::endl;
		return success;
	}
} multigrid_head;

multigrid_level::multigrid_level(int level, const multigrid_params& common) : level(level), common(common)
{
	next = create_level(level+1, common);
}
void multigrid_level::restriction(dof_vector<double>& U)
{
	int n1 = pint::nx(1,level+1), n2 = pint::ny(2,level+1);
	// cycle over u-CV of coarse cells
	for(int iC = 0; iC < n1; ++iC)
		for(int jC = 1; jC < n2; ++jC)
	{
		int iF = 2*iC, jF1 = 2*jC-1, jF2 = 2*jC;
		pint pC(1,iC,jC,level+1), pF1(1,iF,jF1,level), pF2(1,iF,jF2,level);
		U(0,pC) = 0.5*(U(0,pF1)+U(0,pF2));
	}
	// cycle over v-CV of coarse cells
	for(int iC = 1; iC < n1; ++iC)
		for(int jC = 0; jC < n2; ++jC)
	{
		int jF = 2*jC, iF1 = 2*iC-1, iF2 = 2*iC;
		pint pC(2,iC,jC,level+1), pF1(2,iF1,jF,level), pF2(2,iF2,jF,level);
		U(1,pC) = 0.5*(U(1,pF1)+U(1,pF2));
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
void multigrid_level::prolongation(dof_vector<double>& U)
{
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static const double ucoef[4] = {3./8, 3./8, 1./8, 1./8};
	static const double ucoef2[4] = {9./16, 3./32, 3./16, 1./32};//10.1016/j.jcp.2022.111500
	static const double pcoef[4] = {9./16, 3./16, 3./16, 1./16};
	// cycle over u-CV of coarse cells
	for(int iC = 0; iC < n1+1; ++iC)
		for(int jC = 0; jC < n2; ++jC)
	{
		pint dofC11 = pint(1,iC,jC,level), dofC12 = pint(1,iC,jC+1,level), dofC21 = pint(1,iC+1,jC,level), dofC22 = pint(1,iC+1,jC+1,level), dofC01 = pint(1,iC-1,jC,level), dofC02 = pint(1,iC-1,jC+1,level);
		double uC01, uC02, uC11, uC12, uC21, uC22;
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
		double vC10, vC20, vC11, vC21, vC12, vC22;
		if(dofC10.inside()) vC10 = U(1,dofC10);
		if(dofC20.inside()) vC20 = U(1,dofC20);
		if(dofC11.inside()) vC11 = U(1,dofC11);
		if(dofC21.inside()) vC21 = U(1,dofC21);
		if(dofC12.inside()) vC12 = U(1,dofC12);
		if(dofC22.inside()) vC22 = U(1,dofC22);	
		// use du/dn = 0 near boundary
		if(iC == 0)	{vC10 = vC20; vC11 = vC21; vC12 = vC22;}
		if(jC == 0)	{vC10 = vC11; vC20 = vC21;}
		if(iC == n1)	{vC20 = vC10; vC21 = vC11; vC22 = vC12;}
		if(jC == n2-1)	{vC12 = vC11; vC22 = vC21;}
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
		double pC00, pC01, pC10, pC11;
		if(dofC00.inside()) pC00 = U(2,dofC00);
		if(dofC01.inside()) pC01 = U(2,dofC01);
		if(dofC10.inside()) pC10 = U(2,dofC10);
		if(dofC11.inside()) pC11 = U(2,dofC11);
		// use dp/dn = 0 near boundary
		if(iC == 0)	{pC00 = pC10; pC01 = pC11;}
		if(jC == 0)	{pC00 = pC01; pC10 = pC11;}
		if(iC == n1)	{pC10 = pC00; pC11 = pC01;}
		if(jC == n2)	{pC01 = pC00; pC11 = pC10;}
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


#if defined(USE_S3M)
typedef CSRMatrix matrix_type;
#else
typedef Sparse::Matrix matrix_type;
#endif

// fill matrix corresponding to Stokes equations for U,V,P
// -mu \Delta u + \grad p = f
// div u = 0
// TODO RHS0 <- exact() is needed only on finest level?
void fill_stokes(dof_vector<double>& RHS0, multigrid_level* mgl, bool resid = false)
{
	int level = mgl->level;
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	double scale = 1.0;//hx*hy;//pint::hx(0)*pint::hy(0);//hx*hy;
	const gmg_layout& lay1 = *(RHS0.lay);
	dof_vector<double>& UVP = *(mgl->common.UVP);

	bool filled = true;
	matrix_type* A;
	multigrid_coarse_level* mgcl = dynamic_cast<multigrid_coarse_level*>(mgl);
	if(mgcl)
	{
		filled = mgcl->filled;
		if(!filled) mgcl->setup();
		A = mgcl->A;
	}
	if(!resid) std::fill(UVP.U.begin()+UVP.offset(level),UVP.U.begin()+UVP.offset(level+1),0.0);

	pint dofs1[5], dofs2[5], dofs3[2];
	int ndofs;
	// equations for Uh
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(1,i,j,level);
		int k = lay1.dof(pdof, 0);
		laplace(pdof,dofs1,dofs2, -mu_s);
		if(level == 0 && !resid) RHS0(0,pdof) += rhs_exact(test, 1,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() )
			{
				if(!filled && !resid)
				{
#if defined(USE_S3M)
					A->PushBack(lay1.dof(dofs1[q],dofs1[q].grid-1), dofs1[q].coef * scale);
#else
					(*A)[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
#endif
				}
				else if(!mgcl || resid)
					RHS0(0,pdof) -= dofs1[q].coef*scale * UVP(0,dofs1[q]);
			}
			else if(level == 0) //if(test == 0)
				RHS0(0,pdof) -= dofs1[q].coef * scale * exact(test, 1,dofs1[q].x(),dofs1[q].y());
		if(divgrad) // \Delta_h = 2 div_h D_h instead of five-point Laplace:
			for(int q = 0; q < 4; ++q)
				if( dofs2[q].inside() )
				{
					if(!filled && !resid)
					{
#if defined(USE_S3M)
						A->PushBack(lay1.dof(dofs2[q],dofs2[q].grid-1), dofs2[q].coef * scale);
#else
						(*A)[k][lay1.dof(dofs2[q],dofs2[q].grid-1)] += dofs2[q].coef * scale;
#endif
					}
					else if(!mgcl || resid)
						RHS0(0,pdof) -= dofs2[q].coef*scale * UVP(1,dofs2[q]);
				}
				else if(level == 0) //if(test == 0)
					RHS0(0,pdof) -= dofs2[q].coef * scale * exact(test, 2,dofs2[q].x(),dofs2[q].y());
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (i != 0 && i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
			{
				if(!filled && !resid)
				{
#if defined(USE_S3M)
					A->PushBack(lay1.dof(dofs3[q],dofs3[q].grid-1), dofs3[q].coef * scale);
#else
					(*A)[k][lay1.dof(dofs3[q],dofs3[q].grid-1)] += dofs3[q].coef * scale;
#endif
				}
				else if(!mgcl || resid)
					RHS0(0,pdof) -= dofs3[q].coef*scale * UVP(2,dofs3[q]);
			}
#if defined(USE_S3M)
		if(!filled) A->FinalizeRow();
#endif
	}
	// equations for Vh
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		pint pdof(2,i,j,level);
		int k = lay1.dof(pdof, 1);
		laplace(pdof,dofs1,dofs2, -mu_s);
		if(level == 0 && !resid) RHS0(1,pdof) += rhs_exact(test, 2,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs2[q].inside() )
			{
				if(!filled && !resid)
				{
#if defined(USE_S3M)
					A->PushBack(lay1.dof(dofs2[q],dofs2[q].grid-1), dofs2[q].coef * scale);
#else
					(*A)[k][lay1.dof(dofs2[q],dofs2[q].grid-1)] += dofs2[q].coef * scale;
#endif
				}
				else if(!mgcl || resid)
					RHS0(1,pdof) -= dofs2[q].coef*scale * UVP(1,dofs2[q]);
			}
			else if(level == 0) //if(test == 0)
				RHS0(1,pdof) -= dofs2[q].coef * scale * exact(test, 2,dofs2[q].x(),dofs2[q].y());
		if(divgrad) // \Delta_h = 2 div_h D_h instead of five-point Laplace:
			for(int q = 0; q < 4; ++q)
				if( dofs1[q].inside() )
				{
					if(!filled && !resid)
					{
#if defined(USE_S3M)
						A->PushBack(lay1.dof(dofs1[q],dofs1[q].grid-1), dofs1[q].coef * scale);
#else
						(*A)[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
#endif
					}
					else if(!mgcl || resid)
						RHS0(1,pdof) -= dofs1[q].coef*scale * UVP(0,dofs1[q]);
				}
				else if(level == 0) //if(test == 0)
					RHS0(1,pdof) -= dofs1[q].coef * scale * exact(test, 1,dofs1[q].x(),dofs1[q].y());
		ndofs = grad(pdof,dofs3);
		if(test != 0 || (j != 0 && j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
			{
				if(!filled && !resid)
				{
#if defined(USE_S3M)
					A->PushBack(lay1.dof(dofs3[q],dofs3[q].grid-1), dofs3[q].coef * scale);
#else
					(*A)[k][lay1.dof(dofs3[q],dofs3[q].grid-1)] += dofs3[q].coef * scale;
#endif
				}
				else if(!mgcl || resid)
					RHS0(1,pdof) -= dofs3[q].coef*scale * UVP(2,dofs3[q]);
			}
#if defined(USE_S3M)
		if(!filled) A->FinalizeRow();
#endif
	}
	// div u = 0, equations for Ph
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		int k = lay1.dof(pdof, 2);
		ndofs = div(pdof,dofs1);
		if(mgcl && !filled && mgcl->common.integral_constraint)
		{
			int nstokes = lay1.size(level);
			(*A)[k][nstokes] = 1;
			(*A)[nstokes][k] = 1;
		}
		for(int q = 0; q < ndofs; ++q)
			if( dofs1[q].inside() )
			{
				if(!filled && !resid)
				{
#if defined(USE_S3M)
					A->PushBack(lay1.dof(dofs1[q],dofs1[q].grid-1), dofs1[q].coef * scale);
#else
					(*A)[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
#endif
				}
				else if(!mgcl || resid)
					RHS0(2,pdof) -= dofs1[q].coef*scale * UVP(dofs1[q].grid-1,dofs1[q]);
			}
			else if(level == 0) //if(test == 0)
				RHS0(2,pdof) -= dofs1[q].coef * scale * exact(test, q < 2 ? 1 : 2,dofs1[q].x(),dofs1[q].y());
#if defined(USE_S3M)
		if(!filled) A->FinalizeRow();
#endif
	}
	if(!filled)
	{
#if defined(USE_S3M)
		mgcl->A->SortRows();
		mgcl->S.Setup(*A);
#else
		mgcl->S.SetMatrix(*A);
#endif
		mgcl->A->Save("A.mtx");
		mgcl->filled = true;
	}
}

void fill_stokes(multigrid_level* mgl, bool resid = false)
{
	fill_stokes(*(mgl->common.RHS0), mgl, resid);
}
double multigrid::residual()
{
	//fill_stokes(head, true);
	//return l2norm(*(common.RHS0), 0);
	fill_stokes(res, head, true);
	return res.l2norm(0);
}

void multigrid_coarse_level::solve()
{
	dof_vector<double>& UVP = *(common.UVP), &RHS0 = *(common.RHS0);
	fill_stokes(this, level == 0);
	//fill_stokes(RHS0, this, level == 0);
	fill_vector(RHS0, *b, level);
#if !defined(USE_S3M)
	std::fill(x->Begin(),x->End(),0.0);
#else
	std::fill(x->begin(),x->end(),0.0);
#endif
	bool success = S.Solve(*b,*x);
	fill_vector(*x, UVP, level);
#if !defined(USE_S3M)
	b->Save("b.txt");
	x->Save("x.txt");
	int iters = S.Iterations();
	double resid = S.Residual();
	std::string reason = S.ReturnReason();
	std::cout << (!success ? "not " : "") << "converged iters " << iters << " resid " << resid << " reason " << reason << std::endl;
#else
	double resid = Resid(*A,*b,*x);
	std::cout << (!success ? "not " : "") << "converged resid " << resid << std::endl;
#endif
}
void multigrid_level::solve()
{
	dof_vector<double>& UVP = *(common.UVP), &RHS0 = *(common.RHS0);
	dof_vector<double> res(common.UVP->lay);

	common.presmooth_iters = common.postsmooth_iters = 10;
	//smoothing(common.presmooth_iters);
	save_vtk("tmp1.vtk", UVP, *(common.lay2), level);
	fill_stokes(this, true); // RHS0 <- RHS0 - A*UVP
	//fill_stokes(this, level == 0); // RHS0 <- RHS0 - A*UVP
	save_vtk("tmp2.vtk", RHS0, *(common.lay2), level);
	restriction(RHS0);
	save_vtk("tmp3.vtk", RHS0, *(common.lay2), level+1);
	if(level < common.nlevels-1) std::fill(UVP.U.begin()+UVP.offset(level+1), UVP.U.begin()+UVP.offset(level+2), 0.0);
	if(1)
	{
		res.assign(RHS0, level+1);
		fill_stokes(res, next, false);
		double resid = res.l2norm(level+1);
		std::cout << "level " << (level+1) << " restricted residual " << resid << "\n" << std::endl;
	}
	for(int q = 0; q < common.schedule; ++q)
	{
		next->solve();
		save_vtk("tmp4.vtk", UVP, *(common.lay2), level+1);
		next->prolongation(*(common.UVP));
	}
	if(1)
	{
		res.assign(RHS0, level);
		//std::fill(res.U.begin()+res.offset(level), res.U.begin()+res.offset(level+1), 0.0);
		fill_stokes(res, this, true);
		double resid = res.l2norm(level);
		std::cout << "level " << level << " prolongated residual " << resid << "\n" << std::endl;
		save_vtk("tmp5.vtk", res, *(common.lay2), level);
	}
	save_vtk("tmp6.vtk", UVP, *(common.lay2), level);
	//smoothing(common.postsmooth_iters);
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
	dof_vector<double> UVP(&lay1), UVP0(&lay3), RHS(&lay1), RHS0(&lay3);
	dof_vector<symmat2> PSI(&lay2), PSI0(&lay2), PSI1(&lay2), TAU(&lay2), TAU0(&lay2), TAU1(&lay2);

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
	save_vtk("tmp1.vtk", PSI, TAU, UVP);
	UVP0.assign(UVP, 0);

	multigrid_params common; common.UVP = &UVP;
	multigrid gmg(common);
	gmg.head->restriction(UVP);
	save_vtk("tmp2.vtk", PSI, TAU, UVP, 1);
	for(int q = UVP.offset(1); q < UVP.offset(2); ++q) UVP.U[q] = 0.0;
	gmg.head->next->prolongation(UVP);

	for(int q = 0; q < UVP.size(); ++q)
		if(fabs(UVP.U[q] - UVP0.U[q]) > 1.0e-12)
			std::cerr << "component " << q
				<< " difference " << fabs(UVP.U[q] - UVP0.U[q])
				<< std::endl;
	
	save_vtk("tmp3.vtk", PSI, TAU, UVP);

	return 0;
}

#endif //_GMG_H
