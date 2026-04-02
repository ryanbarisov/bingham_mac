#ifndef _GMG_H
#define _GMG_H

#include "vtk.h"

#include <inmost.h>
using namespace INMOST;

void fill_vector(const dof_vector<double>& in, Sparse::Vector& out, int level = 0)
{
	std::copy(in.U.begin()+in.offset(level), in.U.begin()+in.offset(level)+in.size(level), out.Begin());
}
void fill_vector(const Sparse::Vector& in, dof_vector<double>& out, int level = 0)
{
	//std::copy(in.Begin(), in.End(), out.U.begin()+out.offset(level));
	for(int k = 0; k < out.size(level); ++k)
		out.U[out.offset(level)+k] = in[k];
	//std::copy(in.Begin(), in.Begin()+out.offset(level), out.U.begin()+out.offset(level));
}

int level0 = 0;

typedef struct multigrid_params
{
	int nlevels = 2;
	int schedule = 1; // 1 -- V, 2 -- W
	int smooth_iters = 4;
	int maxiters = 10;
	double rtol = 1.0e-6, atol = 1.0e-8;

	std::string slv_type = Solver::INNER_ILU2;
	std::string slv_prefix = "test";

	bool integral_constraint = true;
} multigrid_params;

struct multigrid; //fwd declaration

typedef struct multigrid_level
{
	int level;
	multigrid* ptr;
	struct multigrid_level* next = NULL;
	
	multigrid_level(int level, multigrid* ptr);
	multigrid_level(const multigrid_level& other)
		: level(other.level), ptr(other.ptr)
	{
		if(other.next != NULL) next = new multigrid_level(*(other.next));
	}
	multigrid_level& operator=(const multigrid_level& other)
	{
		level = other.level;
		ptr = other.ptr;
		if(other.next != NULL) next = new multigrid_level(*(other.next));
		return *this;
	}
	virtual ~multigrid_level() {if(next != NULL) delete next;}
	virtual void solve();
} multigrid_level;

struct multigrid_coarse_level : multigrid_level
{
	bool filled = false;
	Solver *S = NULL;
	Sparse::Matrix *A = NULL;
	Sparse::Vector *b = NULL, *x = NULL;
	
	multigrid_coarse_level(int level, multigrid* ptr) : multigrid_level(level, ptr) {}
	virtual ~multigrid_coarse_level()
	{
		if(filled)
		{
			delete A;
			delete b;
			delete x;
			delete S;
		}
	}
	virtual void setup();
	virtual void solve();
};

struct multigrid
{
	multigrid_params common;
	multigrid_level* head;
	
	dof_vector<double> UVP, RHS, RHS0;
	dof_vector<double> RESID;

	multigrid(const multigrid_params& common, const gmg_layout& lay);

	virtual void fill_residual(dof_vector<double>& R, multigrid_level* mgl) = 0;
	virtual void smoothing(multigrid_level* mgl) = 0;
	virtual void restriction(dof_vector<double>& R, multigrid_level* mgl) = 0;
	virtual void restriction(multigrid_level* mgl) {restriction(RHS, mgl);}
	virtual void prolongation(multigrid_level* mgl) = 0;

	void step()
	{
		RHS.assign(RHS0, level0);
		head->solve();
	}
	double residual()
	{
		RESID.zero(level0);//RESID <- 0
		fill_residual(RESID, head);//RESID <- b-Ax
		return RESID.l2norm(level0);
	}
	bool solve()
	{
		bool success = false;
		double resid = 0.0, resid0 = 0.0;
		int iter = 0;
		RHS0.assign(RHS, level0); // save assigned RHS
		
		const int grids2[1] = {3};	// Psi(3), Psi=(Psi_xx, Psi_yy, Psi_xy)
		const int ncomps2[1] = {3};
		gmg_layout lay2(1,grids2,ncomps2,1);

		do
		{
			step();
			resid = residual();
			if(iter == 0) resid0 = resid;
			if(resid/resid0 < common.rtol || resid < common.atol) success = true;
			else if(iter > common.maxiters) {success = false; break;}
			//save_vtk("tmp"+std::to_string(iter+1)+".vtk", UVP, lay2);
			//save_vtk_poisson("tmp"+std::to_string(iter+1)+".vtk", UVP);
			++iter;
			std::cout << "gmg iter " << iter << " resid " << resid << " resid/resid0 " << (resid/resid0) << std::endl;
			//success = true;
		} while(!success);
		std::cout << "gmg: solved? " << success << " iters " << iter << std::endl;
		return success;
	}
};

static multigrid_level* create_level(int level, multigrid* ptr)
{
	if(level == ptr->common.nlevels - 1)
		return new multigrid_coarse_level(level, ptr);
	else if(level < ptr->common.nlevels - 1)
		return new multigrid_level(level, ptr);
	else return NULL;
}
multigrid::multigrid(const multigrid_params& common, const gmg_layout& lay) : common(common), UVP(&lay), RHS(&lay), RHS0(&lay), RESID(&lay) {head = create_level(level0,this);}
multigrid_level::multigrid_level(int level, multigrid* ptr) : level(level), ptr(ptr)
{
	next = create_level(level+1, ptr);
}

void multigrid_coarse_level::setup()
{
	if(!filled)
	{
		int n = ptr->UVP.size(level);
		if(ptr->common.integral_constraint) n += 1; // integral constraint on pressure
		S = new Solver(ptr->common.slv_type, ptr->common.slv_prefix);
		A = new Sparse::Matrix("gmg", 0, n);
		b = new Sparse::Vector("gmg", 0, n);
		x = new Sparse::Vector("gmg", 0, n);
	}
}


void multigrid_coarse_level::solve()
{
	ptr->UVP.zero(level); // x <- 0
	ptr->fill_residual(ptr->RHS, this);//fill b-Ax = b into RHS, x=0
	fill_vector(ptr->RHS, *b, level);
	//if(ptr->common.integral_constraint) (*b)[ptr->RHS.size(level)] = 0.0;
	std::fill(x->Begin(),x->End(),0.0);
	bool success = S->Solve(*b,*x);
	fill_vector(*x, ptr->UVP, level);
	b->Save("b.txt");
	x->Save("x.txt");
	int iters = S->Iterations();
	double resid = S->Residual();
	std::string reason = S->ReturnReason();
	std::cout << (!success ? "not " : "") << "converged iters " << iters << " resid " << resid << " reason " << reason << std::endl;
}
void multigrid_level::solve()
{
	ptr->smoothing(this);
	ptr->fill_residual(ptr->RHS, this); // RHS <- RHS - A*UVP
	ptr->restriction(this);// RHS_coarse <- RHS_fine
	int ntimes = level == level0 ? 1 : ptr->common.schedule;
	//ptr->RESID.assign(ptr->RHS, level+1);
	ptr->UVP.zero(level+1); // x_c <- 0
	for(int q = 0; q < ntimes; ++q)
	{
		next->solve();
		//ptr->RHS.assign(ptr->RESID, level+1);
	}
	ptr->prolongation(next);
	if(level == level0) ptr->smoothing(this);
}

#endif //_GMG_H
