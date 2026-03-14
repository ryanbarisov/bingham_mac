#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>
//#include <petscksp.h>
#include <inmost.h>

using namespace INMOST;

double r = 1.0; // stabilization parameter
int n1;
int n2;
double hx;
double hy;
int test = 0; // 0 -- regularized cavity, 1,2 -- analytical tests for Stokes problem alone
double R1 = 4.2985;
double R2 = 0.1;
// Bingham fluid
double mu = 1.0;
double tau_s = 0.0;

bool divgrad = true; // Delta_h = 2 div_h D_h (true), standard five-point Laplace (false)

// exact solution for analytical tests, including velocity BC: u,v,p
double exact(int grid, double x, double y)
{
	if(test == 1)
	{
		if(grid == 1) // Uh
			return  1.0/(4*M_PI*M_PI) * (1-cos(2*M_PI*x)) * sin(2*M_PI*y);
		else if(grid == 2) // Vh
			return -1.0/(4*M_PI*M_PI) * (1-cos(2*M_PI*y)) * sin(2*M_PI*x);
		else // if(grid == 3) // Ph
		     	return 1.0/M_PI * sin(2*M_PI*x) * sin(2*M_PI*y);
	}
	else if(test == 2)
	{
		double f = (exp(R1*x)-1) / (exp(R1)-1);
		double g = (exp(R2*y)-1) / (exp(R2)-1);
		double df = R1*exp(R1*x) / (exp(R1)-1);
		double dg = R2*exp(R2*y) / (exp(R2)-1);
		if(grid == 1) // Uh
			return  dg/(2*M_PI) * (1-cos(2*M_PI*f)) * sin(2*M_PI*g);
		else if(grid == 2) // Vh
			return  -df/(2*M_PI) * (1-cos(2*M_PI*g)) * sin(2*M_PI*f);
		else // if(grid == 3) // Ph
		    	return df*dg*sin(2*M_PI*f)*sin(2*M_PI*g);
	}
	else // if(test == 0)
	{
		// no analytical solution, it is used to setup velocity BC
		if(grid == 1) // Uh
			return fabs(y-1) < 1.0e-6 ? 16*x*x*(1-x)*(1-x) : 0.0;
		else if(grid == 2) // Vh
			return 0.0;
		else // if(grid == 3) // Ph
			return 0.0;
	}
}

// analytical RHS for Stokes substep
double rhs_exact(int grid, double x, double y)
{
	if(test == 1)
	{
		if(grid == 1) // Uh
			return -r*( 2*cos(2*M_PI*x)-1)*sin(2*M_PI*y) + 2*cos(2*M_PI*x)*sin(2*M_PI*y);
		else //if(grid == 2) // Vh
			return -r*(-2*cos(2*M_PI*y)+1)*sin(2*M_PI*x) + 2*cos(2*M_PI*y)*sin(2*M_PI*x);
	}
	else if(test == 2)
	{
		double f = (exp(R1*x)-1) / (exp(R1)-1);
		double g = (exp(R2*y)-1) / (exp(R2)-1);
		double df = R1*exp(R1*x) / (exp(R1)-1), d2f = R1*df, d3f = R1*d2f;
		double dg = R2*exp(R2*y) / (exp(R2)-1), d2g = R2*dg, d3g = R2*d2g;
		if(grid == 1) // Uh
		{
			double d2u_dx2 = 2*M_PI*dg*sin(2*M_PI*g)*(df*df*cos(2*M_PI*f) + d2f*sin(2*M_PI*f)/(2*M_PI));
			double d2u_dy2 = (1-cos(2*M_PI*f))*( 3*dg*d2g*cos(2*M_PI*g) + (d3g/(2*M_PI) - 2*M_PI*dg*dg*dg)*sin(2*M_PI*g));
			double dp_dx = d2u_dx2;
			return -r*(d2u_dx2 + d2u_dy2) + dp_dx;
		}
		else //if(grid == 2) // Vh
		{
			double d2v_dy2 = -2*M_PI*df*sin(2*M_PI*f)*(dg*dg*cos(2*M_PI*g) + d2g*sin(2*M_PI*g)/(2*M_PI));
			double d2v_dx2 = -(1-cos(2*M_PI*g))*( 3*df*d2f*cos(2*M_PI*f) + (d3f/(2*M_PI) - 2*M_PI*df*df*df)*sin(2*M_PI*f));
			double dp_dy = -d2v_dy2;
			return -r*(d2v_dx2 + d2v_dy2) + dp_dy;
		}
	}
	else // if(test == 0)
	{
		return 0.0;
	}
}

// DOF on staggered grid
typedef struct pint
{
	int i, j;
	double coef;
	int grid; // 1(Uh), 2(Vh), 3(Ph), 4(Qh)
	pint(int grid = -1, int i = -10, int j = -10) : i(i), j(j), grid(grid) {}
	static int nx(int grid) // amount of unknowns in X direction
	{
		if(grid == 1) return n1;
		else if(grid == 2) return n1-1;
		else if(grid == 3) return n1-1;
		else return n1;
	}
	static int ny(int grid) // amount of unknowns in Y direction
	{
		if(grid == 1) return n2-1;
		else if(grid == 2) return n2;
		else if(grid == 3) return n2-1;
		else return n2;
	}
	static int ndofs(int grid) { return nx(grid)*ny(grid); }
	int dof() const
	{
		if(grid == 1) // Uh
			return i*ny(grid) + (j-1);
		if(grid == 2) // Vh
			return j*nx(grid) + (i-1);
		else if(grid == 3) // Ph
			return (i-1)*ny(grid) + (j-1);
		else // if(grid == 4) // Qh
			return i*ny(grid) + j;
	}
	bool inside() const
	{
		int i1 = i, j1 = j;
		if(grid == 1) j1--;
		else if(grid == 2) i1--;
		else if(grid == 3) i1--, j1--;
		return i1 >= 0 && i1 < nx(grid) && j1 >= 0 && j1 < ny(grid);
	}
	double x() const
	{
		if(grid == 1 || grid == 4)
			return std::min(std::max(0.0,(i+0.5)*hx), 1.0);
		else
			return i*hx;
	}
	double y() const
	{
		if(grid == 2 || grid == 4)
			return std::min(std::max(0.0,(j+0.5)*hy), 1.0);
		else
			return j*hy;
	}
} pint;

// controls layout of DOFs in vector
struct layout
{
	int N;
	int * grids;
	int * offsets;
	layout(int N, const int _grids[]) : N(N)
	{
		grids = new int[N];
		offsets = new int[N+1];
		memcpy(grids,_grids,sizeof(int)*N);
		offsets[0] = 0;
		for(int q = 1; q < N+1; ++q) offsets[q] = offsets[q-1] + pint::ndofs(grids[q-1]);
	}
	layout(const layout& other)
	{
		N = other.N;
		memcpy(grids,other.grids,sizeof(int)*N);
		memcpy(offsets,other.offsets,sizeof(int)*(N+1));
	}
	~layout()
	{
		delete[] grids;
		delete[] offsets;
	}
	int offset(int pos) const
	{
		assert(pos < N+1);
		return offsets[pos];
	}
	int dof(const pint& p, int pos) const
	{
		if(p.inside()) return p.dof() + offset(pos);
		throw 1;
		return -1;
	}
};

void laplace(int i, int j, int grid, pint * dofs1, pint * dofs2, double mult)
{
	double coef1 = divgrad ? 2.0 : 1.0;
	if(grid == 1) // \Delta_h, x component
	{
		// Uh
		double h[2] = {hx,hx}, coef[2];
		dofs1[0] = pint(1,i+1,j);
		dofs1[1] = pint(1,i-1,j);
		coef[0] = dofs1[0].inside() ? 1.0 : 0.5;
		coef[1] = dofs1[1].inside() ? 1.0 : 0.5;
		dofs1[0].coef = coef1/(h[0]*h[0]*coef[0]) * mult;
		dofs1[1].coef = coef1/(h[1]*h[1]*coef[1]) * mult;
		dofs1[2] = pint(1,i,j+1); dofs1[2].coef = 1.0/(hx*hy) * mult;
		dofs1[3] = pint(1,i,j-1); dofs1[3].coef = 1.0/(hx*hy) * mult;
		dofs1[4] = pint(1,i,j  ); dofs1[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs1[4].coef -= dofs1[k].coef;
		// Vh
		dofs2[0] = pint(2,i+1,j  ); dofs2[0].coef =  1.0/(hy*hy) * mult;
		dofs2[1] = pint(2,i+1,j-1); dofs2[1].coef = -1.0/(hy*hy) * mult;
		dofs2[2] = pint(2,i,j    ); dofs2[2].coef = -1.0/(hy*hy) * mult;
		dofs2[3] = pint(2,i,j-1  ); dofs2[3].coef =  1.0/(hy*hy) * mult;
	}
	else // \Delta_h, y component
	{
		// Vh
		double h[2] = {hy,hy}, coef[2];
		dofs2[0] = pint(2,i,j+1);
		dofs2[1] = pint(2,i,j-1);
		coef[0] = dofs2[0].inside() ? 1.0 : 0.5;
		coef[1] = dofs2[1].inside() ? 1.0 : 0.5;
		dofs2[0].coef = coef1/(h[0]*h[0]*coef[0]) * mult;
		dofs2[1].coef = coef1/(h[1]*h[1]*coef[1]) * mult;
		dofs2[2] = pint(2,i+1,j); dofs2[2].coef = 1.0/(hx*hx) * mult;
		dofs2[3] = pint(2,i-1,j); dofs2[3].coef = 1.0/(hx*hx) * mult;
		dofs2[4] = pint(2,i,j  ); dofs2[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs2[4].coef -= dofs2[k].coef;
		// Uh
		dofs1[0] = pint(1,i,  j+1); dofs1[0].coef =  1.0/(hx*hy) * mult;
		dofs1[1] = pint(1,i-1,j+1); dofs1[1].coef = -1.0/(hx*hy) * mult;
		dofs1[2] = pint(1,i,j    ); dofs1[2].coef = -1.0/(hx*hy) * mult;
		dofs1[3] = pint(1,i-1,j  ); dofs1[3].coef =  1.0/(hx*hy) * mult;
	}
}

int div(int i, int j, pint * dofs1)
{
	// Uh
	dofs1[0] = pint(1,i,j  ); dofs1[0].coef =  1.0 / hx;
	dofs1[1] = pint(1,i-1,j); dofs1[1].coef = -1.0 / hx;
	// Vh
	dofs1[2] = pint(2,i,j  ); dofs1[2].coef =  1.0 / hy;
	dofs1[3] = pint(2,i,j-1); dofs1[3].coef = -1.0 / hy;
	return 4;
}

int grad(int i, int j, int dir, pint * dofs3)
{
	// Ph
	if(dir == 1) // x
	{
		dofs3[0] = pint(3,i+1,j); dofs3[0].coef =  1.0 / hx;
		dofs3[1] = pint(3,i,j  ); dofs3[1].coef = -1.0 / hx;
	}
	else //if(dir == 2) // y
	{
		dofs3[0] = pint(3,i,j+1); dofs3[0].coef =  1.0 / hy;
		dofs3[1] = pint(3,i,j  ); dofs3[1].coef = -1.0 / hy;
	}
	return 2;
}

// compute implicitly strain rate tensor
int Dh(int grid, int comp, int i, int j, pint * dofs)
{
	int n = 0, m = 0, b;
	double h;
	if(grid == 3)
	{
		if(comp == 0) // xx
		{
			n = 2;
			dofs[0] = pint(1,i,j);
			dofs[1] = pint(1,i-1,j);
			dofs[0].coef =  1.0/hx;
			dofs[1].coef = -1.0/hx;
			for(int q = 0; q < 2; ++q)
				if(!dofs[q].inside() ) dofs[q].coef *= exact(1,dofs[q].x(),dofs[q].y());
		}
		else if(comp == 1) // yy
		{
			n = 2;
			dofs[0] = pint(2,i,j);
			dofs[1] = pint(2,i,j-1);
			dofs[0].coef =  1.0/hy;
			dofs[1].coef = -1.0/hy;
			for(int q = 0; q < 2; ++q)
				if(!dofs[q].inside() ) dofs[q].coef *= exact(2,dofs[q].x(),dofs[q].y());
		}
		else //if(comp == 2) // xy
		{
			// averaging
			if((b = Dh(4,comp,i-1,j-1,dofs+n))) ++m, n += b;
			if((b = Dh(4,comp,i-1,j  ,dofs+n))) ++m, n += b;
			if((b = Dh(4,comp,i  ,j-1,dofs+n))) ++m, n += b;
			if((b = Dh(4,comp,i  ,j  ,dofs+n))) ++m, n += b;
			if(m) for(int q = 0; q < n; ++q)
				dofs[q].coef *= 1.0/m;
			else std::cerr << __LINE__ << std::endl;
		}
	}
	else //if(grid == 4)
	{
		if(comp == 2) // xy
		{
			n = 4;
			dofs[0] = pint(1,i,j+1);
			dofs[1] = pint(1,i,j);
			dofs[2] = pint(2,i+1,j);
			dofs[3] = pint(2,i,j);
			dofs[0].coef =  0.5/hy;
			dofs[1].coef = -0.5/hy;
			dofs[2].coef =  0.5/hx;
			dofs[3].coef = -0.5/hx;
			for(int q = 0; q < n; ++q)
				if(!dofs[q].inside() ) dofs[q].coef *= exact(q < 2 ? 1 : 2, dofs[q].x(),dofs[q].y());
		}
		else //if(comp == 0 || comp == 1) // xx, yy
		{
			// averaging
			if((b = Dh(3,comp,i+1,j+1,dofs+n))) ++m, n += b;
			if((b = Dh(3,comp,i+1,j  ,dofs+n))) ++m, n += b;
			if((b = Dh(3,comp,i  ,j+1,dofs+n))) ++m, n += b;
			if((b = Dh(3,comp,i  ,j  ,dofs+n))) ++m, n += b;
			if(m) for(int q = 0; q < n; ++q)
				dofs[q].coef *= 1.0/m;
			else std::cerr << __LINE__ << std::endl;
		}
	}
	return n;
}
// compute implicitly tau with averaging where necessary
int tauh(int grid, int comp, int i, int j, pint * dofs)
{
	if(grid == 3)
	{
		if(comp == 0 || comp == 1) // xx, yy
		{
			dofs[0] = pint(3,i,j); dofs[0].coef = 1.0;
			return 1;
		}
		else //if(comp == 2) // xy
		{
			// averaging
			int m = 0;
			dofs[m] = pint(4,i-1,j-1); if( dofs[m].inside() ) ++m;
			dofs[m] = pint(4,i-1,j  ); if( dofs[m].inside() ) ++m;
			dofs[m] = pint(4,i  ,j-1); if( dofs[m].inside() ) ++m;
			dofs[m] = pint(4,i  ,j  ); if( dofs[m].inside() ) ++m;
			for(int q = 0; q < m; ++q) dofs[q].coef = m ? 1.0/m : 0.0;
			return m;
		}
	}
	else //if(grid == 4)
	{
		if(comp == 2) // xy
		{
			dofs[0] = pint(4,i,j); dofs[0].coef = 1.0;
			return 1;
		}
		else //if(comp == 0 || comp == 1) // xx, yy
		{
			// averaging
			int m = 0;
			dofs[m] = pint(3,i+1,j+1); if( dofs[m].inside() ) ++m;
			dofs[m] = pint(3,i+1,j  ); if( dofs[m].inside() ) ++m;
			dofs[m] = pint(3,i  ,j+1); if( dofs[m].inside() ) ++m;
			dofs[m] = pint(3,i  ,j  ); if( dofs[m].inside() ) ++m;
			for(int q = 0; q < m; ++q) dofs[q].coef = m ? 1.0/m : 0.0;
			return m;
		}
	}
}

// Frobenius norm of symmetric 2x2 matrix
inline double frobnorm(double gamma[3])
{
	return sqrt(gamma[0]*gamma[0]+gamma[1]*gamma[1]+2*gamma[2]*gamma[2]);
}

void update_gamma(const layout& lay1, const layout& lay2, double * gamma, const double * tau, const Sparse::Vector& x)
{
	pint dofs[32];
	int ndofs;
	// grid 3, diagonal components
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		double tau_ij[3] = {0,0,0}; //xx,yy,xy
		for(int comp = 0; comp < 3; ++comp)
		{
			ndofs = Dh(3,comp,i,j,dofs);
			for(int q = 0; q < ndofs; ++q)
				tau_ij[comp] += 2*r*dofs[q].coef*( dofs[q].inside() ? x[lay1.dof(dofs[q],dofs[q].grid-1)] : 1.0 );
			ndofs = tauh(3,comp,i,j,dofs);
			for(int q = 0; q < ndofs; ++q) if( dofs[q].inside() )
				tau_ij[comp] += dofs[q].coef*( tau[lay2.dof(dofs[q],comp)] );
		}
		double fnorm = frobnorm(tau_ij);
		for(int comp = 0; comp < 2; ++comp)
			gamma[lay2.dof(pint(3,i,j),comp)] = fnorm < tau_s ? 0.0 : (1-tau_s/(fnorm+1.0e-20)) / (2*(r+mu)) * tau_ij[comp];
	}
	// grid 4, mixed component
	for(int i = 0; i < n1; ++i)
		for(int j = 0; j < n2; ++j)
	{
		double tau_ij[3] = {0,0,0}; //xx,yy,xy
		for(int comp = 0; comp < 3; ++comp)
		{
			ndofs = Dh(4,comp,i,j,dofs);
			for(int q = 0; q < ndofs; ++q)
				tau_ij[comp] += 2*r*dofs[q].coef*( dofs[q].inside() ? x[lay1.dof(dofs[q],dofs[q].grid-1)] : 1.0 );
			ndofs = tauh(4,comp,i,j,dofs);
			for(int q = 0; q < ndofs; ++q) if( dofs[q].inside() )
				tau_ij[comp] += dofs[q].coef*( tau[lay2.dof(dofs[q],comp)] );
		}
		double fnorm = frobnorm(tau_ij);
		gamma[lay2.dof(pint(4,i,j),2)] = fnorm < tau_s ? 0.0 : (1-tau_s/(fnorm+1.0e-20)) / (2*(r+mu)) * tau_ij[2];
	}
}

void update_tau(const layout& lay1, const layout& lay2, const double * gamma, const double * tau0, double * tau, const Sparse::Vector& x)
{
	pint dofs[32];
	int ndofs;
	//const double rn = r/2;//(1+sqrt(5))/2*r;
	const double rn = (1+sqrt(5))/2*r;
	// grid 3, diagonal components
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		for(int comp = 0; comp < 2; ++comp)
		{
			int m = lay2.dof(pint(3,i,j),comp);
			tau[m] = tau0[m] - 2*rn*gamma[m];
			ndofs = Dh(3,comp,i,j,dofs);
			for(int q = 0; q < ndofs; ++q)
				tau[m] += 2*rn*dofs[q].coef * ( dofs[q].inside() ? x[lay1.dof(dofs[q],dofs[q].grid-1)] : 1.0 );
		}
	}
	// grid 4, mixed component
	for(int i = 0; i < n1; ++i)
		for(int j = 0; j < n2; ++j)
	{
		int m = lay2.dof(pint(4,i,j),2);
		tau[m] = tau0[m] - 2*rn*gamma[m];
		ndofs = Dh(4,2,i,j,dofs);
		for(int q = 0; q < ndofs; ++q)
			tau[m] += 2*rn*dofs[q].coef * ( dofs[q].inside() ? x[lay1.dof(dofs[q],dofs[q].grid-1)] : 1.0 );
	}
}

double l2norm(double * vec, int n)
{
	double norm = 0.0;
	for(int k = 0; k < n; ++k)
		norm += vec[k]*vec[k];
	return sqrt(norm);
}

double residual(const double* tau0, const double* tau)
{
	double l2norm = 0.0;
	for(int k = 0; k < 2*(n1-1)*(n2-1)+n1*n2; ++k)
	{
		double err = fabs(tau0[k]-tau[k]);
		l2norm += err*err*hx*hy;
	}
	return sqrt(l2norm);
}

int div_tau(int comp, int i, int j, pint * dofs)
{
	double h;
	if(comp == 1) // x
	{
		// xx
		dofs[0] = pint(3,i+1,j);
		dofs[1] = pint(3,i  ,j);
		dofs[0].coef =  1.0 / hx;
		dofs[1].coef = -1.0 / hx;
		// xy
		dofs[2] = pint(4,i,j  );
		dofs[3] = pint(4,i,j-1);
		dofs[2].coef =  1.0 / hx;
		dofs[3].coef = -1.0 / hx;
	}
	else //if(comp == 2) // y
	{
		// yy
		dofs[0] = pint(3,i,j+1);
		dofs[1] = pint(3,i,j  );
		dofs[0].coef =  1.0 / hy;
		dofs[1].coef = -1.0 / hy;
		// xy
		dofs[2] = pint(4,i  ,j);
		dofs[3] = pint(4,i-1,j);
		dofs[2].coef =  1.0 / hy;
		dofs[3].coef = -1.0 / hy;
	}
	return 4;
}

void fill_rhs(const layout& lay1, const layout& lay2, const double * gamma, const double * tau, Sparse::Vector& b, double scale = 1.0)
{
	pint dofs[32];
	int ndofs;
	// grid1
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		int m = lay1.dof(pint(1,i,j),0);
		ndofs = div_tau(1,i,j,dofs);
		for(int q = 0; q < ndofs/2; ++q) if( dofs[2*q].inside() && dofs[2*q+1].inside() )
		for(int k = 0; k < 2; ++k)
		{
			int d = lay2.dof(dofs[2*q+k],2*q+k < 2 ? 0 : 2);
			b[m] += dofs[2*q+k].coef * (tau[d] - 2*r*gamma[d]) * scale;
		}
}
	// grid2
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		int m = lay1.dof(pint(2,i,j),1);
		ndofs = div_tau(2,i,j,dofs);
		for(int q = 0; q < ndofs/2; ++q) if( dofs[2*q].inside() && dofs[2*q+1].inside() )
		for(int k = 0; k < 2; ++k)
		{
			int d = lay2.dof(dofs[2*q+k],2*q+k < 2 ? 1 : 2);
			b[m] += dofs[2*q+k].coef * (tau[d] - 2*r*gamma[d]) * scale;
		}
	}
}

void save_vtk(std::string filename, const layout& lay1, const layout& lay2, const double * gamma, const double * tau, const Sparse::Vector& x)
{
	std::ofstream ofs(filename);

	//version,identifier
	ofs << "# vtk DataFile Version 3.0\n";
	//header
	ofs << "Bingham flow rectangular grid " << n1 << " " << n2 << '\n';
	//file format
	ofs << "ASCII\n";
	//dataset type
	ofs << "DATASET RECTILINEAR_GRID\n";
	ofs << "DIMENSIONS " << (n1+1) << " " << (n2+1) << " 1" << '\n';
	ofs << "X_COORDINATES " << (n1+1) << " double\n";
	for(int i = 0; i < n1+1; ++i)
		ofs << i*hx << ' ';
	ofs << '\n';
	ofs << "Y_COORDINATES " << (n1+1) << " double\n";
	for(int j = 0; j < n2+1; ++j)
		ofs << j*hy << ' ';
	ofs << '\n';
	ofs << "Z_COORDINATES " << 1 << " double\n";
	ofs << "0.0\n";

	//data
	ofs << "POINT_DATA " << ((n1+1)*(n2+1)) << '\n';
	//pressure
	ofs << "SCALARS pressure double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 0; j < n2+1; ++j)
		for(int i = 0; i < n1+1; ++i)
	{
		double p = 0.0;
		if(i > 0 && i < n1 && j > 0 && j < n2)
			p = x[lay1.dof(pint(3,i,j),2)];
		ofs << p << '\n';
	}
	//divergence
	ofs << "SCALARS div_U double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 0; j < n2+1; ++j)
		for(int i = 0; i < n1+1; ++i)
	{
		double divu = 0.0;
		pint dofs[8];
		if(i > 0 && i < n1 && j > 0 && j < n2)
		{
			int ndofs = div(i,j,dofs);
			for(int q = 0; q < ndofs; ++q)
				if( dofs[q].inside() )
					divu += x[lay1.dof(dofs[q],dofs[q].grid-1)] * dofs[q].coef;
				else
					divu += exact(q < 2 ? 1 : 2, dofs[q].x(), dofs[q].y()) * dofs[q].coef;
		}
		ofs << divu << '\n';
	}
	//frobenius norm of tau
	ofs << "SCALARS frob_tau double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 0; j < n2+1; ++j)
		for(int i = 0; i < n1+1; ++i)
	{
		pint dofs[32];
		double tau_ij[3] = {0,0,0}; //xx,yy,xy
		if(i > 0 && i < n1 && j > 0 && j < n2)
		{
			for(int comp = 0; comp < 3; ++comp)
			{
				int ndofs = tauh(3,comp,i,j,dofs);
				for(int q = 0; q < ndofs; ++q) if( dofs[q].inside() )
					tau_ij[comp] += dofs[q].coef*( tau[lay2.dof(dofs[q],comp)] );
			}
		}
		ofs << frobnorm(tau_ij) << '\n';
	}
	// u
	ofs << "SCALARS U double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 0; j < n2+1; ++j)
		for(int i = 0; i < n1+1; ++i)
	{
		double u[2] = {0,0};
		pint dofs[4];
		{
			dofs[0] = pint(1,i,j  );
			dofs[1] = pint(1,i-1,j);
			for(int k = 0; k < 2; ++k)
				if(dofs[k].inside()) u[k] = dofs[k].coef = x[lay1.dof(dofs[k], 0)];
				else u[0] = u[1] = dofs[k].coef = exact(1, dofs[k].x(), dofs[k].y());
		}
		ofs << 0.5*(u[0]+u[1]) << '\n';
	}
	// v
	ofs << "SCALARS V double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 0; j < n2+1; ++j)
		for(int i = 0; i < n1+1; ++i)
	{
		double u[2] = {0,0};
		pint dofs[4];
		{
			dofs[0] = pint(2,i,j  );
			dofs[1] = pint(2,i,j-1);
			for(int k = 0; k < 2; ++k)
				if(dofs[k].inside()) u[k] = dofs[k].coef = x[lay1.dof(dofs[k], 1)];
				else u[0] = u[1] = dofs[k].coef = exact(2, dofs[k].x(), dofs[k].y());
		}
		ofs << 0.5*(u[0]+u[1]) << '\n';
	}

	ofs.close();
}

int main(int argc, char ** argv)
{
	Solver::Initialize(&argc, &argv, "database.xml");

	if(argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " n1 [test=" << test
			<< "] [r=" << r << "] [tau_s=" << tau_s
			<< "] [R1=" << R1 << "] [R2=" << R2
			<< '\n';
		return 0;
	}
	
	n1 = atoi(argv[1]);
	n2 = n1;//argv[2];
	hx = 1.0/n1;
	hy = 1.0/n2;
	if(argc > 2)
		test = atoi(argv[2]);
	if(argc > 3)
		r = atof(argv[3]);
	if(argc > 4)
		tau_s = atof(argv[4]);
	if(argc > 5)
		R1 = atof(argv[5]);
	if(argc > 6)
		R2 = atof(argv[6]);
	
	const int offset1[3] = {1,2,3};
	const int offset2[3] = {3,3,4};
	layout lay1(3,offset1), lay2(3,offset2);
	std::cout << "n1 " << n1 << " n2 " << n2 << " test " << test << " r " << r << " tau_s " << tau_s  << " (tau_s for test=0 only) R1 " << R1 << " R2 " << R2 << " (R1,R2 only for test=2)" << std::endl;

	int ntau = lay2.offset(3);//2*(n1-1)*(n2-1) + n1*n2;
	double * gamma = new double[ntau]; // gamma_11 (Ph), gamma_22 (Ph), gamma_12 (Qh)
	double * tau = new double[ntau]; // tau_11 (Ph), tau_22 (Ph), tau_12 (Qh)
	double * tau0 = new double[ntau]; // prev step
	memset(gamma, 0, ntau*sizeof(double));
	memset(tau0, 0, ntau*sizeof(double));
	memset(tau, 0, ntau*sizeof(double));

	int nstokes = lay1.offset(3);//n1*(n2-1)+n2*(n1-1)+(n1-1)*(n2-1); // Uh,Vh,Ph

	//std::cout << "check ntau1 " << ntau << " ntau2 " << 2*(n1-1)*(n2-1) + n1*n2 << std::endl;
	//std::cout << "check nstokes1 " << nstokes << " nstokes2 " << n1*(n2-1)+n2*(n1-1)+(n1-1)*(n2-1) << std::endl;

	Sparse::Matrix A("stokes", 0, nstokes);
	Sparse::Vector b("stokes", 0, nstokes), x("stokes", 0, nstokes), b0("stokes", 0, nstokes), x0("stokes", 0, nstokes);
	std::fill(b0.Begin(),b0.End(),0.0);
	// equations for velocity
	pint dofs1[5], dofs2[5], dofs3[2];
	int ndofs;
	double scale = hx*hy;
	// equations for Uh
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		int k = lay1.dof(pint(1,i,j), 0);
		laplace(i,j,1,dofs1,dofs2, -r);
		b0[k] = rhs_exact(1,(i+0.5)*hx,j*hy) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() )
				A[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
			else //if(test == 0)
				b0[k] -= dofs1[q].coef * scale * exact(1,dofs1[q].x(),dofs1[q].y());
		if(divgrad) // \Delta_h = 2 div_h D_h instead of five-point Laplace:
			for(int q = 0; q < 4; ++q)
				if( dofs2[q].inside() )
					A[k][lay1.dof(dofs2[q],dofs2[q].grid-1)] += dofs2[q].coef * scale;
				else //if(test == 0)
					b0[k] -= dofs2[q].coef * scale * exact(2,dofs2[q].x(),dofs2[q].y());
		ndofs = grad(i,j,1,dofs3);
		if(i != 0 && i != n1-1) //skipping boundary gives BC dp/dx=0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
				A[k][lay1.dof(dofs3[q],dofs3[q].grid-1)] += dofs3[q].coef * scale;
	}
	// equations for Vh
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		int k = lay1.dof(pint(2,i,j), 1);
		laplace(i,j,2,dofs1,dofs2, -r);
		b0[k] = rhs_exact(2,i*hx,(j+0.5)*hy) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs2[q].inside() )
				A[k][lay1.dof(dofs2[q],dofs2[q].grid-1)] += dofs2[q].coef * scale;
			else //if(test == 0)
				b0[k] -= dofs2[q].coef * scale * exact(2,dofs2[q].x(),dofs2[q].y());
		if(divgrad) // \Delta_h = 2 div_h D_h instead of five-point Laplace:
			for(int q = 0; q < 4; ++q)
				if( dofs1[q].inside() )
					A[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
				else //if(test == 0)
					b0[k] -= dofs1[q].coef * scale * exact(1,dofs1[q].x(),dofs1[q].y());
		ndofs = grad(i,j,2,dofs3);
		if(j != 0 && j != n2-1)//skipping boundary gives BC dp/dy=0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
				A[k][lay1.dof(dofs3[q],dofs3[q].grid-1)] += dofs3[q].coef * scale;
	}
	// div u = 0, equations for Ph
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		int k = lay1.dof(pint(3,i,j), 2);
		ndofs = div(i,j,dofs1);
		for(int q = 0; q < ndofs; ++q)
			if( dofs1[q].inside() )
				A[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
			else //if(test == 0)
				b0[k] -= dofs1[q].coef * scale * exact(q < 2 ? 1 : 2,dofs1[q].x(),dofs1[q].y());
	}

	A.Save("A.mtx");
	b0.Save("b0.txt");
	bool success = false, finish = false;
	double res = 0.0, atol = 1.0e-4;
	int iter = 0, maxiters = 2000;
	Solver S(Solver::INNER_ILU2, "test");
	S.SetMatrix(A);
	do
	{
		std::copy(b0.Begin(),b0.End(),b.Begin());
		std::fill(x.Begin(),x.End(),0.0);
		if(test == 0) fill_rhs(lay1,lay2,gamma, tau, b, scale);
		success = S.Solve(b,x);
		int iters = S.Iterations();
		double resid = S.Residual();
		std::string reason = S.ReturnReason();
		std::cout << (!success ? "not " : "") << "converged iters " << iters << " resid " << resid << " reason " << reason << std::endl;
		if(test)
			finish = true;
		else if(success)
		{
			update_gamma(lay1, lay2, gamma, tau0, x);
			update_tau(lay1, lay2, gamma, tau0, tau, x);
			res = residual(tau0, tau);
			std::copy(tau, tau+ntau, tau0);
			std::copy(x.Begin(), x.End(), x0.Begin());
			std::cout << "iter " << iter << " " << " resid " << res << std::endl;
			//std::cin.ignore();
			++iter;
			finish = res < atol;
			if(iter >= maxiters || !success)
			{
				finish = false;
				break;
			}
		}
		//x.Save("x.txt");
		//b.Save("b.txt");
	}
	while(!finish);
	//if(finish)
	{
		// save vtk
		std::copy(x0.Begin(),x0.End(),x.Begin());
		save_vtk("out.vtk", lay1, lay2, gamma, tau, x);
		// save fields in txt to visualize using plot.py
		{
			std::ofstream ofs("x.txt");
			ofs << n1 << ' ' << n2 << '\n';
			for(int q = 0; q < x.Size(); ++q)
				ofs << std::scientific << x[q] << '\n';
			ofs.close();
		}
		std::ofstream ofs1("tau.txt");
		std::ofstream ofs2("gamma.txt");
		pint dofs[32];
		for(int i = 1; i < n1; ++i)
			for(int j = 1; j < n2; ++j)
		{
			double tau_ij[3] = {0,0,0}, gamma_ij[3] = {0,0,0}; //xx,yy,xy
			for(int comp = 0; comp < 3; ++comp)
			{
				ndofs = Dh(3,comp,i,j,dofs);
				for(int q = 0; q < ndofs; ++q)
					gamma_ij[comp] += dofs[q].coef*( dofs[q].inside() ? x[lay1.dof(dofs[q],dofs[q].grid-1)] : 1.0 );
				ndofs = tauh(3,comp,i,j,dofs);
				for(int q = 0; q < ndofs; ++q) if( dofs[q].inside() )
				{
					//gamma_ij[comp] += dofs[q].coef*( gamma[lay2.dof(dofs[q],comp)] );
					tau_ij[comp] += dofs[q].coef*( tau0[lay2.dof(dofs[q],comp)] );
				}
			}
			double fnorm1 = frobnorm(tau_ij)+1.0e-10;
			ofs1 << fnorm1  << " "<< tau_ij[0] << " " << tau_ij[1] << " " << tau_ij[2] << "\n";
			double fnorm2 = frobnorm(gamma_ij)+1.0e-10;
			ofs2 << fnorm2  << " "<< gamma_ij[0] << " " << gamma_ij[1] << " " << gamma_ij[2] << "\n";
		}
		ofs1.close();
		ofs2.close();

		// compute errors for analytical test
		if(test != 0)
		{
			double errUV_L2 = 0.0, errP_L2 = 0.0;
			double errUV_C = 0.0, errP_C = 0.0;
			double xyUV[2], xyP[2];
			double err;
			for(int i = 0; i < n1; ++i)
				for(int j = 1; j < n2; ++j)
			{
				int k = lay1.dof(pint(1,i,j),0);
				b[k] = exact(1,(i+0.5)*hx,j*hy);
				err = fabs( b[k] - x[k] );
				if(errUV_C < err)
				{
					errUV_C = err;
					xyUV[0] = (i+0.5)*hx;
					xyUV[1] = j*hy;
				}
				errUV_L2 += err*err * hx*hy;
			}
			for(int j = 0; j < n2; ++j)
				for(int i = 1; i < n1; ++i)
			{
				int k = lay1.dof(pint(2,i,j),1);
				b[k] = exact(2,i*hx,(j+0.5)*hy);
				err = fabs( b[k] - x[k] );
				if(errUV_C < err)
				{
					errUV_C = err;
					xyUV[0] = i*hx;
					xyUV[1] = (j+0.5)*hy;
				}
				errUV_L2 += err*err * hx*hy;
			}
			for(int i = 1; i < n1; ++i)
				for(int j = 1; j < n2; ++j)
			{
				int k = lay1.dof(pint(3,i,j),2);
				b[k] = exact(3,i*hx,j*hy);
				err = fabs( b[k] - x[k] );
				if(errP_C < err)
				{
					errP_C = err;
					xyP[0] = i*hx;
					xyP[1] = j*hy;
				}
				errP_L2 += err*err * hx*hy;
			}
			errUV_L2 = sqrt(errUV_L2);
			errP_L2 = sqrt(errP_L2);
			std::cout << "|eh| " << errUV_L2 << " |rh| " << errP_L2 << " |eh|_C " << errUV_C << " at (" << xyUV[0] << " " << xyUV[1] << ") |rh|_C " << errP_C << " at (" << xyP[0] << " " << xyP[1] << ")" << std::endl;
			b.Save("exact.txt");
		}
	}

	delete[] gamma;
	delete[] tau;
	delete[] tau0;

	Solver::Finalize();
	return 0;
}
