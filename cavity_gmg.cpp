#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>
//#include <petscksp.h>
#include <inmost.h>
#include <vector>
#include <array>
#include "analytics.h"

//#define USE_S3M

#if defined(USE_S3M)
#include "s3m/csrmatrix.h"
#include "s3m/bpcg.h"
#include "s3m/pcg.h"
#include "s3m/bicgstab.h"
#include "s3m/amg_ruge_stuben.h"
#include "s3m/gauss_seidel.h"
#include "s3m/ilduc.h"
#include "s3m/vanka.h"
#endif

using namespace INMOST;

//doi: 10.1016/j.jnnfm.2004.12.003

double mu_s = 1.0; // solvent viscosity
double mu_p = mu_s; // polymer viscosity
double We = 1.0;

int glev = 6; // n1 = 2^glev+1
int n1;
int n2;
double hx;
double hy;
double t = 0, T = 1.0;
double dt = 1.0, dt0 = 0.0;
double cfl = 0.5; // Kurganov-Tadmor scheme is stable for cfl <= 0.5
int nframe = 10;

int test = 0; // 0 -- regularized cavity, 1,2 -- analytical tests for Stokes problem alone

bool divgrad = false; // Delta_h = 2 div_h D_h (true), standard five-point Laplace (false)


// DOF on staggered grid
typedef struct pint
{
	int i, j;
	double coef = 0.0;
	int grid; // 1(Uh), 2(Vh), 3(Ph), 4(Qh)
	int level = 0; // multigrid level; level = 0 -- finest grid

	pint(int grid = -1, int i = -10, int j = -10, int level = 0) : i(i), j(j), grid(grid), level(level) {}
	static int nx(int grid, int level) // amount of unknowns in X direction depends on level
	{
		if(grid == 1 || grid == 4) return (1 << (glev-level)) + 1;
		else return 1 << (glev-level) + 0;
	}
	static int ny(int grid, int level) // amount of unknowns in X direction depends on level
	{
		if(grid == 2 || grid == 4) return (1 << (glev-level)) + 1;
		else return 1 << (glev-level) + 0;
	}
	static double hx(int level) { return 1.0/nx(1,level); }
	static double hy(int level) { return 1.0/ny(2,level); }
	static int ndofs(int grid, int level) { return nx(grid,level)*ny(grid,level); }
	int dof() const
	{
		if(grid == 1) // Uh
			return i*ny(grid,level) + (j-1);
		if(grid == 2) // Vh
			return j*nx(grid,level) + (i-1);
		else if(grid == 3) // Ph
			return (i-1)*ny(grid,level) + (j-1);
		else // if(grid == 4) // Qh
			return i*ny(grid,level) + j;
	}
	bool inside() const
	{
		int i1 = i, j1 = j;
		if(grid == 1) j1--;
		else if(grid == 2) i1--;
		else if(grid == 3) i1--, j1--;
		return i1 >= 0 && i1 < nx(grid,level) && j1 >= 0 && j1 < ny(grid,level);
	}
	double x() const
	{
		if(grid == 1 || grid == 4)
			return std::min(std::max(0.0, (i*(1<<level)+0.5)*hx(0)), 1.0);
		else
			return std::min(std::max(0.0, ( (i-0.5)*(1<<level) + 0.5 )*hx(0)), 1.0);
	}
	double y() const
	{
		if(grid == 2 || grid == 4)
			return std::min(std::max(0.0, (j*(1<<level)+0.5)*hy(0)), 1.0);
		else
			return std::min(std::max(0.0, ( (j-0.5)*(1<<level) + 0.5 )*hy(0)), 1.0);
	}
	bool operator==(const pint &other)
	{
		return level == other.level && grid == other.grid && i == other.i && j == other.j;
	}
} pint;

// controls layout of DOFs in vector
struct layout
{
	int N = 1;
	int * grids = NULL;
	int * offsets = NULL;
	int * ncomps = NULL;
	int level = 0;
	layout(int N, const int _grids[], const int _ncomps[], int level) : N(N), level(level)
	{
		grids = new int[N];
		ncomps = new int[N];
		offsets = new int[N+1];
		memcpy(grids,_grids,sizeof(int)*N);
		memcpy(ncomps,_ncomps,sizeof(int)*N);
		offsets[0] = 0;
		for(int q = 1; q < N+1; ++q) offsets[q] = offsets[q-1] + ncomps[q-1] * pint::ndofs(grids[q-1], level);
	}
	layout(const layout& other) : N(other.N), level(other.level)
	{
		grids = new int[N];
		ncomps = new int[N];
		offsets = new int[N+1];
		memcpy(grids,other.grids,sizeof(int)*N);
		memcpy(ncomps,other.ncomps,sizeof(int)*N);
		memcpy(offsets,other.offsets,sizeof(int)*(N+1));
	}
	layout& operator=(const layout& other)
	{
		N = other.N;
		level = other.level;
		if(grids != NULL) delete[] grids;
		grids = new int[N];
		if(ncomps != NULL) delete[] ncomps;
		ncomps = new int[N];
		if(offsets != NULL) delete[] offsets;
		offsets = new int[N+1];
		memcpy(grids,other.grids,sizeof(int)*N);
		memcpy(ncomps,other.ncomps,sizeof(int)*N);
		memcpy(offsets,other.offsets,sizeof(int)*(N+1));
		return *this;
	}
	int size() const { return offsets[N]; }
	~layout()
	{
		delete[] grids;
		delete[] offsets;
		delete[] ncomps;
	}
	int offset(int pos) const
	{
		assert(pos < N+1);
		return offsets[pos];
	}
	int dof(const pint& p, int pos) const
	{
		if(p.inside()) return offset(pos) + ncomps[pos]*p.dof();
		throw "trying to access outside DOF";
		return -1;
	}
};

// multigrid layout
struct gmg_layout
{
	int nlevels = 1;
	layout ** lays = NULL;
	int * offsets = NULL; // offset for each level

	gmg_layout(int N, const int * grids, const int * ncomps, int nlevels)
		: nlevels(nlevels)
	{
		lays = new layout*[nlevels];
		offsets = new int[nlevels+1];
		offsets[0] = 0;
		for(int level = 0; level < nlevels; ++level)
		{
			lays[level] = new layout(N,grids,ncomps, level);
			offsets[level+1] = offsets[level] + lays[level]->size();
		}
	}
	gmg_layout(const gmg_layout& other) : nlevels(other.nlevels)
	{
		if(offsets != NULL) delete[] offsets;
		offsets = new int[nlevels+1];
		memcpy(offsets, other.offsets, (nlevels+1)*sizeof(int));
		if(lays != NULL) delete[] lays;
		lays = new layout*[nlevels];
		for(int level = 0; level < nlevels; ++level)
			lays[level] = new layout(*(other.lays[level]));
	}
	~gmg_layout()
	{
		for(int level = 0; level < nlevels; ++level)
		{
			delete lays[level];
			lays[level] = NULL;
		}
		delete[] lays;
		delete[] offsets;
	}
	int size() const {return offsets[nlevels];}
	int size(int level) const {assert(level < nlevels); return offsets[level+1]-offsets[level];}
	int dof(const pint& p, int pos) const
	{
		return lays[p.level]->dof(p,pos);
	}
};

template<typename T>
struct dof_vector
{
	std::vector<T> U;
	const gmg_layout* lay;

	dof_vector(const gmg_layout* lay) : lay(lay)
	{
		U.resize(lay->size());
	}

	int size() const {return lay->size();}
	int offset(int level) const {return lay->offsets[level];}
	int size(int level) const {return lay->lays[level]->size();}

	void assign(const dof_vector<T>& other, int level)
	{
		assert(level < lay->nlevels);
		std::copy(other.U.begin()+other.offset(level), other.U.begin()+other.offset(level+1), U.begin()+offset(level));
	}
	T& operator()(int igrid, const pint& p)
	{
		//return U[lay->dof(p, igrid)];
		return U[offset(p.level)+lay->dof(p, igrid)];
	}
	const T& operator()(int igrid, const pint& p) const
	{
		//return U[lay->dof(p, igrid)];
		return U[offset(p.level)+lay->dof(p, igrid)];
	}
};

double l2norm(const dof_vector<double>& UVP, int level = 0)
{
	double norm = 0.0;
	double hx = pint::hx(level), hy = pint::hy(level);
	int q = 0;
	for(int k = UVP.offset(level); k < UVP.offset(level+1); ++k)
	//	if(norm < fabs(UVP.U[k])) norm = fabs(UVP.U[k]), q = k;
	//return norm;
		norm += UVP.U[k]*UVP.U[k]*hx*hy;
	return sqrt(norm);
}

// (0 2)
// (2 1)
//typedef double symmat2[3];//xx yy xy
typedef std::array<double,3> symmat2;
// (0 2)
// (1 3)
//typedef double mat2[4];//xx yx xy yy
typedef std::array<double,4> mat2;
// (0 3 6)
// (1 4 7)
// (2 5 8)
//typedef double mat3[9];//xx yx zx  xy yy zy xz yz zz
typedef std::array<double,9> mat3;

inline double det2(const symmat2& A) {return A[0]*A[1]-A[2]*A[2];}
inline double det3(const mat3& A) {return A[0]*(A[4]*A[8]-A[5]*A[7]) - A[3]*(A[1]*A[8]-A[2]*A[7]) + A[6]*(A[1]*A[5]-A[2]*A[4]);}
inline double tr(const symmat2& A) {return A[0]+A[1];}

inline void eigs(const symmat2& A, double eigvals[2], mat2& eigvecs)
{
	double discr = tr(A)*tr(A) - 4*det2(A), d1, d2;
	if(discr < -1.0e-15) throw "Negative discriminant, should not happen!";
	discr = std::max(0.0, discr);
	eigvals[0] = 0.5*(tr(A) + sqrt(discr));
	eigvals[1] = 0.5*(tr(A) - sqrt(discr));
	bool err = fabs(eigvals[0]-eigvals[1]) < 1.0e-12;
	if(!err)
	{
		eigvecs[0] = A[2];
		eigvecs[1] = eigvals[0] - A[0];
		eigvecs[2] = eigvals[1] - A[1];
		eigvecs[3] = A[2];
		d1 = sqrt(eigvecs[0]*eigvecs[0]+eigvecs[1]*eigvecs[1]);
		if(d1)
		{
			eigvecs[0] /= d1;
			eigvecs[1] /= d1;
		}
		else err = true;
		d2 = sqrt(eigvecs[2]*eigvecs[2]+eigvecs[3]*eigvecs[3]);
		if(d2)
		{
			eigvecs[2] /= d2;
			eigvecs[3] /= d2;
		}
		else err = true;
	}
	if(err)
	{
		eigvecs[0] = 1.0; eigvecs[1] = 0.0;
		eigvecs[2] = 0.0; eigvecs[3] = 1.0;
	}
	if(fabs(eigvecs[0]*eigvecs[2] + eigvecs[1]*eigvecs[3]) > 1.0e-7)
		throw "non-orthogonal eigenvectors!";
}
inline void transpose(const mat2& A, mat2& At)
{
	At[0] = A[0];
	At[1] = A[2];
	At[2] = A[1];
	At[3] = A[3];
}
inline void matmul2ss(const symmat2& A, const symmat2& B, symmat2& AB)
{
	AB[0] = A[0]*B[0]+A[2]*B[2]; // xx
	AB[1] = A[2]*B[2]+A[1]*B[1]; // yy
	AB[2] = A[0]*B[2]+A[2]*B[1]; // xy
}
inline void matmul2ns(const mat2& A, const symmat2& B, mat2& AB)
{
	// ( A[0] A[2] ) ( B[0] B[2] )
	// ( A[1] A[3] ) ( B[2] B[1] )
	AB[0] = A[0]*B[0]+A[2]*B[2]; // xx
	AB[1] = A[1]*B[0]+A[3]*B[2]; // yx
	AB[2] = A[0]*B[2]+A[2]*B[1]; // xy
	AB[3] = A[1]*B[2]+A[3]*B[1]; // yy
}
inline void matmul2sn(const symmat2& A, const mat2& B, mat2& AB)
{
	// ( A[0] A[2] ) ( B[0] B[2] )
	// ( A[2] A[1] ) ( B[1] B[3] )
	AB[0] = A[0]*B[0]+A[2]*B[1]; // xx
	AB[1] = A[2]*B[0]+A[1]*B[1]; // yx
	AB[2] = A[0]*B[2]+A[2]*B[3]; // xy
	AB[3] = A[2]*B[2]+A[1]*B[3]; // yy
}
// B = R^T A R, B = B^T when A = A^T
inline void matmul3ns(const mat2& R, const symmat2& A, symmat2& RtAR)
{
	// ( R[0] R[2] ) -> ( r11 r12 ), ( A[0] A[2] )
	// ( R[1] R[3] ) -> ( r21 r22 ), ( A[2] A[1] )
	double r11 = R[0], r21 = R[1], r12 = R[2], r22 = R[3];
	double a11 = A[0], a21 = A[2], a12 = A[2], a22 = A[1];
	RtAR[0] = r11*r11*a11 + r11*r21*(a12+a21) + r21*r21*a22; // xx
	RtAR[1] = r12*r12*a11 + r12*r22*(a12+a21) + r22*r22*a22; // yy
	RtAR[2] = r11*r12*a11 + r11*r22*a12 + r12*r21*a21 + r21*r22*a22; // xy
}
// B = R^T A R
inline void matmul3nn(const mat2& R, const mat2& A, mat2& RtAR)
{
	// ( A[0] A[2] ) -> ( a11 a12 )
	// ( A[1] A[3] ) -> ( a21 a22 )
	double r11 = R[0], r21 = R[1], r12 = R[2], r22 = R[3];
	double a11 = A[0], a21 = A[1], a12 = A[2], a22 = A[3];
	RtAR[0] = r11*r11*a11 + r11*r21*(a12+a21) + r21*r21*a22; // xx
	RtAR[1] = r11*r12*a11 + r12*r21*a12 + r11*r22*a21 + r21*r22*a22; // yx
	RtAR[2] = r11*r12*a11 + r11*r22*a12 + r12*r21*a21 + r21*r22*a22; // xy
	RtAR[3] = r12*r12*a11 + r12*r22*(a12+a21) + r22*r22*a22; // yy
}
inline void solve3(mat3& A, const double b[3], double x[3])
{
	// Cramer method
	double v[3], d = det3(A);
	if(fabs(d) < 1.0e-12) throw "zero determinant";
	for(int i = 0; i < 3; ++i)
	{
		// overwrite column i of A with RHS vector
		for(int q = 0; q < 3; ++q)
		{
			v[q] = A[3*i+q];
			A[3*i+q] = b[q];
		}
		x[i] = det3(A)/d;
		// restore column i
		for(int q = 0; q < 3; ++q)
			A[3*i+q] = v[q];
	}
}
// Frobenius norm of symmetric 2x2 matrix
inline double frobnorm(const symmat2& psi)
{
	return sqrt(psi[0]*psi[0]+psi[1]*psi[1]+2*psi[2]*psi[2]);
}


void laplace(const pint& pdof, pint * dofs1, pint * dofs2, double mult)
{
	double coef1 = divgrad ? 2.0 : 1.0;
	int i = pdof.i, j = pdof.j, grid = pdof.grid, level = pdof.level;
	double hx = pint::hx(level), hy = pint::hy(level), hx0 = pint::hx(0), hy0 = pint::hy(0);
	if(grid == 1) // \Delta_h, x component
	{
		// Uh
		double coef[2], h[2] = {hy,hy};
		dofs1[0] = pint(1,i+1,j,level);
		dofs1[1] = pint(1,i-1,j,level);
		coef[0] = dofs1[0].inside() ? 1.0 : 0.5;
		coef[1] = dofs1[1].inside() ? 1.0 : 0.5;
		dofs1[0].coef = coef1/(hx*hx*coef[0]) * mult;
		dofs1[1].coef = coef1/(hx*hx*coef[1]) * mult;
		dofs1[2] = pint(1,i,j+1,level); if(!dofs1[2].inside()) h[0] = 0.5*(hy+hy0);
		dofs1[2].coef = 1.0/(hy*h[0]) * mult;
		dofs1[3] = pint(1,i,j-1,level); if(!dofs1[3].inside()) h[1] = 0.5*(hy+hy0);
		dofs1[3].coef = 1.0/(hy*h[1]) * mult;
		dofs1[4] = pint(1,i,j  ,level); dofs1[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs1[4].coef -= dofs1[k].coef;
		// Vh
		dofs2[0] = pint(2,i+1,j  ,level); //if(!dofs2[0].inside()) h[0] = 0.5*(hx+hx0);
		dofs2[0].coef =  1.0/(h[0]*hy) * mult;
		dofs2[1] = pint(2,i+1,j-1,level); //if(!dofs2[1].inside()) h[0] = 0.5*(hx+hx0);
		dofs2[1].coef = -1.0/(h[0]*hy) * mult;
		dofs2[2] = pint(2,i,j    ,level); //if(!dofs2[2].inside()) h[1] = 0.5*(hx+hx0);
		dofs2[2].coef = -1.0/(h[1]*hy) * mult;
		dofs2[3] = pint(2,i,j-1  ,level); //if(!dofs2[3].inside()) h[1] = 0.5*(hx+hx0);
		dofs2[3].coef =  1.0/(h[1]*hy) * mult;
	}
	else // \Delta_h, y component
	{
		// Vh
		double coef[2], h[2] = {hx,hx};
		dofs2[0] = pint(2,i,j+1,level);
		dofs2[1] = pint(2,i,j-1,level);
		coef[0] = dofs2[0].inside() ? 1.0 : 0.5;
		coef[1] = dofs2[1].inside() ? 1.0 : 0.5;
		dofs2[0].coef = coef1/(hy*hy*coef[0]) * mult;
		dofs2[1].coef = coef1/(hy*hy*coef[1]) * mult;
		dofs2[2] = pint(2,i+1,j,level); if(!dofs2[2].inside()) h[0] = 0.5*(hx+hx0);
		dofs2[2].coef = 1.0/(hx*h[0]) * mult;
		dofs2[3] = pint(2,i-1,j,level); if(!dofs2[3].inside()) h[1] = 0.5*(hx+hx0);
		dofs2[3].coef = 1.0/(hx*h[1]) * mult;
		dofs2[4] = pint(2,i,j  ,level); dofs2[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs2[4].coef -= dofs2[k].coef;
		// Uh
		dofs1[0] = pint(1,i,  j+1,level); //if(!dofs1[0].inside()) h[0] = 0.5*(hy+hy0);
		dofs1[0].coef =  1.0/(hx*h[0]) * mult;
		dofs1[1] = pint(1,i-1,j+1,level); //if(!dofs1[1].inside()) h[0] = 0.5*(hy+hy0);
		dofs1[1].coef = -1.0/(hx*h[0]) * mult;
		dofs1[2] = pint(1,i,j    ,level); //if(!dofs1[2].inside()) h[1] = 0.5*(hy+hy0);
		dofs1[2].coef = -1.0/(hx*h[1]) * mult;
		dofs1[3] = pint(1,i-1,j  ,level); //if(!dofs1[3].inside()) h[1] = 0.5*(hy+hy0);
		dofs1[3].coef =  1.0/(hx*h[1]) * mult;
	}
}

int div(const pint& pdof, pint * dofs1)
{
	assert(pdof.grid == 3);
	int i = pdof.i, j = pdof.j, level = pdof.level;
	double hx = pint::hx(level), hy = pint::hy(level);
	// Uh
	dofs1[0] = pint(1,i,j  ,level); dofs1[0].coef =  1.0 / hx;
	dofs1[1] = pint(1,i-1,j,level); dofs1[1].coef = -1.0 / hx;
	// Vh
	dofs1[2] = pint(2,i,j  ,level); dofs1[2].coef =  1.0 / hy;
	dofs1[3] = pint(2,i,j-1,level); dofs1[3].coef = -1.0 / hy;
	return 4;
}

// TODO: dp/dn = 0 on boundaries, needs review if otherwise
int grad(const pint& pdof, pint * dofs3)
{
	int i = pdof.i, j = pdof.j, grid = pdof.grid, level = pdof.level;
	double hx = pint::hx(level), hy = pint::hy(level);
	// Ph
	if(grid == 1) // x
	{
		dofs3[0] = pint(3,i+1,j,level); dofs3[0].coef =  1.0 / hx;
		dofs3[1] = pint(3,i,j  ,level); dofs3[1].coef = -1.0 / hx;
	}
	else //if(dir == 2) // y
	{
		assert(grid == 2);
		dofs3[0] = pint(3,i,j+1,level); dofs3[0].coef =  1.0 / hy;
		dofs3[1] = pint(3,i,j  ,level); dofs3[1].coef = -1.0 / hy;
	}
	return 2;
}

double minmod(double a, double b) {return 0.5*((2*(a>0)-1)+(2*(b>0)-1)) * std::min(fabs(a),fabs(b));}

// dir: +x -x +y -y
//       0  1  2  3
void reconstruct(const dof_vector<symmat2>& psi, const pint& pdof, symmat2& PSI_PLUS, symmat2& PSI_MINUS, int dir)
{
	if(dir < 0 || dir >= 4) throw "invalid direction, see reconstruct function";
	int i = pdof.i, j = pdof.j, level = pdof.level;
	int iLL = i, iL = i, iR = i, iRR = i;
	int jLL = j, jL = j, jR = j, jRR = j;
	if(dir < 2) // x
	{
		iLL = i-dir-1; iL = i-dir; iR = i-dir+1; iRR = i-dir+2;
	}
	else // y
	{
		dir -= 2;
		jLL = j-dir-1; jL = j-dir; jR = j-dir+1; jRR = j-dir+2;
	}
	pint pLL = pint(3,iLL,jLL,level), pL = pint(3,iL,jL,level), pR = pint(3,iR,jR,level), pRR = pint(3,iRR,jRR,level);

	if(!pL.inside() || !pR.inside())
	{
		pint p = pL.inside() ? pL : pR;
		const symmat2& PSI = psi(0,p);
		for(int q = 0; q < 3; ++q) PSI_PLUS[q] = PSI_MINUS[q] = PSI[q];
	}
	else
	{
		const symmat2& PSIL = psi(0,pL);
		const symmat2& PSIR = psi(0,pR);
		if(!pRR.inside())
			for(int q = 0; q < 3; ++q) PSI_PLUS[q] = PSIR[q];
		else
		{
			const symmat2& PSIRR = psi(0,pRR);
			for(int q = 0; q < 3; ++q)
				PSI_PLUS[q] = PSIR[q] - 0.5*minmod(PSIRR[q]-PSIR[q],PSIR[q]-PSIL[q]);
		}
		if(!pLL.inside())
			for(int q = 0; q < 3; ++q) PSI_MINUS[q] = PSIL[q];
		else
		{
			const symmat2& PSILL = psi(0,pLL);
			for(int q = 0; q < 3; ++q)
				PSI_MINUS[q] = PSIL[q] + 0.5*minmod(PSIR[q]-PSIL[q],PSIL[q]-PSILL[q]);
		}
	}
}

void flux_advection(const dof_vector<symmat2>& psi, const dof_vector<double>& x, const pint& pdof, symmat2& Nc)
{
	const double c = 0.5; // smoothing factor
	symmat2 PSI_PLUS, PSI_MINUS, H1, H0;
	double u1, u0, v1, v0;
	int i = pdof.i, j = pdof.j, level = pdof.level;
	pint pu1(1,i,j,level), pu0(1,i-1,j,level), pv1(2,i,j,level), pv0(2,i,j-1,level);
	double hx = pint::hx(level), hy = pint::hy(level);

	Nc[0] = Nc[1] = Nc[2] = 0.0;
	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 0); // +x
	u1 = x(0,pu1);//x(0,i,j,level);
	for(int q = 0; q < 3; ++q)
		H1[q] = u1*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(u1)*(PSI_PLUS[q] - PSI_MINUS[q]);
	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 1); // -x
	u0 = x(0,pu0);//x(0,i-1,j,level);
	for(int q = 0; q < 3; ++q)
		H0[q] = u0*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(u0)*(PSI_PLUS[q] - PSI_MINUS[q]);
	for(int q = 0; q < 3; ++q)
		Nc[q] += (H1[q]-H0[q])/hx;

	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 2); // +y
	v1 = x(1,pv1);//x(1,i,j,level);
	for(int q = 0; q < 3; ++q)
		H1[q] = v1*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(v1)*(PSI_PLUS[q] - PSI_MINUS[q]);
	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 3); // -y
	v0 = x(1,pv0);//x(1,i,j-1,level);
	for(int q = 0; q < 3; ++q)
		H0[q] = v0*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(v0)*(PSI_PLUS[q] - PSI_MINUS[q]);
	for(int q = 0; q < 3; ++q)
		Nc[q] += (H1[q]-H0[q])/hy;
}

//OLDROYD-B:
double g(const symmat2& PSI) {return 1.0;}
void P(const symmat2& PSI, symmat2& vP)
{
	// I - PSI (OLDROYD-B)
	vP[0] = 1.0 - PSI[0];
	vP[1] = 1.0 - PSI[1];
	vP[2] = 0.0 - PSI[2];
}

void decompose_grad(const mat2& GU, const symmat2& PSI, symmat2& B, mat2& OMEGA, symmat2& Nr)
{
	double lmbd[2];
	symmat2 diagM, eLMBD, ePSI, iePSI, PePSI, Nr1;
	mat2 M, R, Rt, OMG;
	//get eigenvalues and eigenvectors
	eigs(PSI, lmbd, R); // PSI = R * LMBD * R^T
	transpose(R, Rt); // should be R^T*R = I
	eLMBD[0] = exp(lmbd[0]); eLMBD[1] = exp(lmbd[1]); eLMBD[2] = 0.0;
	matmul3ns(Rt, eLMBD, ePSI);// exp(PSI) = R*exp(LMBD)*R^T
	eLMBD[0] = exp(-lmbd[0]); eLMBD[1] = exp(-lmbd[1]); eLMBD[2] = 0.0;
	matmul3ns(Rt, eLMBD, iePSI);// exp(-PSI) = R*exp(-LMBD)*R^T

	// Nr = (exp(-PSI) - I) / We
	double g1 = g(ePSI);
	P(ePSI, PePSI);
	//matmul2ss(iePSI, PePSI, Nr1);
	Nr[0] = g1/We * (iePSI[0]-1);//Nr1[0];
	Nr[1] = g1/We * (iePSI[1]-1);//Nr1[1];
	Nr[2] = g1/We * (iePSI[2]-0);//Nr1[2];

	if(fabs(lmbd[1]-lmbd[0]) < 1.0e-12)
	{
		// B = (GU+GU^T)/2, OMEGA = 0
		B[0] = GU[0];
		B[1] = GU[3];
		B[2] = 0.5*(GU[1]+GU[2]);
		OMEGA[0] = OMEGA[1] = OMEGA[2] = OMEGA[3] = 0.0; 
	}
	else
	{
		matmul3nn(R, GU, M); // M = R^T*GU*R
		// OMEGA = R*(0  w)*R^T = (0  w)
		//           (-w 0)       (-w 0)
		double w = (exp(lmbd[1])*M[2]+exp(lmbd[0])*M[1])/(exp(lmbd[1]) - exp(lmbd[0]));
		OMG[0] = OMG[3] = 0.0;
		OMG[1] = -w; OMG[2] = w;
		matmul3nn(Rt, OMG, OMEGA);
		// B = R*diag(M)*R^T
		diagM[0] = M[0]; diagM[1] = M[3]; diagM[2] = 0.0;
		matmul3ns(Rt, diagM, B);
	}
}

void fill_grad(const dof_vector<double>& U, const pint& pdof, mat2& GU)
{
	double u1, u0, v1, v0;
	int i = pdof.i, j = pdof.j, level = pdof.level;
	pint p11(1,i,j,level), p10(1,i-1,j,level), p01(2,i,j,level), p00(2,i,j-1,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	u1 = U(0,p11);
	u0 = U(0,p10);
	v1 = U(1,p01);
	v0 = U(1,p00);
	GU[0] = (u1-u0)/hx; //u_x
	GU[3] = (v1-v0)/hy; //v_y
	//u_y 
	p11 = pint(1,i  ,j+1,level);
	p10 = pint(1,i-1,j+1,level);
	p01 = pint(1,i  ,j-1,level);
	p00 = pint(1,i-1,j-1,level);
	if(p10.inside() && p11.inside())
		u1 = 0.5*( U(0,p11) + U(0,p10) );
	else
		u1 = exact(test, 1, 0.5*(p11.x()+p10.x()), 0.5*(p11.y()+p10.y()) );
	if(p00.inside() && p01.inside())
		u0 = 0.5*( U(0,p01) + U(0,p00) );
	else
		u0 = exact(test, 1, 0.5*(p01.x()+p00.x()), 0.5*(p01.y()+p00.y()) );
	GU[2] = (u1-u0)/(2*hy); // u_y
	//v_x
	p11 = pint(2,i+1,j  ,level);
	p10 = pint(2,i+1,j-1,level);
	p01 = pint(2,i-1,j  ,level);
	p00 = pint(2,i-1,j-1,level);
	if(p10.inside() && p11.inside())
		v1 = 0.5*( U(1,p11) + U(1,p10) );
	else
		v1 = exact(test, 2, 0.5*(p11.x()+p10.x()), 0.5*(p11.y()+p10.y()) );
	if(p00.inside() && p01.inside())
		v0 = 0.5*( U(1,p01) + U(1,p00) );
	else
		v0 = exact(test, 2, 0.5*(p01.x()+p00.x()), 0.5*(p01.y()+p00.y()) );
	GU[1] = (v1-v0)/(2*hx); // v_x
}

void update_psi(const dof_vector<double>& x, const dof_vector<double>& x0, dof_vector<symmat2>& psi, dof_vector<symmat2>& psi0, dof_vector<symmat2>& psi1)
{
	double time_rel = fabs(dt0) < 1.0e-12 ? 0.0 : dt/dt0;
	double a0 = (2*time_rel+1)/(time_rel+1), a1 = -(1+time_rel), a2 = time_rel*time_rel/(1+time_rel);
	double b1 = 1+time_rel, b2 = -time_rel;
	int level = 0; // TODO
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		mat2 OMEGA, OMEGA0, GU, GU0;
		symmat2 B, B0, Nr, Nr0, Nc, Nc0;
		const symmat2& PSI = psi(0,pdof);
		const symmat2& PSI0 = psi0(0,pdof);
		mat3 A = {0,0,0,0,0,0,0,0,0};
		double rhs[3] = {0,0,0}, sol[3];
		// time derivative: a0*phi(n+1) + a1*psi(n) + a2*psi(n-1)
		A[0+0*3] = A[1+1*3] = A[2+2*3] = a0;
		for(int q = 0; q < 3; ++q)
			rhs[q] -= (a1*PSI[q] + a2*PSI0[q]);
		fill_grad(x, pdof, GU);
		fill_grad(x0, pdof, GU0);
		decompose_grad(GU, PSI, B, OMEGA, Nr);
		decompose_grad(GU0, PSI0, B0, OMEGA0, Nr0);
		// OMEGA*PSI(n+1) - PSI(n+1)*OMEGA = (0  w) (p[0] p[2]) - (p[0] p[2]) (0  w)
		//                                   (-w 0) (p[2] p[1]) - (p[2] p[1]) (-w 0)
		double w = OMEGA[2], w0 = OMEGA0[2];
		A[0+2*3] -= dt*(b1 * 2*w + b2 * 2*w0); // xx: -2*w*PSI_xy
		A[1+2*3] += dt*(b1 * 2*w + b2 * 2*w0); // yy:  2*w*PSI_xy
		A[2+0*3] += dt*(b1 * w   + b2 * w0);   // xy:    w*PSI_xx
		A[2+1*3] -= dt*(b1 * w   + b2 * w0);   // xy:  - w*PSI_yy
		
		// Nc(UV, PSI) = (UV * nabla) PSI
		flux_advection(psi, x, pdof, Nc);
		flux_advection(psi0, x0, pdof, Nc0);
		// Nr(PSI) + 2*B - Nc(UV, PSI) into rhs
		for(int q = 0; q < 3; ++q)
		{
			rhs[q] += dt*(b1*2*B[q] + b2*2*B0[q]);
			rhs[q] += dt*(b1*Nr[q] + b2*Nr0[q]);
			rhs[q] -= dt*(b1*Nc[q] + b2*Nc0[q]);
		}
		// solve 3x3 system
		solve3(A, rhs, sol);
		for(int q = 0; q < 3; ++q)
			psi1(0,pdof)[q] = sol[q];
	}
	
	psi0.assign(psi , level);// psi(n-1) <- psi(n)
	psi.assign (psi1, level);// psi(n)   <- psi(n+1)
}

void update_tau(const dof_vector<symmat2>& psi, dof_vector<symmat2>& tau1, dof_vector<symmat2>& tau)
{
	int level = 0; // TODO
	// remember tau before it is overwritten
	tau1.assign(tau,level);
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	// update tau: tau = g(exp(psi)) / We * (exp(psi) -I)
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		symmat2 ePSI, eLMBD;
		mat2 R, Rt;
		double lmbd[2];
		// compute exp(psi)
		eigs(psi(0,pdof), lmbd, R);
		transpose(R, Rt);
		eLMBD[0] = exp(lmbd[0]); eLMBD[1] = exp(lmbd[1]); eLMBD[2] = 0.0;
		matmul3ns(Rt, eLMBD, ePSI);
		double g1 = g(ePSI);
		tau(0,pdof)[0] = g1/We * (ePSI[0] - 1);
		tau(0,pdof)[1] = g1/We * (ePSI[1] - 1);
		tau(0,pdof)[2] = g1/We * (ePSI[2] - 0);
	}
}

double residual(const dof_vector<symmat2>& tau0, const dof_vector<symmat2>& tau)
{
	int level = 0; // TODO
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	double l2norm = 0.0;
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(3,i,j,level);
		double err = 0.0;
		for(int q = 0; q < 3; ++q)
			err = std::max(fabs(tau0(0,pdof)[q]-tau(0,pdof)[q]), err);
		l2norm += err*err*hx*hy;
	}
	return sqrt(l2norm);
}

int div_tau(const pint& pdof, std::array<pint,6>& dofs)
{
	int i = pdof.i, j = pdof.j, grid = pdof.grid, level = pdof.level;
	double hx = pint::hx(level), hy = pint::hy(level);
	if(grid == 1) // x
	{
		// xx
		dofs[0] = pint(3,i+1,j,level);
		dofs[1] = pint(3,i  ,j,level);
		dofs[0].coef =  1.0 / hx;
		dofs[1].coef = -1.0 / hx;
		// xy
		dofs[2] = pint(3,i  ,j+1,level);
		dofs[3] = pint(3,i+1,j+1,level);
		dofs[2].coef =  0.5 / hy;
		dofs[3].coef =  0.5 / hy;
		if(!dofs[2].inside() || !dofs[3].inside())
		{
			dofs[2] = pint(3,i  ,j,level);dofs[2].coef *= 2;
			dofs[3] = pint(3,i+1,j,level);dofs[3].coef *= 2;
		}
		dofs[4] = pint(3,i  ,j-1,level);
		dofs[5] = pint(3,i+1,j-1,level);
		dofs[4].coef = -0.5 / hy;
		dofs[5].coef = -0.5 / hy;
		if(!dofs[4].inside() || !dofs[5].inside())
		{
			dofs[4] = pint(3,i  ,j,level);dofs[4].coef *= 2;
			dofs[5] = pint(3,i+1,j,level);dofs[5].coef *= 2;
		}
	}
	else //if(grid == 2) // y
	{
		assert(grid == 2);
		// yy
		dofs[0] = pint(3,i,j+1,level);
		dofs[1] = pint(3,i,j  ,level);
		dofs[0].coef =  1.0 / hy;
		dofs[1].coef = -1.0 / hy;
		// xy
		dofs[2] = pint(3,i+1,j  ,level);
		dofs[3] = pint(3,i+1,j+1,level);
		dofs[2].coef =  0.5 / hx;
		dofs[3].coef =  0.5 / hx;
		if(!dofs[2].inside() || !dofs[3].inside())
		{
			dofs[2] = pint(3,i,j  ,level);dofs[2].coef *= 2;
			dofs[3] = pint(3,i,j+1,level);dofs[3].coef *= 2;
		}
		dofs[4] = pint(3,i-1,j  ,level);
		dofs[5] = pint(3,i-1,j+1,level);
		dofs[4].coef = -0.5 / hx;
		dofs[5].coef = -0.5 / hx;
		if(!dofs[4].inside() || !dofs[5].inside())
		{
			dofs[4] = pint(3,i,j  ,level);dofs[4].coef *= 2;
			dofs[5] = pint(3,i,j+1,level);dofs[5].coef *= 2;
		}
	}
	return 6;
}

void fill_rhs(const dof_vector<symmat2>& tau, dof_vector<double>& b)
{
	std::array<pint,6> dofs;
	int ndofs = 6;
	int level = 0; // TODO
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	double scale = pint::hx(level)*pint::hy(level);
	// grid1
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(1,i,j,level);
		div_tau(pdof,dofs);
		for(int q = 0; q < ndofs/2; ++q)
			//if(dofs[q].inside())
			if( dofs[2*q].inside() && dofs[2*q+1].inside() )
			for(int k = 0; k < 2; ++k)
		{
			int q1 = 2*q+k;
			int d = q1 < 2 ? 0 : 2;//xx(0) xy(2)
			b(0,pdof) += dofs[q1].coef * mu_p*tau(0,dofs[q1])[d] * scale;
		}
	}
	// grid2
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		pint pdof(2,i,j,level);
		div_tau(pdof,dofs);
		for(int q = 0; q < ndofs/2; ++q)
			//if(dofs[q].inside())
			if( dofs[2*q].inside() && dofs[2*q+1].inside() )
			for(int k = 0; k < 2; ++k)
		{
			int q1 = 2*q+k;
			int d = q1 < 2 ? 1 : 2;//yy(1) xy(2)
			b(1,pdof) += dofs[q1].coef * mu_p*tau(0,dofs[q1])[d] * scale;
		}
	}
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

	double alpha_uv = 0.8;
	for(int q = 0; q < 4; ++q) D[q] /= alpha_uv;

	bool fail = solve_bordered<5>(D,U,L,RHS,SOL);
	UVP(0,puvp[0]) = SOL[0];
	UVP(0,puvp[1]) = SOL[1];
	UVP(1,puvp[2]) = SOL[2];
	UVP(1,puvp[3]) = SOL[3];
	if(!fail) UVP(2,puvp[4]) = SOL[4];
}

void vanka_smoother(int iters, dof_vector<double>& UVP, const dof_vector<double>& RHS0, int level)
{
	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
	static bool fwd = true;
	for(int iter = 0; iter < iters; ++iter)
	{
		if(fwd)
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

void save_vtk(std::string filename, const dof_vector<symmat2>& psi, const dof_vector<symmat2>& tau, const dof_vector<double>& x, int level = 0)
{
	std::ofstream ofs(filename);

	int n1 = pint::nx(1,level), n2 = pint::ny(2,level);
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
		ofs << i*pint::hx(level) << ' ';
	ofs << '\n';
	ofs << "Y_COORDINATES " << (n1+1) << " double\n";
	for(int j = 0; j < n2+1; ++j)
		ofs << j*pint::hy(level) << ' ';
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
		pint pdof(3,i,j,level);
		double p = 0.0;
		if(i > 0 && i < n1 && j > 0 && j < n2)
			p = x(2,pdof);
		ofs << p << '\n';
	}
	//divergence
	ofs << "SCALARS div_U double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 0; j < n2+1; ++j)
		for(int i = 0; i < n1+1; ++i)
	{
		double divu = 0.0;
		pint dofs[8], dof = pint(3,i,j,level);
		if(i > 0 && i < n1 && j > 0 && j < n2)
		{
			int ndofs = div(dof,dofs);
			for(int q = 0; q < ndofs; ++q)
				if( dofs[q].inside() )
					divu += x(dofs[q].grid-1,dofs[q]) * dofs[q].coef;
				else
					divu += exact(test, q < 2 ? 1 : 2, dofs[q].x(), dofs[q].y()) * dofs[q].coef;
		}
		ofs << divu << '\n';
	}
	// tau
	const char* tau_name[3] = {"tau_xx","tau_yy","tau_xy"};
	for(int k = 0; k < 3; ++k) if(level == 0)
	{
		ofs << "SCALARS " << tau_name[k] << " double 1\n";
		ofs << "LOOKUP_TABLE default\n";
		for(int j = 0; j < n2+1; ++j)
			for(int i = 0; i < n1+1; ++i)
		{
			pint pdof(3,i,j,level);
			double tau_ij = 0.0;
			if(i > 0 && i < n1 && j > 0 && j < n2)
				tau_ij = tau(0,pdof)[k];
			if(std::isnan(tau_ij) || std::isinf(tau_ij))
				tau_ij = -1e10;
			ofs << tau_ij << '\n';
		}
	}
	// psi
	const char* psi_name[3] = {"psi_xx","psi_yy","psi_xy"};
	for(int k = 0; k < 3; ++k) if(level == 0)
	{
		ofs << "SCALARS " << psi_name[k] << " double 1\n";
		ofs << "LOOKUP_TABLE default\n";
		for(int j = 0; j < n2+1; ++j)
			for(int i = 0; i < n1+1; ++i)
		{
			pint pdof(3,i,j,level);
			double psi_ij = 0.0;
			if(i > 0 && i < n1 && j > 0 && j < n2)
				psi_ij = psi(0,pdof)[k];
			if(std::isnan(psi_ij) || std::isinf(psi_ij))
				psi_ij = -1e10;
			ofs << psi_ij << '\n';
		}
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
			dofs[0] = pint(1,i,j  ,level);
			dofs[1] = pint(1,i-1,j,level);
			for(int k = 0; k < 2; ++k)
				if(dofs[k].inside()) u[k] = x(0, dofs[k]);
				else u[0] = u[1] = exact(test, 1, dofs[k].x(), dofs[k].y());
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
			dofs[0] = pint(2,i,j  ,level);
			dofs[1] = pint(2,i,j-1,level);
			for(int k = 0; k < 2; ++k)
				if(dofs[k].inside()) u[k] = x(1, dofs[k]);
				else u[0] = u[1] = exact(test, 2, dofs[k].x(), dofs[k].y());
		}
		ofs << 0.5*(u[0]+u[1]) << '\n';
	}

	ofs.close();
}
void save_vtk(std::string filename, const dof_vector<double>& x, const gmg_layout& lay2, int level = 0, const dof_vector<symmat2>* psi = NULL, const dof_vector<symmat2>* tau = NULL)
{
	if(psi == NULL || tau == NULL)
	{
		dof_vector<symmat2> psi(&lay2);
		save_vtk(filename, psi, psi, x, level);
	}
	else save_vtk(filename, *psi, *tau, x, level);
}

#if 1
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
#endif


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

int main2()
{
	symmat2 A1;
	A1[0] = 1e-2; A1[1] = -1e-2; A1[2] = 1e-5;
	mat2 eigvecs;
	double eigvals[2];
	eigs(A1, eigvals, eigvecs);
	std::cout << "eigvals: b1 " << eigvals[0] << " b2 " << eigvals[1] << std::endl;
	std::cout << "eigvecs: v1 " << eigvecs[0] << " " << eigvecs[1] << std::endl;
	std::cout << "eigvecs: v2 " << eigvecs[2] << " " << eigvecs[3] << std::endl;
	mat2 AV;
	matmul2sn(A1, eigvecs, AV);
	std::cout << "matvec: A*v1 " << AV[0] << " " << AV[1]
		<< " A*v2 " << AV[2] << " " << AV[3] << std::endl;
	std::cout << "matvec: A*v1-b1*v1 " << (AV[0]-eigvals[0]*eigvecs[0]) << " " << (AV[1]-eigvals[0]*eigvecs[1]) << std::endl;
	std::cout << "matvec: A*v2-b2*v2 " << (AV[2]-eigvals[1]*eigvecs[2]) << " " << (AV[3]-eigvals[1]*eigvecs[3]) << std::endl;

	mat2 A, AB, C, G, AGAt, H, At;
	symmat2 B, D, E, DB, ABAt, F;
	A[0] = 1; A[1] = 3; A[2] = 2; A[3] = 4;
	B[0] = 5; B[1] = 8; B[2] = 6;
	G[0] = 5; G[1] = 6; G[2] = 7; G[3] = 8;
	C[0] = 17; C[1] = 39; C[2] = 22; C[3] = 50;
	F[0] = 61; F[1] = 317; F[2] = 139;
	D[0] = 1; D[1] = 1; D[2] = 2;
	E[0] = 17; E[1] = 20; E[2] = 22;
	H[0] = 63; H[1] = 145; H[2] = 143; H[3] = 329;
	matmul2ns(A,B,AB);
	std::cout << "A*B - AB = " << (AB[0]-C[0]) << " " << (AB[1]-C[1]) << " " << (AB[2]-C[2]) << " " << (AB[3]-C[3]) << std::endl;
	matmul2ss(D,B,DB);
	std::cout << "D*B - DB = " << (DB[0]-E[0]) << " " << (DB[1]-E[1]) << " " << (DB[2]-E[2]) << std::endl;

	transpose(A,At);
	matmul3nn(At,G,AGAt);
	std::cout << "A*G*At - AGAt = " << (AGAt[0]-H[0]) << " " << (AGAt[1]-H[1]) << " " << (AGAt[2]-H[2]) << " " << (AGAt[3]-H[3]) << std::endl;

	matmul3ns(At,B,ABAt);
	std::cout << "A*B*At - ABAt = " << (ABAt[0]-F[0]) << " " << (ABAt[1]-F[1]) << " " << (ABAt[2]-F[2]) << std::endl;

	std::cout << "det2 is ok? " << (fabs(det2(B)-4) < 1.0e-12) << std::endl;

	mat3 A3 = {1,0,0,2,5,6,3,6,8};
	std::cout << "det3 is ok? " << (fabs(det3(A3)-4) < 1.0e-12) << std::endl;

	double rhs[3] = {14,28,36}, x[3] = {1,2,3}, sol[3];
	solve3(A3, rhs, sol);
	std::cout << "solve3: sol " << sol[0] << " " << sol[1] << " " << sol[2] << "; expected " << x[0] << " " << x[1] << " " << x[2] << std::endl;

	double D1[4] = {1,2,3,4};
	double U1[4] = {5,6,7,8};
	double L1[4] = {9,10,11,12};
	double RHS1[5] = {-26,-26,-44,-24,26};
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
		//if(!resid) std::fill(UVP.U.begin()+UVP.offset(level),UVP.U.begin()+UVP.offset(level+1),0.0);
	}

	pint dofs1[5], dofs2[5], dofs3[2];
	int ndofs;
	// equations for Uh
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		pint pdof(1,i,j,level);
		int k = lay1.dof(pdof, 0);
		laplace(pdof,dofs1,dofs2, -mu_s);
		if(level == 0) RHS0(0,pdof) = rhs_exact(test, 1,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() )
			{
				if(mgcl && !filled && !resid)
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
					if(mgcl && !filled && !resid)
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
				if(mgcl && !filled && !resid)
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
		if(mgcl && !filled) A->FinalizeRow();
#endif
	}
	// equations for Vh
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		pint pdof(2,i,j,level);
		int k = lay1.dof(pdof, 1);
		laplace(pdof,dofs1,dofs2, -mu_s);
		if(level == 0) RHS0(1,pdof) = rhs_exact(test, 2,pdof.x(),pdof.y()) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs2[q].inside() )
			{
				if(mgcl && !filled && !resid)
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
					if(mgcl && !filled && !resid)
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
				if(mgcl && !filled && !resid)
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
		if(mgcl && !filled) A->FinalizeRow();
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
				if(mgcl && !filled && !resid)
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
		if(mgcl && !filled) A->FinalizeRow();
#endif
	}
	if(mgcl && !filled)
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
	return l2norm(res, 0);
}

void multigrid_coarse_level::solve()
{
	dof_vector<double>& UVP = *(common.UVP), &RHS0 = *(common.RHS0);
	//std::fill(RHS0.U.begin()+RHS0.offset(level), RHS0.U.begin()+RHS0.offset(level+1), 0.0);
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
	fill_stokes(this, level == 0); // RHS0 <- RHS0 - A*UVP
	save_vtk("tmp2.vtk", RHS0, *(common.lay2), level);
	restriction(RHS0);
	if(1)
	{
		res.assign(RHS0, level+1);
		fill_stokes(res, next);
		double resid = l2norm(res, level+1);
		std::cout << "level " << (level+1) << " restricted residual " << resid << "\n" << std::endl;
	}
	save_vtk("tmp3.vtk", RHS0, *(common.lay2), level+1);
	if(level < common.nlevels-1) std::fill(UVP.U.begin()+UVP.offset(level+1), UVP.U.begin()+UVP.offset(level+2), 0.0);
	for(int q = 0; q < common.schedule; ++q)
	{
		next->solve();
		save_vtk("tmp4.vtk", UVP, *(common.lay2), level+1);
		next->prolongation(*(common.UVP));
	}
	if(1)
	{
		res.assign(RHS0, level);
		fill_stokes(res, this, true);
		double resid = l2norm(res, level);
		std::cout << "level " << level << " prolongated residual " << resid << "\n" << std::endl;
	}
	save_vtk("tmp5.vtk", UVP, *(common.lay2), level);
	//smoothing(common.postsmooth_iters);
}


double func(int grid, double x, double y)
{
	if(grid == 1)
		return x+y+2;
	else if(grid == 2)
		return x-y-1;
	else if(grid == 3)
		return x;
	else return 0.0;
}

int main3()
{
	n1 = pow(2, glev)+1;//atoi(argv[1]);
	n2 = n1;//argv[2];
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

int main(int argc, char ** argv)
{
#if !defined(USE_S3M)
	Solver::Initialize(&argc, &argv, "database.xml");
#endif
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
	n1 = pow(2, glev)+1;//atoi(argv[1]);
	n2 = n1;//argv[2];
	hx = 1.0/n1;
	hy = 1.0/n2;
	if(argc > 2)
		T = atof(argv[2]);
	if(argc > 3)
		nframe = atoi(argv[3]);
	if(argc > 4)
		We = atof(argv[4]);

	test = 1;
	
	int nlevels = 1;
	int level = 0; // TODO
	double hx = pint::hx(level), hy = pint::hy(level);
	dt = cfl * std::min(hx,hy);
	std::cout << "n1 " << pint::nx(1,level) << " n2 " << pint::ny(2,level) << " test " << test << " dt " << dt << std::endl;
	std::cout << "nframe " << nframe << " T " << T << " We " << We << std::endl;
	
	const int grids1[3] = {1,2,3};	// Uh(1), Vh(2), Ph(3)
	const int ncomps1[3] = {1,1,1};
	const int grids2[1] = {3};	// Psi(3), Psi=(Psi_xx, Psi_yy, Psi_xy)
	const int ncomps2[1] = {3};

	gmg_layout lay1(3,grids1,ncomps1,nlevels), lay2(1,grids2,ncomps2,1), lay3(3,grids1,ncomps1,1);
	dof_vector<double> UVP(&lay1), UVP0(&lay3), RHS(&lay1), RHS0(&lay3);
	dof_vector<symmat2> PSI(&lay2), PSI0(&lay2), PSI1(&lay2), TAU(&lay2), TAU0(&lay2), TAU1(&lay2);

	int nstokes = lay1.size(level);//n1*(n2-1)+n2*(n1-1)+(n1-1)*(n2-1); // Uh,Vh,Ph
	int ntau = lay2.size(level);//3*(n1-1)*(n2-1);
	std::cout << "nstokes " << nstokes << " ntau " << ntau << std::endl;
	// psi = (psi_xx (Ph), psi_yy (Ph), psi_xy (Ph))
	// psi -- psi(n), psi0 -- psi(n-1), psi1 -- psi(n+1)
	// tau -- tau(n+1), tau0 -- tau(n-1), tau1 -- tau(n)

	//if(integral_constraint) nstokes += 1; // integral constraint on pressure

	if(test) setup_exact(UVP);
	save_vtk("tmp0.vtk", PSI, TAU, UVP);
	UVP.assign(UVP0, 0); //save_vtk("tmp1.vtk", PSI, TAU, UVP);

	bool filled = false;
	bool success = false, finish = false;
	double res = 0.0, atol = 1.0e-4;
	int iter = 0, maxiters = 10000;

	multigrid_params common; common.nlevels = nlevels; common.lay2 = &lay2;
	common.UVP = &UVP; common.RHS0 = &RHS;
	multigrid gmg(common);
	do
	{
		std::fill(RHS.U.begin()+RHS.offset(0),RHS.U.begin()+RHS.offset(1),0.0);
		fill_stokes(RHS, gmg.head, false);
		//RHS.assign(RHS0, 0); // RHS <- RHS0
		if(test == 0) fill_rhs(TAU, RHS); // RHS <- RHS + div(TAU)
		success = gmg.solve();
		if(test)
		{
			finish = true;
			save_vtk("tmp6.vtk", PSI, TAU, UVP);
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

#if !defined(USE_S3M)
	Solver::Finalize();
#endif
	return 0;
}
