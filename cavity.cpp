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
	int ncomps = 1;
	pint(int grid = -1, int i = -10, int j = -10, int ncomps = 1) : i(i), j(j), grid(grid), ncomps(ncomps) {}
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
			return ncomps*(i*ny(grid) + (j-1));
		if(grid == 2) // Vh
			return ncomps*(j*nx(grid) + (i-1));
		else if(grid == 3) // Ph
			return ncomps*((i-1)*ny(grid) + (j-1));
		else // if(grid == 4) // Qh
			return ncomps*(i*ny(grid) + j);
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
	int * ncomps;
	layout(int N, const int _grids[], const int _ncomps[]) : N(N)
	{
		grids = new int[N];
		ncomps = new int[N];
		offsets = new int[N+1];
		memcpy(grids,_grids,sizeof(int)*N);
		memcpy(ncomps,_ncomps,sizeof(int)*N);
		offsets[0] = 0;
		for(int q = 1; q < N+1; ++q) offsets[q] = offsets[q-1] + ncomps[q-1] * pint::ndofs(grids[q-1]);
	}
	layout(const layout& other) : N(other.N)
	{
		grids = new int[N];
		ncomps = new int[N];
		offsets = new int[N+1];
		memcpy(grids,other.grids,sizeof(int)*N);
		memcpy(ncomps,other.ncomps,sizeof(int)*N);
		memcpy(offsets,other.offsets,sizeof(int)*(N+1));
	}
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
		if(p.ncomps != ncomps[pos]) throw "inconsistent amount of DOF components: pint::ncomps vs. layout.ncomps";
		if(p.inside()) return offset(pos) + p.dof();
		throw "trying to access outside DOF";
		return -1;
	}
};

template<typename T>
struct dof_vector
{
	std::vector<T> U;
	layout* lay;

	dof_vector(layout* lay) : lay(lay)
	{
		int N = lay->offset(lay->N);
		U.resize(N);
	}
	//dof_vector(const dof_vector& other) : lay(other.lay), U(other.U) {}

	int ndofs() const {return lay->offset(lay->N);}
	T& operator()(int igrid, int i, int j)
	{
		return U[lay->dof(pint(lay->grids[igrid],i,j,lay->ncomps[igrid]), igrid)];
	}
	const T& operator()(int igrid, int i, int j) const
	{
		return U[lay->dof(pint(lay->grids[igrid],i,j,lay->ncomps[igrid]), igrid)];
	}
	T& operator[](unsigned pos) {return U[pos];}
	const T& operator[](unsigned pos) const {return U[pos];}
};

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

double minmod(double a, double b) {return 0.5*((2*(a>0)-1)+(2*(b>0)-1)) * std::min(fabs(a),fabs(b));}

// dir: +x -x +y -y
//       0  1  2  3
void reconstruct(dof_vector<symmat2>& psi, int i, int j, symmat2& PSI_PLUS, symmat2& PSI_MINUS, int dir)
{
	if(dir < 0 || dir >= 4) throw "invalid direction, see reconstruct function";
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
	pint pLL = pint(3,iLL,jLL,3), pL = pint(3,iL,jL,3), pR = pint(3,iR,jR,3), pRR = pint(3,iRR,jRR,3);

	if(!pL.inside() || !pR.inside())
	{
		pint p = pL.inside() ? pL : pR;
		const symmat2 &PSI = psi(0,p.i,p.j);
		for(int q = 0; q < 3; ++q) PSI_PLUS[q] = PSI_MINUS[q] = PSI[q];
	}
	else
	{
		const symmat2& PSIL = psi(0,iL,jL);
		const symmat2& PSIR = psi(0,iR,jR);
		if(!pRR.inside())
			for(int q = 0; q < 3; ++q) PSI_PLUS[q] = PSIR[q];
		else
		{
			const symmat2& PSIRR = psi(0,iRR,jRR);
			for(int q = 0; q < 3; ++q)
				PSI_PLUS[q] = PSIR[q] - 0.5*minmod(PSIRR[q]-PSIR[q],PSIR[q]-PSIL[q]);
		}
		if(!pLL.inside())
			for(int q = 0; q < 3; ++q) PSI_MINUS[q] = PSIL[q];
		else
		{
			const symmat2& PSILL = psi(0,iLL,jLL);
			for(int q = 0; q < 3; ++q)
				PSI_MINUS[q] = PSIL[q] + 0.5*minmod(PSIR[q]-PSIL[q],PSIL[q]-PSILL[q]);
		}
	}
}

void flux_advection(const layout& lay1, dof_vector<symmat2>& psi, const dof_vector<double>& x, int i, int j, symmat2& Nc)
{
	const double c = 0.5; // smoothing factor
	symmat2 PSI_PLUS, PSI_MINUS, H1, H0;
	double u1, u0, v1, v0;

	Nc[0] = Nc[1] = Nc[2] = 0.0;
	reconstruct(psi, i, j, PSI_PLUS, PSI_MINUS, 0); // +x
	u1 = x[lay1.dof(pint(1,i,j), 0)];
	for(int q = 0; q < 3; ++q)
		H1[q] = u1*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(u1)*(PSI_PLUS[q] - PSI_MINUS[q]);
	reconstruct(psi, i, j, PSI_PLUS, PSI_MINUS, 1); // -x
	u0 = x[lay1.dof(pint(1,i-1,j), 0)];
	for(int q = 0; q < 3; ++q)
		H0[q] = u0*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(u0)*(PSI_PLUS[q] - PSI_MINUS[q]);
	for(int q = 0; q < 3; ++q)
		Nc[q] += (H1[q]-H0[q])/hx;

	reconstruct(psi, i, j, PSI_PLUS, PSI_MINUS, 2); // +y
	v1 = x[lay1.dof(pint(2,i,j), 1)];
	for(int q = 0; q < 3; ++q)
		H1[q] = v1*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(v1)*(PSI_PLUS[q] - PSI_MINUS[q]);
	reconstruct(psi, i, j, PSI_PLUS, PSI_MINUS, 3); // -y
	v0 = x[lay1.dof(pint(2,i,j-1), 1)];
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

void fill_grad(const layout& lay1, const dof_vector<double>& x, int i, int j, mat2& GU)
{
	double u1, u0, v1, v0;
	pint p11, p10, p01, p00;
	u1 = x[lay1.dof(pint(1,i  ,j), 0)];
	u0 = x[lay1.dof(pint(1,i-1,j), 0)];
	v1 = x[lay1.dof(pint(2,i  ,j), 1)];
	v0 = x[lay1.dof(pint(2,i,j-1), 1)];
	GU[0] = (u1-u0)/hx; //u_x
	GU[3] = (v1-v0)/hy; //v_y
	//u_y 
	p11 = pint(1,i  ,j+1);
	p10 = pint(1,i-1,j+1);
	p01 = pint(1,i  ,j-1);
	p00 = pint(1,i-1,j-1);
	if(p10.inside() && p11.inside())
		u1 = 0.5*( x[lay1.dof(p11,0)] + x[lay1.dof(p10,0)] );
	else
		u1 = exact(test, 1, 0.5*(p11.x()+p10.x()), 0.5*(p11.y()+p10.y()) );
	if(p00.inside() && p01.inside())
		u0 = 0.5*( x[lay1.dof(p01,0)] + x[lay1.dof(p00,0)] );
	else
		u0 = exact(test, 1, 0.5*(p01.x()+p00.x()), 0.5*(p01.y()+p00.y()) );
	GU[2] = (u1-u0)/(2*hy); // u_y
	//v_x
	p11 = pint(2,i+1,j  );
	p10 = pint(2,i+1,j-1);
	p01 = pint(2,i-1,j  );
	p00 = pint(2,i-1,j-1);
	if(p10.inside() && p11.inside())
		v1 = 0.5*( x[lay1.dof(p11,1)] + x[lay1.dof(p10,1)] );
	else
		v1 = exact(test, 2, 0.5*(p11.x()+p10.x()), 0.5*(p11.y()+p10.y()) );
	if(p00.inside() && p01.inside())
		v0 = 0.5*( x[lay1.dof(p01,1)] + x[lay1.dof(p00,1)] );
	else
		v0 = exact(test, 2, 0.5*(p01.x()+p00.x()), 0.5*(p01.y()+p00.y()) );
	GU[1] = (v1-v0)/(2*hx); // v_x
}

//void update_psi(const layout& lay1, const Sparse::Vector& x, const Sparse::Vector& x0, dof_vector<symmat2>& psi, dof_vector<symmat2>& psi0, dof_vector<symmat2>& psi1)
void update_psi(const layout& lay1, const dof_vector<double>& x, const dof_vector<double>& x0, dof_vector<symmat2>& psi, dof_vector<symmat2>& psi0, dof_vector<symmat2>& psi1)
{
	double time_rel = fabs(dt0) < 1.0e-12 ? 0.0 : dt/dt0;
	double a0 = (2*time_rel+1)/(time_rel+1), a1 = -(1+time_rel), a2 = time_rel*time_rel/(1+time_rel);
	double b1 = 1+time_rel, b2 = -time_rel;
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		mat2 OMEGA, OMEGA0, GU, GU0;
		symmat2 B, B0, Nr, Nr0, Nc, Nc0;
		const symmat2& PSI = psi(0,i,j);
		const symmat2& PSI0 = psi0(0,i,j);
		mat3 A = {0,0,0,0,0,0,0,0,0};
		double rhs[3] = {0,0,0}, sol[3];
		// time derivative: a0*phi(n+1) + a1*psi(n) + a2*psi(n-1)
		A[0+0*3] = A[1+1*3] = A[2+2*3] = a0;
		for(int q = 0; q < 3; ++q)
			rhs[q] -= (a1*PSI[q] + a2*PSI0[q]);
		fill_grad(lay1, x, i, j, GU);
		fill_grad(lay1, x0, i, j, GU0);
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
		flux_advection(lay1, psi, x, i, j, Nc);
		flux_advection(lay1, psi0, x0, i, j, Nc0);
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
			psi1(0,i,j)[q] = sol[q];
	}
	
	psi0 = psi;// psi(n-1) <- psi(n)
	psi = psi1;// psi(n)   <- psi(n+1)
}

void update_tau(const layout& lay1, const dof_vector<symmat2>& psi, dof_vector<symmat2>& tau1, dof_vector<symmat2>& tau)
{
	// remember tau before it is overwritten
	tau1 = tau;
	// update tau: tau = g(exp(psi)) / We * (exp(psi) -I)
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		symmat2 ePSI, eLMBD;
		mat2 R, Rt;
		double lmbd[2];
		// compute exp(psi)
		eigs(psi(0,i,j), lmbd, R);
		transpose(R, Rt);
		eLMBD[0] = exp(lmbd[0]); eLMBD[1] = exp(lmbd[1]); eLMBD[2] = 0.0;
		matmul3ns(Rt, eLMBD, ePSI);
		double g1 = g(ePSI);
		tau(0,i,j)[0] = g1/We * (ePSI[0] - 1);
		tau(0,i,j)[1] = g1/We * (ePSI[1] - 1);
		tau(0,i,j)[2] = g1/We * (ePSI[2] - 0);
	}
}

double residual(const dof_vector<symmat2>& tau0, const dof_vector<symmat2>& tau)
{
	double l2norm = 0.0;
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		double err = 0.0;
		for(int q = 0; q < 3; ++q)
			err = std::max(fabs(tau0(0,i,j)[q]-tau(0,i,j)[q]), err);
		l2norm += err*err*hx*hy;
	}
	return sqrt(l2norm);
}

int div_tau(int comp, int i, int j, std::array<pint,6>& dofs)
{
	if(comp == 1) // x
	{
		// xx
		dofs[0] = pint(3,i+1,j,3);
		dofs[1] = pint(3,i  ,j,3);
		dofs[0].coef =  1.0 / hx;
		dofs[1].coef = -1.0 / hx;
		// xy
		dofs[2] = pint(3,i  ,j+1,3);
		dofs[3] = pint(3,i+1,j+1,3);
		dofs[2].coef =  0.5 / hy;
		dofs[3].coef =  0.5 / hy;
		if(!dofs[2].inside() || !dofs[3].inside())
		{
			dofs[2] = pint(3,i  ,j,3);dofs[2].coef *= 2;
			dofs[3] = pint(3,i+1,j,3);dofs[3].coef *= 2;
		}
		dofs[4] = pint(3,i  ,j-1,3);
		dofs[5] = pint(3,i+1,j-1,3);
		dofs[4].coef = -0.5 / hy;
		dofs[5].coef = -0.5 / hy;
		if(!dofs[4].inside() || !dofs[5].inside())
		{
			dofs[4] = pint(3,i  ,j,3);dofs[4].coef *= 2;
			dofs[5] = pint(3,i+1,j,3);dofs[5].coef *= 2;
		}
	}
	else //if(comp == 2) // y
	{
		// yy
		dofs[0] = pint(3,i,j+1,3);
		dofs[1] = pint(3,i,j  ,3);
		dofs[0].coef =  1.0 / hy;
		dofs[1].coef = -1.0 / hy;
		// xy
		dofs[2] = pint(3,i+1,j  ,3);
		dofs[3] = pint(3,i+1,j+1,3);
		dofs[2].coef =  0.5 / hx;
		dofs[3].coef =  0.5 / hx;
		if(!dofs[2].inside() || !dofs[3].inside())
		{
			dofs[2] = pint(3,i,j  ,3);dofs[2].coef *= 2;
			dofs[3] = pint(3,i,j+1,3);dofs[3].coef *= 2;
		}
		dofs[4] = pint(3,i-1,j  ,3);
		dofs[5] = pint(3,i-1,j+1,3);
		dofs[4].coef = -0.5 / hx;
		dofs[5].coef = -0.5 / hx;
		if(!dofs[4].inside() || !dofs[5].inside())
		{
			dofs[4] = pint(3,i,j  ,3);dofs[4].coef *= 2;
			dofs[5] = pint(3,i,j+1,3);dofs[5].coef *= 2;
		}
	}
	return 6;
}

void fill_rhs(const layout& lay1, const dof_vector<symmat2>& psi, const dof_vector<symmat2>& tau, dof_vector<double>& b, double scale = 1.0)
{
	std::array<pint,6> dofs;
	int ndofs = 6;
	// grid1
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		int m = lay1.dof(pint(1,i,j),0);
		div_tau(1,i,j,dofs);
		for(int q = 0; q < ndofs/2; ++q)
			//if(dofs[q].inside())
			if( dofs[2*q].inside() && dofs[2*q+1].inside() )
			for(int k = 0; k < 2; ++k)
		{
			int q1 = 2*q+k;
			int d = q1 < 2 ? 0 : 2;//xx(0) xy(2)
			b[m] += dofs[q1].coef * mu_p*tau(0,dofs[q1].i,dofs[q1].j)[d] * scale;
		}
}
	// grid2
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		int m = lay1.dof(pint(2,i,j),1);
		div_tau(2,i,j,dofs);
		for(int q = 0; q < ndofs/2; ++q)
			//if(dofs[q].inside())
			if( dofs[2*q].inside() && dofs[2*q+1].inside() )
			for(int k = 0; k < 2; ++k)
		{
			int q1 = 2*q+k;
			int d = q1 < 2 ? 1 : 2;//yy(1) xy(2)
			b[m] += dofs[q1].coef * mu_p*tau(0,dofs[q1].i,dofs[q1].j)[d] * scale;
		}
	}
}

void fill_vector(const dof_vector<double>& in, Sparse::Vector& out)
{
	std::copy(in.U.begin(), in.U.end(), out.Begin());
}
void fill_vector(const Sparse::Vector& in, dof_vector<double>& out)
{
	std::copy(in.Begin(), in.End(), out.U.begin());
}
void fill_vector(const dof_vector<double>& in, std::vector<double>& out)
{
	std::copy(in.U.begin(), in.U.end(), out.begin());
}
void fill_vector(const std::vector<double>& in, dof_vector<double>& out)
{
	std::copy(in.begin(), in.end(), out.U.begin());
}

void save_vtk(std::string filename, const layout& lay1, const dof_vector<symmat2>& psi, const dof_vector<symmat2>& tau, const dof_vector<double>& x)
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
					divu += exact(test, q < 2 ? 1 : 2, dofs[q].x(), dofs[q].y()) * dofs[q].coef;
		}
		ofs << divu << '\n';
	}
	// tau
	const char* tau_name[3] = {"tau_xx","tau_yy","tau_xy"};
	for(int k = 0; k < 3; ++k)
	{
		ofs << "SCALARS " << tau_name[k] << " double 1\n";
		ofs << "LOOKUP_TABLE default\n";
		for(int j = 0; j < n2+1; ++j)
			for(int i = 0; i < n1+1; ++i)
		{
			double tau_ij = 0.0;
			if(i > 0 && i < n1 && j > 0 && j < n2)
				tau_ij = tau(0,i,j)[k];
			if(std::isnan(tau_ij) || std::isinf(tau_ij))
				tau_ij = -1e10;
			ofs << tau_ij << '\n';
		}
	}
	// psi
	const char* psi_name[3] = {"psi_xx","psi_yy","psi_xy"};
	for(int k = 0; k < 3; ++k)
	{
		ofs << "SCALARS " << psi_name[k] << " double 1\n";
		ofs << "LOOKUP_TABLE default\n";
		for(int j = 0; j < n2+1; ++j)
			for(int i = 0; i < n1+1; ++i)
		{
			double psi_ij = 0.0;
			if(i > 0 && i < n1 && j > 0 && j < n2)
				psi_ij = psi(0,i,j)[k];
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
			dofs[0] = pint(1,i,j  );
			dofs[1] = pint(1,i-1,j);
			for(int k = 0; k < 2; ++k)
				if(dofs[k].inside()) u[k] = dofs[k].coef = x[lay1.dof(dofs[k], 0)];
				else u[0] = u[1] = dofs[k].coef = exact(test, 1, dofs[k].x(), dofs[k].y());
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
				else u[0] = u[1] = dofs[k].coef = exact(test, 2, dofs[k].x(), dofs[k].y());
		}
		ofs << 0.5*(u[0]+u[1]) << '\n';
	}

	ofs.close();
}

double cfl_time_step(double t, double dt)
{
	double maxu = max_exact(test, 1);
	double h = std::min(hx,hy);
	double dt_cfl = cfl/(maxu+1.0e-12) * h;
	return dt_cfl;
}

void update_time()
{
	dt0 = dt;
	t += dt;
	dt = std::min(dt, cfl_time_step(t, dt));
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
void fill_stokes(dof_vector<double>& RHS0, matrix_type& A, const layout& lay1, double scale, bool filled)
{
	std::fill(RHS0.U.begin(),RHS0.U.end(),0.0);
	// equations for velocity
	pint dofs1[5], dofs2[5], dofs3[2];
	int ndofs;
	// equations for Uh
	for(int i = 0; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		int k = lay1.dof(pint(1,i,j), 0);
		laplace(i,j,1,dofs1,dofs2, -mu_s);
		RHS0[k] = rhs_exact(test, 1,(i+0.5)*hx,j*hy) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs1[q].inside() )
			{
				if(!filled)
				{
#if defined(USE_S3M)
					A.PushBack(lay1.dof(dofs1[q],dofs1[q].grid-1), dofs1[q].coef * scale);
#else
					A[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
#endif
				}
			}
			else //if(test == 0)
				RHS0[k] -= dofs1[q].coef * scale * exact(test, 1,dofs1[q].x(),dofs1[q].y());
		if(divgrad) // \Delta_h = 2 div_h D_h instead of five-point Laplace:
			for(int q = 0; q < 4; ++q)
				if( dofs2[q].inside() )
				{
					if(!filled)
					{
#if defined(USE_S3M)
						A.PushBack(lay1.dof(dofs2[q],dofs2[q].grid-1), dofs2[q].coef * scale);
#else
						A[k][lay1.dof(dofs2[q],dofs2[q].grid-1)] += dofs2[q].coef * scale;
#endif
					}
				}
				else //if(test == 0)
					RHS0[k] -= dofs2[q].coef * scale * exact(test, 2,dofs2[q].x(),dofs2[q].y());
		ndofs = grad(i,j,1,dofs3);
		if(test != 0 || (i != 0 && i != n1-1)) //skipping boundary gives BC dp/dx=0 (test = 0), some analytics may have dp/dx != 0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
				if(!filled)
				{
#if defined(USE_S3M)
					A.PushBack(lay1.dof(dofs3[q],dofs3[q].grid-1), dofs3[q].coef * scale);
#else
					A[k][lay1.dof(dofs3[q],dofs3[q].grid-1)] += dofs3[q].coef * scale;
#endif
				}
#if defined(USE_S3M)
		if(!filled) A.FinalizeRow();
#endif
	}
	// equations for Vh
	for(int j = 0; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		int k = lay1.dof(pint(2,i,j), 1);
		laplace(i,j,2,dofs1,dofs2, -mu_s);
		RHS0[k] = rhs_exact(test, 2,i*hx,(j+0.5)*hy) * scale;
		for(int q = 0; q < 5; ++q)
			if( dofs2[q].inside() )
			{
				if(!filled)
				{
#if defined(USE_S3M)
					A.PushBack(lay1.dof(dofs2[q],dofs2[q].grid-1), dofs2[q].coef * scale);
#else
					A[k][lay1.dof(dofs2[q],dofs2[q].grid-1)] += dofs2[q].coef * scale;
#endif
				}

			}
			else //if(test == 0)
				RHS0[k] -= dofs2[q].coef * scale * exact(test, 2,dofs2[q].x(),dofs2[q].y());
		if(divgrad) // \Delta_h = 2 div_h D_h instead of five-point Laplace:
			for(int q = 0; q < 4; ++q)
				if( dofs1[q].inside() )
				{
					if(!filled)
					{
#if defined(USE_S3M)
						A.PushBack(lay1.dof(dofs1[q],dofs1[q].grid-1), dofs1[q].coef * scale);
#else
						A[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
#endif
					}
				}
				else //if(test == 0)
					RHS0[k] -= dofs1[q].coef * scale * exact(test, 1,dofs1[q].x(),dofs1[q].y());
		ndofs = grad(i,j,2,dofs3);
		if(test != 0 || (j != 0 && j != n2-1))//skipping boundary gives BC dp/dy=0 (test = 0), some analytics may have dp/dy != 0
			for(int q = 0; q < ndofs; ++q) if( dofs3[q].inside() )
				if(!filled)
				{
#if defined(USE_S3M)
					A.PushBack(lay1.dof(dofs3[q],dofs3[q].grid-1), dofs3[q].coef * scale);
#else
					A[k][lay1.dof(dofs3[q],dofs3[q].grid-1)] += dofs3[q].coef * scale;
#endif
				}
#if defined(USE_S3M)
		if(!filled) A.FinalizeRow();
#endif
	}
	// div u = 0, equations for Ph
	for(int i = 1; i < n1; ++i)
		for(int j = 1; j < n2; ++j)
	{
		int k = lay1.dof(pint(3,i,j), 2);
		ndofs = div(i,j,dofs1);
		/*if(integral_constraint && !filled)
		{
			A[k][nstokes-1] = scale;
			A[nstokes-1][k] = scale;
		}*/
		for(int q = 0; q < ndofs; ++q)
			if( dofs1[q].inside() )
			{
				if(!filled)
				{
#if defined(USE_S3M)
					A.PushBack(lay1.dof(dofs1[q],dofs1[q].grid-1), dofs1[q].coef * scale);
#else
					A[k][lay1.dof(dofs1[q],dofs1[q].grid-1)] += dofs1[q].coef * scale;
#endif
				}
			}
			else //if(test == 0)
				RHS0[k] -= dofs1[q].coef * scale * exact(test, q < 2 ? 1 : 2,dofs1[q].x(),dofs1[q].y());
#if defined(USE_S3M)
		if(!filled) A.FinalizeRow();
#endif
	}
}


int main(int argc, char ** argv)
{
#if !defined(USE_S3M)
	Solver::Initialize(&argc, &argv, "database.xml");
#endif
	if(argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " n1 [T=" << T
			<< "] [nframe=" << nframe
			<< "] [We=" << We
			<< "]\n";
		return 0;
	}
	
	n1 = atoi(argv[1]);
	n2 = n1;//argv[2];
	hx = 1.0/n1;
	hy = 1.0/n2;
	if(argc > 2)
		T = atof(argv[2]);
	if(argc > 3)
		nframe = atoi(argv[3]);
	if(argc > 4)
		We = atof(argv[4]);
	dt = cfl * std::min(hx,hy);
	std::cout << "n1 " << n1 << " n2 " << n2 << " test " << test << " dt " << dt << std::endl;
	std::cout << "nframe " << nframe << " T " << T << " We " << We << std::endl;
	
	const int grids1[3] = {1,2,3};	// Uh(1), Vh(2), Ph(3)
	const int ncomps1[3] = {1,1,1};
	const int grids2[1] = {3};	// Psi(3), Psi=(Psi_xx, Psi_yy, Psi_xy)
	const int ncomps2[1] = {3};
	layout lay1(3,grids1,ncomps1), lay2(1,grids2,ncomps2);
	dof_vector<double> UVP(&lay1), UVP0(&lay1), RHS(&lay1), RHS0(&lay1);
	dof_vector<symmat2> PSI(&lay2), PSI0(&lay2), PSI1(&lay2), TAU(&lay2), TAU0(&lay2), TAU1(&lay2);

	int nstokes = lay1.offset(3);//n1*(n2-1)+n2*(n1-1)+(n1-1)*(n2-1); // Uh,Vh,Ph
	int ntau = lay2.offset(1);//3*(n1-1)*(n2-1);
	std::cout << "nstokes " << nstokes << " ntau " << ntau << std::endl;
	// psi = (psi_xx (Ph), psi_yy (Ph), psi_xy (Ph))
	// psi -- psi(n), psi0 -- psi(n-1), psi1 -- psi(n+1)
	// tau -- tau(n+1), tau0 -- tau(n-1), tau1 -- tau(n)

	//bool integral_constraint = false;
	//if(integral_constraint) nstokes += 1; // integral constraint on pressure

	bool filled = false;
	bool success = false, finish = false;
	double res = 0.0, atol = 1.0e-4;
	int iter = 0, maxiters = 10000;
#if defined(USE_S3M)
	//BramblePasciakCG< AMGRugeStuben< GaussSeidel, BICGSTAB<ILDUC> > > S;
	BICGSTAB< Vanka > S;
	S.GetParameters().Load("params.txt");
	S.GetParameters().Set<int>("B_block_end", lay1.offset(2));//Uh+Vh
	CSRMatrix A;
	A.ReserveSize(nstokes);
	std::vector<double> x(nstokes), b(nstokes);
#else
	Sparse::Matrix A("stokes", 0, nstokes);
	Sparse::Vector b("stokes", 0, nstokes), x("stokes", 0, nstokes);//, b0("stokes", 0, nstokes), x0("stokes", 0, nstokes);
	Solver S(Solver::INNER_ILU2, "test");
	//Solver S(Solver::PETSc, "test");
#endif	
	do
	{
		double scale = hx*hy;
		fill_stokes(RHS0, A, lay1, scale, filled);
		if(!filled)
		{
#if defined(USE_S3M)
			A.SortRows();
			A.Save("A.mtx");
			S.Setup(A);
#else
			A.Save("A.mtx");
			S.SetMatrix(A);
#endif
			filled = true;
		}
		//b0.Save("b0.txt");
		//std::copy(b0.Begin(),b0.End(),b.Begin());
		RHS = RHS0;
		if(test == 0) fill_rhs(lay1, PSI, TAU, RHS, scale);
		fill_vector(RHS, b);
#if !defined(USE_S3M)
		std::fill(x.Begin(),x.End(),0.0);
#else
		std::fill(x.begin(),x.end(),0.0);
#endif
		success = S.Solve(b,x);
		fill_vector(x, UVP);
#if !defined(USE_S3M)
		b.Save("b.txt");
		x.Save("x.txt");
		int iters = S.Iterations();
		double resid = S.Residual();
		std::string reason = S.ReturnReason();
		std::cout << (!success ? "not " : "") << "converged iters " << iters << " resid " << resid << " reason " << reason << std::endl;
#else
		double resid = Resid(A,b,x);
		std::cout << (!success ? "not " : "") << "converged resid " << resid << std::endl;
#endif
		if(test)
			finish = true;
		else if(success)
		{
			//update_psi(lay1, x, x0, PSI, PSI0, PSI1);
			update_psi(lay1, UVP, UVP0, PSI, PSI0, PSI1);
			//std::copy(UVP.U.begin(), UVP.U.end(), UVP0.U.begin());
			UVP0 = UVP;
			//update_tau(lay1, PSI, TAU1, TAU);
			update_tau(lay1, PSI, TAU1, TAU);
			res = residual(TAU0, TAU);
			TAU0 = TAU1;
			std::cout << "t " << t << "/" << T << " iter " << iter << " " << " resid " << res << std::endl;
			//std::cin.ignore();
			if(iter % nframe == 0)
			{
				int frame = iter/nframe;
				save_vtk("out" + std::to_string(frame) + ".vtk", lay1, PSI, TAU, UVP);
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
