#ifndef _DISCR_H
#define _DISCR_H

#include "analytics.h"
#include "grid.h"

extern double We;
extern int test;
extern double dt, dt0;
extern bool divgrad;
extern int level0;

void laplace(const pint& pdof, pint * dofs1, pint * dofs2, double mult)
{
	int i = pdof.i, j = pdof.j, grid = pdof.grid, level = pdof.level;
	double hx = pint::hx(level), hy = pint::hy(level);
	if(grid == 1) // \Delta_h, x component
	{
		// Uh
		double h[4] = {hx,hx,hy,hy};
		dofs1[0] = pint(1,i+1,j,level);
		dofs1[1] = pint(1,i-1,j,level);
		dofs1[0].coef = 1./(hx*hx) * mult;
		dofs1[1].coef = 1./(hx*hx) * mult;
		dofs1[2] = pint(1,i,j+1,level); if(!dofs1[2].inside()) h[2] = 0.5*hy;
		dofs1[3] = pint(1,i,j-1,level); if(!dofs1[3].inside()) h[3] = 0.5*hy;
		dofs1[2].coef = 2./(h[2]*(h[2]+h[3])) * mult;
		dofs1[3].coef = 2./(h[3]*(h[2]+h[3])) * mult;
		dofs1[4] = pint(1,i,j  ,level); dofs1[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs1[4].coef -= dofs1[k].coef;
	}
	else if(grid == 2) // \Delta_h, y component
	{
		// Vh
		double h[4] = {hy,hy,hx,hx};
		dofs2[0] = pint(2,i,j+1,level);
		dofs2[1] = pint(2,i,j-1,level);
		dofs2[0].coef = 1./(hy*hy) * mult;
		dofs2[1].coef = 1./(hy*hy) * mult;
		dofs2[2] = pint(2,i+1,j,level); if(!dofs2[2].inside()) h[2] = 0.5*hx;
		dofs2[3] = pint(2,i-1,j,level); if(!dofs2[3].inside()) h[3] = 0.5*hx;
		dofs2[2].coef = 2./(h[2]*(h[2]+h[3])) * mult;
		dofs2[3].coef = 2./(h[3]*(h[2]+h[3])) * mult;
		dofs2[4] = pint(2,i,j  ,level); dofs2[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs2[4].coef -= dofs2[k].coef;
	}
	else if(grid == 3) // \Delta_h on Ph
	{
		double h[4] = {hx,hx,hy,hy};
		dofs1[0] = pint(3,i+1,j,level); if(!dofs1[0].inside()) h[0] = 0.5*hx;
		dofs1[1] = pint(3,i-1,j,level); if(!dofs1[1].inside()) h[1] = 0.5*hx;
		dofs1[0].coef = 2./(h[0]*(h[0]+h[1])) * mult;
		dofs1[1].coef = 2./(h[1]*(h[0]+h[1])) * mult;
		dofs1[2] = pint(3,i,j+1,level); if(!dofs1[2].inside()) h[2] = 0.5*hy;
		dofs1[3] = pint(3,i,j-1,level); if(!dofs1[3].inside()) h[3] = 0.5*hy;
		dofs1[2].coef = 2./(h[2]*(h[2]+h[3])) * mult;
		dofs1[3].coef = 2./(h[3]*(h[2]+h[3])) * mult;
		dofs1[4] = pint(3,i,j  ,level); dofs1[4].coef = 0.0;
		for(int k = 0; k < 4; ++k) dofs1[4].coef -= dofs1[k].coef;
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

// dp/dn = 0 on boundaries, needs review if otherwise
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
	else //if(grid == 2) // y
	{
		assert(grid == 2);
		dofs3[0] = pint(3,i,j+1,level); dofs3[0].coef =  1.0 / hy;
		dofs3[1] = pint(3,i,j  ,level); dofs3[1].coef = -1.0 / hy;
	}
	return 2;
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
		bool one_sided = false;
		if(!dofs[2].inside() || !dofs[3].inside())
		{
			dofs[2] = pint(3,i  ,j,level);
			dofs[3] = pint(3,i+1,j,level);
			one_sided = true;
		}
		dofs[4] = pint(3,i  ,j-1,level);
		dofs[5] = pint(3,i+1,j-1,level);
		dofs[4].coef = -0.5 / hy;
		dofs[5].coef = -0.5 / hy;
		if(!dofs[4].inside() || !dofs[5].inside())
		{
			dofs[4] = pint(3,i  ,j,level); dofs[4].coef *= 2;
			dofs[5] = pint(3,i+1,j,level); dofs[5].coef *= 2;
			one_sided = true;
		}
		if(one_sided)
		{
			dofs[2].coef = dofs[3].coef =  1.0 / hy;
			dofs[4].coef = dofs[5].coef = -1.0 / hy;
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
		bool one_sided = false;
		if(!dofs[2].inside() || !dofs[3].inside())
		{
			dofs[2] = pint(3,i,j  ,level);
			dofs[3] = pint(3,i,j+1,level);
			one_sided = true;
		}
		dofs[4] = pint(3,i-1,j  ,level);
		dofs[5] = pint(3,i-1,j+1,level);
		dofs[4].coef = -0.5 / hx;
		dofs[5].coef = -0.5 / hx;
		if(!dofs[4].inside() || !dofs[5].inside())
		{
			dofs[4] = pint(3,i,j  ,level);
			dofs[5] = pint(3,i,j+1,level);
			one_sided = true;
		}
		if(one_sided)
		{
			dofs[2].coef = dofs[3].coef =  1.0 / hx;
			dofs[4].coef = dofs[5].coef = -1.0 / hx;
		}
	}
	return 6;
}

void fill_rhs(const dof_vector<symmat2>& tau, dof_vector<double>& b, double mult)
{
	std::array<pint,6> dofs;
	int ndofs = 6;
	int level = 0; // TODO
	int nx1 = pint::nx(1,level), ny1 = pint::ny(1,level);
	int nx2 = pint::nx(2,level), ny2 = pint::ny(2,level);
	double scale = pint::hx(level)*pint::hy(level);
	// grid1
	for(int i = 1; i < nx1+1; ++i)
		for(int j = 1; j < ny1+1; ++j)
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
			b(0,pdof) += dofs[q1].coef * mult*tau(0,dofs[q1])[d] * scale;
		}
	}
	// grid2
	for(int j = 1; j < ny2+1; ++j)
		for(int i = 1; i < nx2+1; ++i)
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
			b(1,pdof) += dofs[q1].coef * mult*tau(0,dofs[q1])[d] * scale;
		}
	}
}



static double minmod(double a, double b) {return 0.5*((2*(a>0)-1)+(2*(b>0)-1)) * std::min(fabs(a),fabs(b));}

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
	u1 = pu1.inside() ? x(0,pu1) : exact(test,1,pu1.x(),pu1.y());//x(0,i,j,level);
	for(int q = 0; q < 3; ++q)
		H1[q] = u1*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(u1)*(PSI_PLUS[q] - PSI_MINUS[q]);
	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 1); // -x
	u0 = pu0.inside() ? x(0,pu0) : exact(test,1,pu0.x(),pu0.y());//x(0,i-1,j,level);
	for(int q = 0; q < 3; ++q)
		H0[q] = u0*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(u0)*(PSI_PLUS[q] - PSI_MINUS[q]);
	for(int q = 0; q < 3; ++q)
		Nc[q] += (H1[q]-H0[q])/hx;

	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 2); // +y
	v1 = pv1.inside() ? x(1,pv1) : exact(test,2,pv1.x(),pv1.y());//x(1,i,j,level);
	for(int q = 0; q < 3; ++q)
		H1[q] = v1*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(v1)*(PSI_PLUS[q] - PSI_MINUS[q]);
	reconstruct(psi, pdof, PSI_PLUS, PSI_MINUS, 3); // -y
	v0 = pv0.inside() ? x(1,pv0) : exact(test,2,pv0.x(),pv0.y());//x(1,i,j-1,level);
	for(int q = 0; q < 3; ++q)
		H0[q] = v0*0.5*(PSI_PLUS[q]+PSI_MINUS[q]) - c * fabs(v0)*(PSI_PLUS[q] - PSI_MINUS[q]);
	for(int q = 0; q < 3; ++q)
		Nc[q] += (H1[q]-H0[q])/hy;
}

//OLDROYD-B:
static double g(const symmat2& PSI) {return 1.0;}
static void P(const symmat2& PSI, symmat2& vP)
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
	double hx = pint::hx(level), hy = pint::hy(level), h;
	u1 = p11.inside() ? U(0,p11) : exact(test,1,p11.x(),p11.y());
	u0 = p10.inside() ? U(0,p10) : exact(test,1,p10.x(),p10.y());
	v1 = p01.inside() ? U(1,p01) : 0;//exact(test,2,p01.x(),p01.y());
	v0 = p00.inside() ? U(1,p00) : exact(test,2,p00.x(),p00.y());
	GU[0] = (u1-u0)/hx; //u_x
	GU[3] = (v1-v0)/hy; //v_y
	//u_y
	p11 = pint(1,i  ,j+1,level);
	p10 = pint(1,i-1,j+1,level);
	p01 = pint(1,i  ,j-1,level);
	p00 = pint(1,i-1,j-1,level);

	h = 2*hy;
	if(p10.inside() && p11.inside())
		u1 = 0.5*( U(0,p11) + U(0,p10) );
	else
	{
		u1 = exact(test, 1, 0.5*(p11.x()+p10.x()), 0.5*(p11.y()+p10.y()) );
		h = 1.5*hy;
	}
	if(p00.inside() && p01.inside())
		u0 = 0.5*( U(0,p01) + U(0,p00) );
	else
	{
		u0 = exact(test, 1, 0.5*(p01.x()+p00.x()), 0.5*(p01.y()+p00.y()) );
		h = 1.5*hy;
	}
	GU[2] = (u1-u0)/h; // u_y
	//v_x
	p11 = pint(2,i+1,j  ,level);
	p10 = pint(2,i+1,j-1,level);
	p01 = pint(2,i-1,j  ,level);
	p00 = pint(2,i-1,j-1,level);
	h = 2*hx;
	if(p10.inside() && p11.inside())
		v1 = 0.5*( U(1,p11) + U(1,p10) );
	else
	{
		v1 = exact(test, 2, 0.5*(p11.x()+p10.x()), 0.5*(p11.y()+p10.y()) );
		h = 1.5*hx;
	}
	if(p00.inside() && p01.inside())
		v0 = 0.5*( U(1,p01) + U(1,p00) );
	else
	{
		v0 = exact(test, 2, 0.5*(p01.x()+p00.x()), 0.5*(p01.y()+p00.y()) );
		h = 1.5*hx;
	}
	GU[1] = (v1-v0)/h; // v_x
}

void update_psi(const dof_vector<double>& x, const dof_vector<double>& x0, dof_vector<symmat2>& psi, dof_vector<symmat2>& psi0, dof_vector<symmat2>& psi1)
{
	double time_rel = fabs(dt0) < 1.0e-12 ? 0.0 : dt/dt0;
	double a0 = (2*time_rel+1)/(time_rel+1), a1 = -(1+time_rel), a2 = time_rel*time_rel/(1+time_rel);
	double b1 = 1+time_rel, b2 = -time_rel;
	int level = 0; // TODO
	int nx3 = pint::nx(3,level), ny3 = pint::ny(3,level);
	for(int i = 1; i < nx3+1; ++i)
		for(int j = 1; j < ny3+1; ++j)
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
	int nx3 = pint::nx(3,level), ny3 = pint::ny(3,level);
	// update tau: tau = g(exp(psi)) / We * (exp(psi) -I)
	for(int i = 1; i < nx3+1; ++i)
		for(int j = 1; j < ny3+1; ++j)
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
	int nx3 = pint::nx(3,level), ny3 = pint::ny(3,level);
	double hx = pint::hx(level), hy = pint::hy(level);
	double l2norm = 0.0;
	for(int i = 1; i < nx3+1; ++i)
		for(int j = 1; j < ny3+1; ++j)
	{
		pint pdof(3,i,j,level);
		double err = 0.0;
		for(int q = 0; q < 3; ++q)
			err = std::max(fabs(tau0(0,pdof)[q]-tau(0,pdof)[q]), err);
		l2norm += err*err*hx*hy;
	}
	return sqrt(l2norm);
}

#endif //_DISCR_H
