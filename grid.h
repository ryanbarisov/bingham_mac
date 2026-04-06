#ifndef _GRID_H
#define _GRID_H

#include <cassert>
#include <cmath>
#include <vector>

// DOF on staggered grid
typedef struct pint
{
	int i, j;
	double coef = 0.0;
	int grid; // 1(Uh), 2(Vh), 3(Ph), 4(Qh)
	int level = 0; // multigrid level; level = 0 -- finest grid
	static int glev; // finest grid has 2^glev+1 nodes in each direction

	pint(int grid = -1, int i = -10, int j = -10, int level = 0) : i(i), j(j), grid(grid), level(level) {}
	static int nx(int grid, int level) // amount of unknowns in X direction depends on level
	{
		if(grid == 1 || grid == 4) return (1 << (glev-level)) - 1;
		else return (1 << (glev-level)) + 0;
	}
	static int ny(int grid, int level) // amount of unknowns in Y direction depends on level
	{
		if(grid == 2 || grid == 4) return (1 << (glev-level)) - 1;
		else return (1 << (glev-level)) + 0;
	}
	static double hx(int level) { return 1.0/nx(3,level); }
	static double hy(int level) { return 1.0/ny(3,level); }
	static int ndofs(int grid, int level) { return nx(grid,level)*ny(grid,level); }
	int dof() const
	{
		//if(grid == 1) // Uh
		//	return i*ny(grid,level) + (j-1);
		if(grid == 2) // Vh
			return (j-1)*nx(grid,level) + (i-1);
		//else if(grid == 3) // Ph
		else
			return (i-1)*ny(grid,level) + (j-1);
		//else // if(grid == 4) // Qh
		//	return i*ny(grid,level) + j;
	}
	bool inside() const
	{
		int i1 = i-1, j1 = j-1;
		//if(grid == 1) i1--;
		//else if(grid == 2) j1--;
		//else if(grid == 3) i1--, j1--;
		return i1 >= 0 && i1 < nx(grid,level) && j1 >= 0 && j1 < ny(grid,level);
	}
	double x() const
	{
		if(grid == 1 || grid == 4)
			return std::min(std::max(0.0, i*hx(level)), 1.0);
		else
			return std::min(std::max(0.0, (i-0.5)*hx(level)), 1.0);
	}
	double y() const
	{
		if(grid == 2 || grid == 4)
			return std::min(std::max(0.0, j*hy(level)), 1.0);
		else
			return std::min(std::max(0.0, (j-0.5)*hy(level)), 1.0);
	}
	bool operator==(const pint &other)
	{
		return level == other.level && grid == other.grid && i == other.i && j == other.j;
	}
} pint;
int pint::glev = -1;


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
	void zero(int level)
	{
		std::fill(U.begin()+offset(level), U.begin()+offset(level+1), 0.0);
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

	double l2norm(int level = 0) const
	{
		double norm = 0.0;
		double hx = pint::hx(level), hy = pint::hy(level);
		for(int k = offset(level); k < offset(level+1); ++k)
			norm += U[k]*U[k]*hx*hy;
		return sqrt(norm);
	}
};


#endif //_GRID_H
