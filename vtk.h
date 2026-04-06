#ifndef _VTK_H
#define _VTK_H

#include "mat2.h"
#include "discr.h"

void save_vtk(std::string filename, const dof_vector<symmat2>& psi, const dof_vector<symmat2>& tau, const dof_vector<double>& x, int level = 0)
{
	std::ofstream ofs(filename);

	int nx1 = pint::nx(1,level), ny1 = pint::ny(1,level);
	int nx2 = pint::nx(2,level), ny2 = pint::ny(2,level);
	int nx3 = pint::nx(3,level), ny3 = pint::ny(3,level);
	int n1 = nx3+1, n2 = ny3+1; // amount of points
	//version,identifier
	ofs << "# vtk DataFile Version 3.0\n";
	//header
	ofs << "Bingham flow rectangular grid " << n1 << " " << n2 << '\n';
	//file format
	ofs << "ASCII\n";
	//dataset type
	ofs << "DATASET RECTILINEAR_GRID\n";
	ofs << "DIMENSIONS " << n1 << " " << n2 << " 1" << '\n';
	ofs << "X_COORDINATES " << n1 << " double\n";
	for(int i = 0; i < n1; ++i)
	{
		pint pdof(1,i,1,level);
		ofs << pdof.x() << ' ';
	}
	ofs << '\n';
	ofs << "Y_COORDINATES " << n2 << " double\n";
	for(int j = 0; j < n2; ++j)
	{
		pint pdof(2,1,j,level);
		ofs << pdof.y() << ' ';
	}
	ofs << '\n';
	ofs << "Z_COORDINATES " << 1 << " double\n";
	ofs << "0.0\n";

	//data
	ofs << "CELL_DATA " << ((n1-1)*(n2-1)) << '\n';
	//pressure
	ofs << "SCALARS pressure double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 1; j < ny3+1; ++j)
		for(int i = 1; i < nx3+1; ++i)
			ofs << x(2,pint(3,i,j,level)) << '\n';
	//divergence
	ofs << "SCALARS div_U double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 1; j < ny3+1; ++j)
		for(int i = 1; i < nx3+1; ++i)
	{
		double divu = 0.0;
		pint dofs[8], dof = pint(3,i,j,level);
		int ndofs = div(dof,dofs);
		for(int q = 0; q < ndofs; ++q)
			if( dofs[q].inside() )
				divu += x(dofs[q].grid-1,dofs[q]) * dofs[q].coef;
			else
				divu += exact(test, dofs[q].grid, dofs[q].x(), dofs[q].y()) * dofs[q].coef;
		ofs << divu << '\n';
	}
	// tau
	const char* tau_name[3] = {"tau_xx","tau_yy","tau_xy"};
	for(int k = 0; k < 3; ++k) if(level == 0)
	{
		ofs << "SCALARS " << tau_name[k] << " double 1\n";
		ofs << "LOOKUP_TABLE default\n";
		for(int j = 1; j < ny3+1; ++j)
			for(int i = 1; i < nx3+1; ++i)
		{
			pint pdof(3,i,j,level);
			double tau_ij = tau(0,pdof)[k];
			if(std::isnan(tau_ij) || std::isinf(tau_ij)) tau_ij = -1e10;
			ofs << tau_ij << '\n';
		}
	}
	// psi
	const char* psi_name[3] = {"psi_xx","psi_yy","psi_xy"};
	for(int k = 0; k < 3; ++k) if(level == 0)
	{
		ofs << "SCALARS " << psi_name[k] << " double 1\n";
		ofs << "LOOKUP_TABLE default\n";
		for(int j = 1; j < ny3+1; ++j)
			for(int i = 1; i < nx3+1; ++i)
		{
			pint pdof(3,i,j,level);
			double psi_ij = psi(0,pdof)[k];
			if(std::isnan(psi_ij) || std::isinf(psi_ij)) psi_ij = -1e10;
			ofs << psi_ij << '\n';
		}
	}
	// u
	ofs << "SCALARS U double 2\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 1; j < ny3+1; ++j)
		for(int i = 1; i < nx3+1; ++i)
		{
			for(int q = 0; q < 2; ++q)
			{
				pint pdof(1,i-1+q,j  ,level);
				double u;
				if(pdof.inside()) u = x(0, pdof);
				else u = exact(test, 1, pdof.x(), pdof.y());
				ofs << u << ' ';
			}
			ofs << '\n';
		}
	// v
	ofs << "SCALARS V double 2\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 1; j < ny3+1; ++j)
		for(int i = 1; i < nx3+1; ++i)
		{
			for(int q = 0; q < 2; ++q)
			{
				pint pdof(2,i ,j-1+q,level);
				double v;
				if(pdof.inside()) v = x(1, pdof);
				else v = exact(test, 2, pdof.x(), pdof.y());
				ofs << v << ' ';
			}
			ofs << '\n';
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


void save_vtk_poisson(std::string filename, const dof_vector<double>& x, int level = 0)
{
	std::ofstream ofs(filename);

	int nx1 = pint::nx(1,level), ny1 = pint::ny(1,level);
	int nx2 = pint::nx(2,level), ny2 = pint::ny(2,level);
	int nx3 = pint::nx(3,level), ny3 = pint::ny(3,level);
	int n1 = nx3+1, n2 = ny3+1; // amount of points
	//version,identifier
	ofs << "# vtk DataFile Version 3.0\n";
	//header
	ofs << "Bingham flow rectangular grid " << n1 << " " << n2 << '\n';
	//file format
	ofs << "ASCII\n";
	//dataset type
	ofs << "DATASET RECTILINEAR_GRID\n";
	ofs << "DIMENSIONS " << n1 << " " << n2 << " 1" << '\n';
	ofs << "X_COORDINATES " << n1 << " double\n";
	for(int i = 0; i < n1; ++i)
	{
		pint pdof(1,i,1,level);
		ofs << pdof.x() << ' ';
	}
	ofs << '\n';
	ofs << "Y_COORDINATES " << n2 << " double\n";
	for(int j = 0; j < n2; ++j)
	{
		pint pdof(2,1,j,level);
		ofs << pdof.y() << ' ';
	}
	ofs << '\n';
	ofs << "Z_COORDINATES " << 1 << " double\n";
	ofs << "0.0\n";

	//data
	ofs << "CELL_DATA " << ((n1-1)*(n2-1)) << '\n';
	//pressure
	ofs << "SCALARS pressure double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int i = 1; i < nx3+1; ++i)
		for(int j = 1; j < ny3+1; ++j)
			ofs << x(0,pint(3,i,j,level)) << '\n';
	ofs.close();
}

#endif
