#ifndef _VTK_H
#define _VTK_H

#if 0
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
#else
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
	for(int j = 1; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
			ofs << x(2,pint(3,i,j,level)) << '\n';
	//divergence
	ofs << "SCALARS div_U double 1\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 1; j < n2; ++j)
		for(int i = 1; i < n1; ++i)
	{
		double divu = 0.0;
		pint dofs[8], dof = pint(3,i,j,level);
		int ndofs = div(dof,dofs);
		for(int q = 0; q < ndofs; ++q)
			if( dofs[q].inside() )
				divu += x(dofs[q].grid-1,dofs[q]) * dofs[q].coef;
			else
				divu += exact(test, q < 2 ? 1 : 2, dofs[q].x(), dofs[q].y()) * dofs[q].coef;
		ofs << divu << '\n';
	}
	// tau
	const char* tau_name[3] = {"tau_xx","tau_yy","tau_xy"};
	if(0) for(int k = 0; k < 3; ++k) if(level == 0)
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
	if(0) for(int k = 0; k < 3; ++k) if(level == 0)
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
	ofs << "SCALARS U double 2\n";
	ofs << "LOOKUP_TABLE default\n";
	for(int j = 1; j < n2; ++j)
		for(int i = 0; i < n1-1; ++i)
		{
			for(int q = 0; q < 2; ++q)
			{
				pint pdof(1,i+q,j  ,level);
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
	for(int j = 0; j < n2-1; ++j)
		for(int i = 1; i < n1; ++i)
		{
			for(int q = 0; q < 2; ++q)
			{
				pint pdof(2,i ,j+q,level);
				double v;
				if(pdof.inside()) v = x(1, pdof);
				else v = exact(test, 2, pdof.x(), pdof.y());
				ofs << v << ' ';
			}
			ofs << '\n';
		}
	ofs.close();
}
#endif

#endif
