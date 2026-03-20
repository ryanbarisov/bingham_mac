#include "analytics.h"
#include <cmath>

extern double mu_s;
extern double t;

// exact solution for analytical tests, including velocity BC: u,v,p
double exact(int test, int grid, double x, double y)
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
	else // if(test == 0)
	{
		// no analytical solution, it is used to setup velocity BC
		if(grid == 1) // Uh
			//return fabs(y-1) < 1.0e-6 ? 16*x*x*(1-x)*(1-x) : 0.0;
			return fabs(y-1) < 1.0e-6 ? 8*(1+tanh(8*(t-0.5)))*x*x*(1-x)*(1-x) : 0.0;
		else if(grid == 2) // Vh
			return 0.0;
		else // if(grid == 3) // Ph
			return 0.0;
	}
}

double max_exact(int test, int grid)
{
	if(test == 1)
	{
		if(grid == 1 || grid == 2) // Uh,Vh
			return  1.0/(4*M_PI*M_PI) * sqrt(2);
		else // if(grid == 3) // Ph
		     	return 1.0/M_PI;
	}
	else
	{
		if(grid == 1 || grid == 2) // Uh,Vh
			return 8*(1+tanh(t-0.5))/16;
		else
			return 1.0;
	}
}

// analytical RHS for Stokes substep
double rhs_exact(int test, int grid, double x, double y)
{
	if(test == 1)
	{
		if(grid == 1) // Uh
			return -mu_s*( 2*cos(2*M_PI*x)-1)*sin(2*M_PI*y) + 2*cos(2*M_PI*x)*sin(2*M_PI*y);
		else //if(grid == 2) // Vh
			return -mu_s*(-2*cos(2*M_PI*y)+1)*sin(2*M_PI*x) + 2*cos(2*M_PI*y)*sin(2*M_PI*x);
	}
	else // if(test == 0)
	{
		return 0.0;
	}
}
