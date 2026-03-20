#ifndef _ANALYTICS_H
#define _ANALYTICS_H

// exact solution for analytical tests, including velocity BC: u,v,p
double exact(int test, int grid, double x, double y);
// analytical RHS for Stokes substep
double rhs_exact(int test, int grid, double x, double y);
// for CFL:
double max_exact(int test, int grid);

#endif
