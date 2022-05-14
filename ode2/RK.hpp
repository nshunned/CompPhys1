// generic solvers for linear first order ODE
// for convenience these function return a graph of y(x)

#include <TGraph.h>

// generic fcn pointer for a linear differential equation, f = y'(x,y)
// double (*f)(double x, double y)
// the equation may depend on one independent variable (x) and a
// dependent variable (y)
// y0: initial condition
// nsteps, x0, xmax are used to se t the range and h, the step size
TGraph RK1Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax);


TGraph RK2Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax);

TGraph RK4Solve(double (*f)(double x, double y), double y0,
                int nsteps, double x0, double xmax);
