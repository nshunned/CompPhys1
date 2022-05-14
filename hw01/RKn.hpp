// generic solvers for linear first order ODE
// for convenience these function return a graph of the solution for y(x)
#pragma once 
#include <TGraph.h>
#include <vector>
using std::vector;

/**********************************************
RK4 solver for set of N 1st order ODEs (1 independent parameter)
The list of eqns is stored in a vector of fcn pointers
initial conditions are given in a vector of doubles

B/c the 1st order ODEs may be coupled, each equation will have 
access to the entire set of dependent parameters, thus we use 
the interface: double func_t(double X, const vector<double> &Y),
where X represents the independent paramater and Y is a vector of
the dependent parameters 

The solver will return a set of graphs for each dependent variable
***********************************************/
typedef double func_t(double, const vector<double>&);
typedef func_t* pfunc_t;   // function pointer matching our ODE interface

vector<double> RK4SolveNAProblem_1(vector<pfunc_t> &fnlist, vector<double> &y0, int nsteps, double x0, double xmax, pfunc_t fstop=0, double errdef=1e-9, int maxrep=5);

vector<TGraph> RK4SolveNAProblem_2(vector<pfunc_t> &fnlist, vector<double> &y0, int nsteps, double x0, double xmax, pfunc_t fstop=0, double errdef=1e-9, int maxrep=5);