#include "RKn.hpp"
#include <assert.h>
#include <algorithm>
#include <cmath>

#include <stdio.h>

using namespace std;

// called by RK4SolveN, this function performs one step using the RK4 algorithm and returns the new vector of y values
vector<double> RK4StepN(vector<pfunc_t> &f, vector<double> y, double x, double h){
  int nFcn=f.size();
  vector<double> k1(nFcn), k2(nFcn), k3(nFcn), k4(nFcn);
  vector<double> ytmp(nFcn);

  for (int i=0; i<nFcn; i++){
    k1[i] = h * f[i](x, y);
    ytmp[i] = y[i] + k1[i]/2;
  }
  for (int i=0; i<nFcn; i++){
    k2[i] = h * f[i](x+h/2, ytmp);
    ytmp[i] = y[i] + k2[i]/2;
  }
  for (int i=0; i<nFcn; i++){
    k3[i] = h * f[i](x+h/2, ytmp);
    ytmp[i] = y[i] + k3[i];
  }
  for (int i=0; i<nFcn; i++){
    k4[i] = h * f[i](x+h, ytmp);
  }
  // calculate next step
  for (int i=0; i<nFcn; i++) {
    y[i] = y[i] + k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6;
  }
  return y;
}

// make the graphs that will return information on the solution
vector<TGraph> SetupGraphs(int nFcn){
  vector<TGraph> tg(nFcn);
  tg[0].SetTitle("y vs. x");
  tg[0].SetName("yx");
  tg[1].SetTitle("z vs. x");
  tg[1].SetName("zx");
  tg[2].SetTitle("v vs. t");
  tg[2].SetName("vt");
  return tg;
}

vector<double> RK4SolveNAProblem_1(vector<pfunc_t> &f, vector<double> &y0, int nsteps, double x0, double xmax, pfunc_t fstop, double errdef, int maxrep){
  // **
  // The following could be passed as parameters, to allow more control
  // when calling the function but here they are hardcoded to (hopefully)
  // resonable values to simplify the function interface
  const double hmin=errdef*4;
  const double hmax=(xmax-x0)/4;
  // **

  assert(f.size() == y0.size());   // need one initial condition per function

  vector<double> result;         // stores results for x_end and t_max

  int nFcn=f.size();
  double h=(xmax-x0)/nsteps;     // step size
  double h_last=0, h_new=0;      // for tracking adapted step size
  double x=x0;                   // independent variable
  vector<double> y=y0;           // start w/ init. conditions on dependent vars
  vector<double> y1,y2;          // temporary storage for adapting steps
  
  while (x<xmax){ // loop over steps in the x range
    int nreps=0;
    
    while (nreps<maxrep){ // take the steps, adapting step size as needed
      //full step
      y1=RK4StepN(f, y, x, h);
      // two half steps
      y2=RK4StepN(f, y, x, h/2);
      y2=RK4StepN(f, y2, x+h/2, h/2);
      h_last=h;           // last step size used to evaluate the functions
      
      // Calculate truncation error, comparing the full step to two half steps
      double errAbs=0, errRel=0;
      for (int i=0; i<nFcn; i++){
	      errAbs = max(errAbs,fabs(y1[i] - y2[i]));
	      errRel = max(errRel,fabs(errAbs / y2[i]));
      }
      double err=max(errAbs,errRel);  // error estimate
      err=max(err,1e-15);             // protect against unlikely case of 0 error
      h_new = h * pow (fabs (errdef / err), 0.2);
      h_new=max(h_new,hmin);          // restrict adapted step size to above limits
      h_new=min(h_new,hmax);
      if (err<errdef) {
	      h=h_new;
	      break;       // this step satisfies the requested errdef
      }
      nreps++;
      h=h_new;       // if not try again
    }

    if (fstop && fstop(x+h_last,y2)) break;
    if (abs(y2[0]-18.5)<=0.01 && abs(y2[2]-0.9)<=0.01) {
      result.push_back(y2[0]);
      result.push_back(x);
      break;
    }

    // advance to next step
    y=y2;
    x+=h_new;
  }
  return result;
}

vector<TGraph> RK4SolveNAProblem_2(vector<pfunc_t> &f, vector<double> &y0, int nsteps, double x0, double xmax, pfunc_t fstop, double errdef, int maxrep){
  // **
  // The following could be passed as parameters, to allow more control
  // when calling the function but here they are hardcoded to (hopefully)
  // resonable values to simplify the function interface
  const double hmin=errdef*4;
  const double hmax=(xmax-x0)/4;
  // **

  assert(f.size()==y0.size());   // need one initial condition per function

  int nFcn=f.size();
  double h=(xmax-x0)/nsteps;     // step size
  double h_last=0, h_new=0;      // for tracking adapted step size
  double x=x0;                   // independent variable
  vector<double> y=y0;           // start w/ init. conditions on dependent vars
  vector<double> y1,y2;          // temporary storage for adapting steps
  
  vector<TGraph> tg = SetupGraphs(3);
  tg[0].SetPoint(0,y0[0],y0[2]);                                        // tg[0] is the y vs x graph
  tg[1].SetPoint(0,y0[0],y0[4]);                                        // tg[1] is the z vs x graph
  tg[2].SetPoint(0,x0,sqrt(y0[1]*y0[1]+y0[3]*y0[3]+y0[5]*y0[5]));       // tg[2] is the v vs t graph

  while (x<xmax){ // loop over steps in the x range
    int nreps=0;
    
    while (nreps<maxrep){ // take the steps, adapting step size as needed
      //full step
      y1 = RK4StepN(f, y, x, h);
      // two half steps
      y2 = RK4StepN(f, y, x, h/2);
      y2 = RK4StepN(f, y2, x+h/2, h/2);
      h_last=h;           // last step size used to evaluate the functions
      
      // Calculate truncation error, comparing the full step to two half steps
      double errAbs=0, errRel=0;
      for (int i=0; i<nFcn; i++){
	      errAbs = max(errAbs,fabs(y1[i] - y2[i]));
	      errRel = max(errRel,fabs(errAbs / y2[i]));
      }
      double err=max(errAbs,errRel);  // error estimate
      err = max(err,1e-15);             // protect against unlikely case of 0 error
      h_new = h * pow (fabs (errdef / err), 0.2);
      h_new = max(h_new,hmin);          // restrict adapted step size to above limits
      h_new = min(h_new,hmax);
      if (err<errdef) {
	      h=h_new;
	      break;       // this step satisfies the requested errdef
      }
      nreps++;
      h=h_new;       // if not try again
    }
    
    if (fstop && fstop(x+h_last,y2)) break;  // fix me
    // advance to next step and store results in graphs
    y=y2;

    tg[0].SetPoint(tg[0].GetN(),y[0],y[2]);
    tg[1].SetPoint(tg[1].GetN(),y[0],y[4]);
    tg[2].SetPoint(tg[2].GetN(),x+h_last,sqrt(y[1]*y[1]+y[3]*y[3]+y[5]*y[5]));

    x+=h_new;
  }
  return tg;
}