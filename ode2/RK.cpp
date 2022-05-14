#include "RK.hpp"
#include <vector>

TGraph RK1Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax) {
  
  double h = (xmax-x0) / nsteps;     // step size
  double x = x0;                     // independent variable
  
  double *y = new double[nsteps];    // dependent variable to plot vs x
  TGraph tg;

  y[0] = y0;
  tg.SetPoint(0, x0, y0);
	      
  double k1;
  for (int i = 0; i < nsteps-1; i++){
    k1 = f(x, y[i]);
    y[i+1] = y[i] + h*k1;
    x += h;
    tg.SetPoint(i+1, x, y[i+1]);
  }
  delete[] y;
  return tg;
}


TGraph RK2Solve(double (*f)(double x, double y), double y0,
		 int nsteps, double x0, double xmax) {

  double h = (xmax-x0)/nsteps;     // step size
  double x = x0;                   // independent variable
  
  double *y = new double[nsteps];  // dependent variable to plot vs x
  TGraph tg;

  y[0] = y0;
  tg.SetPoint(0, x0, y0);
	      
  double ytmp, k1, k2;
  for (int i = 0; i < nsteps-1; i++){
    k1 = f(x, y[i]);
    ytmp = y[i] + h*k1 / 2;
    k2 = f(x + h/2, ytmp);
    y[i+1] = y[i] + h*k2;
    x += h;
    tg.SetPoint(i+1, x, y[i+1]);
  }
  delete[] y;
  return tg;
}

TGraph RK4Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax) {

  double h = (xmax-x0) / nsteps;     // step size                                                   
  double x = x0;                     // independent variable                                        

  double *y = new double[nsteps];    // dependent variable to plot vs x
  TGraph tg;

  y[0] = y0;
  tg.SetPoint(x, x0, y0);

  double ytmp1, ytmp2, ytmp3, k1, k2, k3, k4;
  for(int i = 0; i < nsteps-1; i++){
    k1 = f(x, y[i]);
    ytmp1 = y[i] + h*k1 / 2;
    k2 = f(x + h/2, ytmp1);
    ytmp2 = y[i] + h*k2 / 2;
    k3 = f(x + h/2, ytmp2);
    ytmp3 = y[i] + h*k3;
    k4 = f(x + h, ytmp3);
    
    y[i+1] = y[i] + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
    x += h;
    tg.SetPoint(i+1, x, y[i+1]);
  }
  delete[] y;
  return tg;
}
