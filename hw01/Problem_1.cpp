#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "math.h"
#include <iostream>

using namespace std;

// first we define our vaiables
// x    = time
// y[0] = position along x axis   ; f_ri = dri/dt => velocity along x axis  
// y[1] = velocity along x axis   ; f_vi = dvi/dt => acceleration along x axis
// y[2] = position along z axis   ; f_rj = drj/dt => velocity along z axis  
// y[3] = velocity along z axis   ; f_vj = dvj/dt => acceleration along z axis

const double g=9.81;            // [m/s^2] gravity
const double m=0.145;           // [kg] mass
const double d=0.075;           // [m] diameter of baseball
const double b=1.6e-4*d;        // [m] air resistance constant
const double c=0.25*d*d;        // [m^2] air resistance constant
const double theta=M_PI/180;    // [rad] angle of strike (above horizon)

double f_ri(double x, const vector<double> &y){  // change in position along x axis
  (void) x;
  return y[1];
}
double f_vi(double x, const vector<double> &y){  // change in velocity along x axis
  (void) x;
  return -y[1] * (b + c * sqrt(y[1]*y[1] + y[3]*y[3])) / m ;
}
double f_rj(double x, const vector<double> &y){  // change in position along z axis
  (void) x;
  return y[3];
}
double f_vj(double x, const vector<double> &y){  // change in velocity along z axis
  (void) x;
  return -y[3] * (b + c * sqrt(y[1]*y[1] + y[3]*y[3])) / m - g;
}
double f_stop(double x, const vector<double> &y){
  (void) x;
  if (y[2]<0) return 1;
  return 0;  // continue calculation
}

int main(int argc, char **argv){
  // 4 element vector of function pointers
  vector<pfunc_t> v_fun(4);
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;

  // initial conditions for position
  vector<double> y0(4);
  y0[0]=0;     // x(0)
  y0[2]=1.4;   // z(0)

  double v_0 = 1;
  while(v_0 <= 100) {
    // initial conditions for velocity
    y0[1]=v_0*cos(theta); // vx(0)
    y0[3]=v_0*sin(theta); // vz(0)

    auto result = RK4SolveNAProblem_1(v_fun, y0, 200, 0, 20, f_stop);
    if(result.size() == 2) {
      cout << "v_0 = " << v_0 << " m/s = " << v_0*2.23694 << " mi/h" << endl;
      cout << "x_end = " << result[0] << " m" << endl;
      cout << "t_max = " << result[1] << " s" << endl;
      break;
    }

    v_0 += 0.1;
  }
}