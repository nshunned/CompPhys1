#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"
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
// y[2] = position along y axis   ; f_rj = drj/dt => velocity along y axis  
// y[3] = velocity along y axis   ; f_vj = dvj/dt => acceleration along y axis
// y[4] = position along z axis   ; f_ri = dri/dt => velocity along z axis
// y[5] = velocity along z axis   ; f_vi = dvi/dt => acceleration along z axis

const double g=9.81;            // [m/s^2] gravity
const double B=4.1e-4;          // [] dimensionless quantity for a typical baseball

// curveball
const double v_0=85*0.44704;      // [m/s] initial velocity v_0 = 85 mph
const double theta=M_PI/180;      // [rad] angle of strike (above horizon) θ = 1°
const double omega=1800*M_PI/30;  // [rad/s] angular speed ω = 1800 rpm
const double phi=45*M_PI/180;     // [rad] baseball tilt φ = 45°

/*
// fastball
const double v_0=95*0.44704;      // [m/s] initial velocity v_0 = 95 mph
const double theta=M_PI/180;      // [rad] angle of strike (above horizon) θ = 1°
const double omega=1800*M_PI/30;  // [rad/s] angular speed ω = 1800 rpm
const double phi=225*M_PI/180;     // [rad] baseball tilt φ = 225°
*/

// time steps
const double h = 1e-4;            // time step width
const double t_min = 0;           // 0 < t < 5s
const double t_max = 1;

double f_ri(double x, const vector<double> &y){  // change in position along x axis
  (void) x;
  return y[1];
}
double f_vi(double x, const vector<double> &y){  // change in velocity along x axis
  (void) x;
  double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
  double f = 0.0039 + 0.0058/(1 + exp((v-35)/5));
  return -f*v*y[1] + B*omega*(y[5]*sin(phi) - y[3]*cos(phi));
}
double f_rj(double x, const vector<double> &y){  // change in position along y axis
  (void) x;
  return y[3];
}
double f_vj(double x, const vector<double> &y){  // change in velocity along y axis
  (void) x;
  double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
  double f = 0.0039 + 0.0058/(1 + exp((v-35)/5));
  return -f*v*y[3] + B*omega*y[1]*cos(phi);
}
double f_rk(double x, const vector<double> &y){  // change in position along z axis
  (void) x;
  return y[5];
}
double f_vk(double x, const vector<double> &y){  // change in velocity along z axis
  (void) x;
  double v = sqrt(y[1]*y[1] + y[3]*y[3] + y[5]*y[5]);
  double f = 0.0039 + 0.0058/(1 + exp((v-35)/5));
  return -g - f*v*y[5] - B*omega*y[1]*sin(phi);
}
double f_stop(double x, const vector<double> &y){
  (void) x;
  if (y[0]>18.5) return 1; // the plate is 60 ft 6 in (~ 18.5 m) away
  return 0;  // continue calculation
}

int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // 4 element vector of function pointers
  vector<pfunc_t> v_fun(6);
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;
  v_fun[4]=f_rk;
  v_fun[5]=f_vk;

  // initial conditions for position
  vector<double> y0(6);
  y0[0]=0;              // x(0)
  y0[2]=0;              // y(0)
  y0[4]=0;              // z(0)
  y0[1]=v_0*cos(theta); // vx(0)
  y0[3]=0;              // vy(0)
  y0[5]=v_0*sin(theta); // vz(0)
  
  auto result = RK4SolveNAProblem_2(v_fun, y0, (t_min-t_max)/h, t_min, t_max, f_stop);

  // plots the curveball model
  TCanvas* canvas = new TCanvas("canvas","Curve Ball",600,600);
  gStyle->SetOptStat("0");
  result[0].GetXaxis()->SetLimits(-1,20);
  result[0].GetYaxis()->SetRangeUser(-1.3,0.5);
  result[0].GetXaxis()->SetTitle("x (m)");
  result[0].GetYaxis()->SetTitle("y (m) / z (m)");
  result[0].GetXaxis()->SetTitleOffset(1.2);
  result[0].GetYaxis()->SetTitleOffset(1.2);
  result[0].SetLineStyle(7);
  result[1].SetLineStyle(1);

  result[0].Draw("AC");
  result[1].Draw("SAME");

  // plots the legend
  TLegend* legend = new TLegend(0.5, 0.875, 0.16, 0.75, "");
  legend->SetTextSize(0.05);
  legend->AddEntry((TGraph*)&result[0], "y vs. x");
  legend->AddEntry((TGraph*)&result[1], "z vs. x");
  legend->Draw();
  
  // prints graph
  canvas->Print("Curveball.png","png");
  // canvas->Print("Fastball.png","png");

  // writes plots to .root file
  TFile *tf=new TFile("Problem_2.root","recreate");
  for (unsigned i=0; i<result.size(); i++){
    result[i].Write();
  }
  tf->Close();

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}