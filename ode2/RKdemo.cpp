#include "RK.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <iostream>

using namespace std;

// the differential equation to be solved
double fun1(double x, double y) {
  (void)x;              // prevent unused variable warning
  return -2*y;          // f = y'(x,y) = -2 * y(x)  
}                       // solution: y(x) = 3 * exp(-2*x)

double fun2(double x, double y) {
  return -y/x-2/(x*x);  // prevent unused variable warning
                        // f = y'(x,y) = -y/x - 2/(x^2)
                        // solution: y(x) = -2*log(x)/x + 2/x
}

double fun3(double x, double y) {
  return y + sin(x);    // f = y'(x,y) = y + sin(x)
                        // solution: y(x) = 1.5*e^x - sin(x)/2 - cos(x)/2
}

int main(int argc, char **argv) {
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // solve our DEQ using RK1 or RK2 methods!
  // Two examples are given.  Choose a fucntion for testing
  // TGraph tg1 = RK1Solve(fun1, 3, 30, 0, 3);   // initial condition y(0)=3
  // TGraph tg2 = RK2Solve(fun1, 3, 30, 0, 3);
  // TGraph tg4 = RK4Solve(fun1, 3, 30, 0, 3);                                                    
  // TF1 fun_sol = TF1("fun_sol", "3*exp(-2*x)", 0, 3);   // exact solution

  // TGraph tg1 = RK1Solve(fun2, 2, 100, 1, 100);   // initial condition y(1)=2
  // TGraph tg2 = RK2Solve(fun2, 2, 100, 1, 100);
  // TGraph tg4 = RK4Solve(fun2, 2, 100, 1, 100);
  // TF1 fun_sol = TF1("fun_sol", "-2*log(x)/x + 2/x", 1, 100);   // exact solution

  TGraph tg1 = RK1Solve(fun3, 1, 100, 0, 10);   // initial condition y(0)=1
  TGraph tg2 = RK2Solve(fun3, 1, 100, 0, 10);
  TGraph tg4 = RK4Solve(fun3, 1, 100, 0, 10);
  TF1 fun_sol = TF1("fun_sol", "1.5*exp(x) - sin(x)/2 - cos(x)/2", 0, 50);   // exact solution

   // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  tg1.SetMarkerSize(0.03*dh/8);
  tg2.SetMarkerSize(0.03*dh/8);
  tg4.SetMarkerSize(0.03*dh/8);
  TCanvas *c1 = new TCanvas("c1","DEQ solutions", dw, dh);
  // ******************************************************************************

  tg1.SetMarkerColor(kRed);
  tg2.SetMarkerColor(kGreen-2);
  tg4.SetMarkerColor(kBlue);
  fun_sol.SetLineColor(kBlack);
  fun_sol.SetLineStyle(2);

  tg1.GetXaxis()->SetTitle("x");
  tg1.GetYaxis()->SetTitle("y");
  tg1.SetTitle("1.5*e^x - sin(x)/2 - cos(x)/2");
 
  // plot the results
  tg1.Draw("A*");
  tg2.Draw("*");
  tg4.Draw("*");
  fun_sol.Draw("SAME");
  
  TLegend *tl = new TLegend(0.6, 0.7, 0.9, 0.9);
  tl->AddEntry(&tg1, "RK1 Solution", "p");
  tl->AddEntry(&tg2, "RK2 Solution", "p");
  tl->AddEntry(&tg4, "RK4 Solution", "p");
  tl->AddEntry(&fun_sol, "Exact Solution", "l");
  tl->Draw();
  c1->Draw();
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30, ".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
