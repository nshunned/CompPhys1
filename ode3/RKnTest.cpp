#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

// a single 1st order ODE, written w/ old and new interfaces  
double fun(double x, double y){
  return (3*x*x+4*x+2)/2/(y-1);     // y0 = -1
}                                   // solution: 1-sqrt(x^3 + 2*x^2 +2*x+4)

// same function, but use interface that supports a coupled array of ODEs
double fun_vy(double x, const vector<double> &y){
  return (3*x*x+4*x+2)/2/(y[0]-1);   // y0[0] = -1
}                                    // solution: 1-sqrt(x^3 + 2*x^2 + 2*x+ 4)

 
// functions to describe simple projectile motion

int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  TCanvas *c1 = new TCanvas("c1","ODE solutions",dw,dh);
  // ******************************************************************************
  

  // *** test 1: Compare RK4Solve method with RK4SolveN
  // RK4SolveN should produce the same result for a single ODE

  // Solve using single equation solver
  TGraph tg4=RK4Solve(fun,-1,30,0,10);  // y0=-1
  tg4.SetMarkerStyle(kFullStar);
  tg4.SetMarkerSize(0.022*dh/8);
  tg4.SetMarkerColor(kBlue);
 
  TF1 fun_sol=TF1("fun_sol","1-sqrt(x*x*x+2*x*x+2*x+4)",0,10);   // exact solution
  fun_sol.SetLineColor(kBlack);
  fun_sol.SetLineWidth(3);
  fun_sol.SetLineStyle(2);

  // Solve using equation array solver
  vector<pfunc_t> v_fun(1);   // 1 element vector of function pointers
  v_fun[0]=fun_vy;
  vector<double> y0(1);
  y0[0]=-1;
  auto tgN = RK4SolveN(v_fun, y0, 30, 0, 10);
  tgN[0].SetMarkerStyle(kFullDiamond);
  tgN[0].SetMarkerSize(0.015*dh/8);
  tgN[0].SetMarkerColor(kRed);
  
  // plot the results
  tg4.SetTitle("solution to dy/dx=(3x^2+4x+2)/2/(y-1);x;y");
  tg4.Draw("Ap");
  fun_sol.Draw("same");
  tgN[0].Draw("p");
  
  TLegend *tl = new TLegend(0.5,0.7,0.9,0.9);
  tl->AddEntry(&tg4,"RK4 Solution","p");
  tl->AddEntry(&tgN[0],"RK4 Solution (array code)","p");
  tl->AddEntry(&fun_sol,"Exact Solution","l");
  tl->Draw();
  
  c1->Update();
  
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

