// A simple C++ program to illustrate the use of ROOT class TMinuit 
// for function minimization.  The example shows a Maximum Likelihood
// fit for the mean of an exponential pdf in which TMinuit 
// minimizes - 2*log(L).   The user must specify what to minimize in the 
// function fcn, shown in this file below.

// fcn passes back f = -2*ln L by reference; this is the function to minimize.
// The factor of -2 allows MINUIT to get the errors using the same
// recipe as for least squares, i.e., go up from the minimum by 1.

// TMinuit does not allow fcn to be a member function, and the function
// arguments are fixed, so one of the only ways to bring the data  
// into fcn is to declare a pointer to the data (xVecPtr) as global.

// For more info on TMinuit see root.cern.ch/root/html/TMinuit.html .

// Based on example by 
// Glen Cowan
// RHUL Physics
// 4 December 2006

// Update Bob Hirosky: Sep 2013


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <TMinuit.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TGraph.h>

#include "functions.h"

using namespace std;

// Declare pointer to data as global (not elegant, but TMinuit needs this).
TH1F *hdata;
// also passing the fit function, to make the Minuit fcn more generic
TF1 *fparam;


//-------------------------------------------------------------------------
// fcn passes back f = - 2*ln(L), the function to be minimized.

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  for (int i=0; i<npar; i++){
    fparam->SetParameter(i,par[i]);
  }


  f = calcCHI2(hdata,fparam);

 
}                         // end of fcn

//-------------------------------------------------------------------------

int main(int argc, char **argv) {

  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this can be ignored
  TApplication theApp("App", &argc, argv);

  TCanvas* canvas = new TCanvas();

  // create a histogem and fill it w/ randomly distributed data
  // based on our function
  double xmin = 0.0;
  double xmax = 200;

  
  TFile *f=new TFile("distros.root");
  TH1F *hexp=(TH1F*)f->Get("dist1");

  
  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 3;              // the number of parameters
 
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);              // the fcn to be minized 
  // (eg calculates chi^2 or NLL, not the
  // parameterization of the data!


  TF1* myfunc = new TF1("myfunc", gumbel2, xmin, xmax, npar);  // use this to fit the data and get results
  myfunc->SetParameters(19.9, 77.5, .152);


  TString parName[npar] = {"C","mu","a"};

  for (int i=0; i<npar; i++){
    minuit.DefineParameter(i, parName[i].Data(), 
			   myfunc->GetParameter(i), myfunc->GetParameter(i)/10, 0, 0);
  }


  // here we define the pointers to pass information to Minuit's fcn
  // not pretty, but works well
  hdata=hexp;
  fparam=myfunc;

  // Do the minimization!

  minuit.Migrad();       // Minuit's best minimization algorithm

  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
    minuit.GetParameter(i,outpar[i],err[i]);
  }
  // get dof
  int dof;
  double chi2=calcCHI2(hexp, myfunc, &dof);
  cout << "dof: " << dof-npar << endl;

  // do minos analysis
  Int_t ierflg = 0;

  gMinuit->mncomd("MINOS",ierflg);
  // or
  // gMinuit->mnmnos();

  // Plot the result.
  hdata->Draw("e");
  myfunc->SetParameters(outpar);
  myfunc->Draw("same");

  myfunc->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
  myfunc->SetLineColor(1);             //  black (default)
  myfunc->SetLineWidth(1);

  myfunc->GetXaxis()->SetTitle("x");
  myfunc->GetYaxis()->SetTitle("f(x;gaus2)");

  cout << "Prob: " << TMath::Prob(chi2,dof) << endl;

  TCanvas* c2 = new TCanvas();

  double x[101],y[101];
  double a=outpar[2];
  double da=err[2]*1.5/50;  // scan 1.5 sigma
  for (int i=-50;i<=50;i++){
    myfunc->SetParameter(2,a+da*i);
    x[i+50]=a+da*i;
    y[i+50]=calcCHI2(hexp, myfunc);
  }
  TGraph tg(100,x,y);
  tg.SetTitle("Scan of par[2];a;chi^2");
  tg.Draw("AC*");


  cout << "To exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program
  theApp.Run(true);
  canvas->Close();

  return 0;

}
