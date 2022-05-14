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
#include <TRandom2.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>

using namespace std;

// declare pointer to data as global (not elegant but TMinuit needs this).
TH1F *hdata;
// also passing the fit function, to make the Minuit fcn more generic
TF1 *fparam;

//-------------------------------------------------------------------------
// the pdf to be fitted, here an exponential
// this is a convenient interface to associating the model with a TF1 later
double expPdf(double* xPtr, double par[]) {        
  double x = *xPtr;
  double A = par[0];         // normalization
  double lam = par[1];       // mean of x
  return A*(1.0/lam) * exp(x/lam);
}

double gumPdf(double* xPtr, double par[]) {
  double x = *xPtr;
  double a = par[0];
  double b = par[1];
  x = x-par[2];
  return par[3]*a*b * exp(-(b * exp(-a*x) + a*x));
}

double gausPdf(double* xPtr, double par[]) {
  double x = *xPtr;
  double a1 = par[0];
  double b1 = par[1];
  double c1 = par[2];
  double a2 = par[3];
  double b2 = par[4];
  double c2 = par[5];
  return a1*exp(-(x-b1)*(x-b1)/(2*c1*c1)) + a2*exp(-(x-b2)*(x-b2)/(2*c2*c2));
}

//-------------------------------------------------------------------------
// return NLL given a histogram and function
double calcNLL(TH1F* h, TF1* f) {
  double nll=0;
  for (int i=1; i<=h->GetNbinsX(); i++) {
    double x = h->GetBinCenter(i);
    int n = (int)(h->GetBinContent(i));
    double mu = f->Eval(x);
    if (mu<1e-10) mu = 1e-10;    // avoid log(0) problems!
    nll -= n * TMath::Log(mu) - mu  - TMath::LnGamma(n+1);
  }
  // cout << "nll " << nll << endl;
  return 2*nll;   // factor of -2 so minuit gets the errors right
}

double calcChi2(TH1F* h, TF1* f){
  double chi2=0;
  for (int i=1; i<=h->GetNbinsX(); i++) {
    int y = (int)(h->GetBinContent(i));
    double ey = h->GetBinError(i);
    if (ey == 0)
      continue;
    double yFit = f->Eval(h->GetBinCenter(i));
    chi2 += (y-yFit)/ey * (y-yFit)/ey;
  }
  return chi2;
}

//-------------------------------------------------------------------------
// Minuit fcn: calculates value of the function to be minimized.
// this is the interface for our objective function
void fcn(int& npar, double* deriv, double& f, double par[], int flag) {
  for (int i=0; i<npar; i++) {
    fparam->SetParameter(i,par[i]);
  }
  f = calcNLL(hdata,fparam);
}

void chi2Fcn(int& npar, double* deriv, double& f, double par[], int flag) {
for (int i=0; i<npar; i++) {
    fparam->SetParameter(i,par[i]);
  }
  f = calcChi2(hdata,fparam);
}

//-------------------------------------------------------------------------
// styling canvas
void setCanvas(TCanvas* canvas) {
  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);        
  canvas->SetFrameBorderMode(0);  
  gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();
  gROOT->ForceStyle();

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();
}

//-------------------------------------------------------------------------
// the example expFit
void example() {
  // setup canvas
  TCanvas* canvas = new TCanvas();
  setCanvas(canvas);

  // generate some data in a histogram
  TRandom2 r(0);

  // create a histogram and fill it w/ randomly distributed data based on our function
  double xmin = 0.0;
  double xmax = 5.0;
  TH1F *hexp = new TH1F("hexp","exponential distribution",100,xmin,xmax);
  for (int i=0; i<1000; i++) {
    hexp->Fill(r.Exp(3.14159));
  }
  
  // initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);              // the fcn to be minized 
  // (eg calculates chi^2 or NLL, not the parameterization of the data!

  TF1* myfunc = new TF1("myfunc", expPdf, xmin, xmax, npar);  // use this to fit the data and get results
  double par[npar];               // the start values
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  TString parName[npar];

  par[0] = hexp->GetMaximum();       // guesses for starting the fit
  par[1] = -2;                       // this MUST be done by some means to get things started
  stepSize[0] = TMath::Abs(par[0]*0.1);   // usually 10% is OK for an initial step size, YMMV
  stepSize[1] = TMath::Abs(par[1]*0.1);   // step size MUST be positive!
  minVal[0] = 0;      // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 0;
  minVal[1] = -100; 
  maxVal[1] = -.1;

  parName[0] = "A";       // let's give our fit parameters useful names
  parName[1] = "lambda";

  // initialize the parameters
  for (int i=0; i<npar; i++) {
    minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], minVal[i], maxVal[i]);
  }

  // here we define the pointers to pass information to Minuit's fcn
  // not pretty, but simple
  hdata=hexp;
  fparam=myfunc;

  // Do the minimization!
  minuit.Migrad();       // Minuit's best minimization algorithm
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++) {
    minuit.GetParameter(i,outpar[i],err[i]);
  }

  // Plot the result.
  hdata->Draw("e");
  myfunc->SetParameters(outpar);
  myfunc->Draw("same");

  myfunc->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
  myfunc->SetLineColor(1);             //  black (default)
  myfunc->SetLineWidth(1);

  myfunc->GetXaxis()->SetTitle("x");
  myfunc->GetYaxis()->SetTitle("f(x;#lambda)");
}

//-------------------------------------------------------------------------
// homework problem
void homework() {
  // retrieve historgram to be fitted
  TFile* file = new TFile("distros.root");
  TH1F* hist = new TH1F(*(TH1F*)file->Get("dist1"));
  double xmin = hist->GetBinCenter(hist->FindFirstBinAbove()-1);
  double xmax = hist->GetBinCenter(hist->FindLastBinAbove()+1);
  
  const int npar1 = 6;
  const int npar2 = 4;

  TMinuit minuit1(npar1);
  TMinuit minuit2(npar2);
  minuit1.SetFCN(chi2Fcn);
  minuit2.SetFCN(chi2Fcn);

  // functions to fit the data and get results
  TF1* func1 = new TF1("func1", gausPdf, xmin, xmax, npar1);
  TF1* func2 = new TF1("func2", gumPdf, xmin, xmax, npar2);

  // set up arrays for minimizer 1
  TString parName1[npar1];
  double par1[npar1];
  double stepSize1[npar1];
  double minVal1[npar1];
  double maxVal1[npar1];

  parName1[0] = "A1"; parName1[1] = "mu1"; parName1[2] = "sigma1";
  parName1[3] = "A2"; parName1[4] = "mu2"; parName1[5] = "sigma2";
  par1[0] = 1000; par1[1] = 80; par1[2] = 10;
  par1[3] = -4000; par1[4] = 70; par1[5] = 5;
  for (int i=0; i<npar1; i++) {
    stepSize1[i] = TMath::Abs(par1[i]*0.01);
    minVal1[i] = 0;
    maxVal1[i] = 0;
  }

  // set up arrays for minimizer 2
  TString parName2[npar2];
  double par2[npar2];
  double stepSize2[npar2];
  double minVal2[npar2];
  double maxVal2[npar2];

  parName2[0] = "A";
  parName2[1] = "B";
  parName2[2] = "x0";
  parName2[3] = "C";
  par2[0] = sqrt(hist->GetMaximum());
  par2[1] = sqrt(hist->GetMaximum());
  par2[2] = 80;
  par2[3] = 100;
  for (int i=0; i<npar2; i++) {
    stepSize2[i] = TMath::Abs(par2[i]*0.01);
    minVal2[i] = 0;
    maxVal2[i] = 0;
  }

  // initializeParameters
  for (int i=0; i<npar1; i++) {
    minuit1.DefineParameter(i, parName1[i].Data(), par1[i], stepSize1[i], minVal1[i], maxVal1[i]);
  }
  for (int i=0; i<npar2; i++) {
    minuit2.DefineParameter(i, parName2[i].Data(), par2[i], stepSize2[i], minVal2[i], maxVal2[i]);
  }

  hdata=hist;

  // minimize using func1
  fparam=func1;
  minuit1.Migrad();
  double outpar1[npar1], err1[npar1];
  for (int i=0; i<npar1; i++) {
    minuit1.GetParameter(i,outpar1[i],err1[i]);
  }

  // minimize using func2
  fparam=func2;
  minuit2.Migrad();
  double outpar2[npar2], err2[npar2];
  for (int i=0; i<npar2; i++) {
    minuit2.GetParameter(i,outpar2[i],err2[i]);
  }

  // setup canvas
  TCanvas* canvas = new TCanvas();
  setCanvas(canvas);
  canvas->Divide(2,1);

  // fitted using sum of two Gaussian
  canvas->cd(1);
  hdata->Draw("e");
  func1->SetParameters(outpar1);
  func1->Draw("same");

  func1->SetLineStyle(1); // 1 = solid, 2 = dashed, 3 = dotted
  func1->SetLineColor(1); // black (default)
  func1->SetLineWidth(1);
  func1->GetXaxis()->SetTitle("x");
  func1->GetYaxis()->SetTitle("f(x;#mu_{1},#sigma_{1},#mu_{2},#sigma_{2})");

  // fitted using Gumbel
  canvas->cd(2);
  hdata->Draw("e");
  func2->SetParameters(outpar2);
  func2->Draw("same");

  func2->SetLineStyle(1); // 1 = solid, 2 = dashed, 3 = dotted
  func2->SetLineColor(1); // black (default)
  func2->SetLineWidth(1);
  func2->GetXaxis()->SetTitle("x");
  func2->GetYaxis()->SetTitle("f(x;a,b)");
}

//-------------------------------------------------------------------------
int main(int argc, char **argv) {
  // this allows you to view ROOT-based graphics in your C++ program
  // if you don't want view graphics (eg just read/process/write data files), this line can be removed
  TApplication theApp("App", &argc, argv);

  //example();
  homework();


  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program
  theApp.Run(true);
  //canvas->Close();

  return 0;
}
