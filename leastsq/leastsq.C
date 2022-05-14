#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"

using TMath::Log;

//parms
const double xmin=1;
const double xmax=20;
const int npoints=12;
const double sigma=0.2;

double f(double x){
  const double a=0.5;
  const double b=1.3;
  const double c=0.5;
  return a+b*Log(x)+c*Log(x)*Log(x);
}


void getX(double *x){
  double step=(xmax-xmin)/npoints;
  for (int i=0; i<npoints; i++){
    x[i]=xmin+i*step;
  }
}

void getY(const double *x, double *y, double *ey){
  static TRandom2 tr(0);
  for (int i=0; i<npoints; i++){
    y[i]=f(x[i])+tr.Gaus(0,sigma);
    ey[i]=sigma;
  }
}

double chiSquared(const double *y, const double *ey, const double *yFit){
  double chi2=0;
  for (int i=0; i<npoints; i++)
    chi2 += (y[i] - yFit[i])/ey[i] * (y[i] - yFit[i])/ey[i];
  return chi2;
}

void getYFit(const double *x, const double *y, const double *ey, double *fitParam, double *yFit){
  double Sy=0, S=0, Sx=0, Sxx=0, Syx=0, Sxxx=0, Syxx=0, Sxxxx=0;
  double si2=0, lxi=0;
  for (int i=0; i<npoints; i++){
    si2 = ey[i]*ey[i];
    lxi = Log(x[i]);

    Sy += y[i]/si2;
    S += 1/si2;
    Sx += lxi/si2;
    Sxx += lxi*lxi/si2;
    Syx += y[i]*lxi/si2;
    Sxxx += lxi*lxi*lxi/si2;
    Syxx += y[i]*lxi*lxi/si2;
    Sxxxx += lxi*lxi*lxi*lxi/si2;
  }

  double det = S*(Sxx*Sxxxx-Sxxx*Sxxx) - Sx*(Sx*Sxxxx-Sxx*Sxxx) + Sxx*(Sx*Sxxx - Sxx*Sxx);
  fitParam[0] = (Sy*(Sxx*Sxxxx-Sxxx*Sxxx) + Syx*(Sxx*Sxxx-Sx*Sxxxx) + Syxx*(Sx*Sxxx-Sxx*Sxx))/det;
  fitParam[1] = (Sy*(Sxx*Sxxx-Sx*Sxxxx) + Syx*(S*Sxxxx-Sxx*Sxx) + Syxx*(Sx*Sxx-S*Sxxx))/det;
  fitParam[2] = (Sy*(Sx*Sxxx-Sxx*Sxx) + Syx*(Sx*Sxx-S*Sxxx) + Syxx*(S*Sxx-Sx*Sx))/det;

  for (int i=0; i<npoints; i++)
    yFit[i] = fitParam[0] + fitParam[1]*Log(x[i]) + fitParam[2]*Log(x[i])*Log(x[i]);
}

void setAxes(TH1F* h) {
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);

  h->GetYaxis()->SetTitle("Count");
}

void leastsq(){
  // canvas setup
  TCanvas* canvas = new TCanvas("Four Plots", "Four Plots", 1400, 700);
  canvas->Divide(2,2);
  gStyle->SetOptStat("M");
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.25);
  gStyle->SetTitleSize(0.075,"t");
  gPad->SetLeftMargin(0.125);
  gPad->SetBottomMargin(0.125);

  // data set arrays
  double x[npoints];
  double y[npoints];
  double ey[npoints];

  // least square fit arrays
  double fitParam[3];
  double yFit[npoints];
  double chi2;

  // historgrams
  TH1F* hA = new TH1F("a", "Fit Parameter a",60,-0.5,1.5);
  TH1F* hB = new TH1F("b", "Fit Parameter b",50,0,3);
  TH1F* hC = new TH1F("c", "Fit Parameter c",50,0,1);
  TH1F* hChi2 = new TH1F("chi2", "Fit Chi Squared",60,0,30);

  // make the histograms prettier
  setAxes(hA);
  setAxes(hB);
  setAxes(hC);
  setAxes(hChi2);
  hA->GetXaxis()->SetTitle("a");
  hB->GetXaxis()->SetTitle("b");
  hC->GetXaxis()->SetTitle("c");
  hChi2->GetXaxis()->SetTitle("Chi2");

  // least square fit
  for (int i=0; i<5000; i++) {
    getX(x);
    getY(x,y,ey);
    getYFit(x,y,ey,fitParam,yFit);
    chi2 = chiSquared(y,ey,yFit);

    hA->Fill(fitParam[0]);
    hB->Fill(fitParam[1]);
    hC->Fill(fitParam[2]);
    hChi2->Fill(chi2);
  }

  // draw histograms
  canvas->cd(1);
  hA->Draw();
  canvas->cd(2);
  hB->Draw();
  canvas->cd(3);
  hC->Draw();
  canvas->cd(4);
  hChi2->Draw();

  //auto tg = new TGraphErrors(npoints,x,y,0,ey);
  //tg->Draw("alp");
}

