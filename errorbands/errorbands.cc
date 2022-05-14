// Perform a fit and calculate the error bands based on 
// the error matrix
// usage:
// .L errorbands.cc+
// errorbands()

#include "TH1F.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGClient.h"

using namespace ROOT::Math;


Double_t fcn(Double_t *x, Double_t *par){
  //return par[0] + par[1]*TMath::Log(x[0]) + par[2]*TMath::Log(x[0])*TMath::Log(x[0]);
  return par[0] + par[1]*TMath::Sin(par[2]*x[0]);
}

// calculate 1 sigma uncertainy of fcn model wrt parameter uncertainties
double sigmaFit(Double_t *grad, TMatrixD &COV, Int_t npar){
  double sigma=0;

  // s^2 = sum_i sum_j df/dp_i df/dp_j COV_ij
  for (int pi=0; pi<npar; pi++){
    for (int pj=0; pj<npar; pj++){
      sigma+=grad[pi]*grad[pj]*COV[pi][pj];
    }
  }
  //return sigma;
  return TMath::Sqrt(sigma);
}

// alternate example code for calculating 1 sigma error
double sigmaFit2(Double_t *grad, TMatrixD &COV, Int_t npar){
  double var;
  TVectorD g(npar,grad);
  var = g*(COV*g);
  return TMath::Sqrt(var);
}

// here we can do error propagation for an arbitrary function
// return sigma_f(x;par)
double OneSD(WrappedTF1 &f, Double_t x, Double_t *pars, TMatrixD &COV, Int_t npar){
  Double_t *grad=new Double_t[npar];
  f.ParameterGradient(x,pars,grad);
  TVectorD g(npar,grad);
  double var = g*(COV*g);
  delete[] grad;
  return TMath::Sqrt(var);
}


void errorbands(){
  UInt_t dh = gClient->GetDisplayHeight();

  TCanvas *c1=new TCanvas("c1","fit",50,60,dh/2,dh/3);

  TFile *tf=new TFile("mydata.root");
  TH1F *data=(TH1F*)tf->Get("data");
  const Int_t nPars=3;
  Double_t pars[nPars], grad[nPars];
  TF1 model("model",fcn,1,500,nPars);
  model.SetParameter(0,1.23);
  model.SetParameter(1,.5);
  model.SetParameter(2,.2);
  data->Fit(&model);
  data->Draw();

  model.GetParameters(pars);
  WrappedTF1 wmodel(model);

  // Covariance matrix
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  //  int nPars = fitter->GetNumberFreeParameters();
  TMatrixD* COV = new TMatrixD( nPars, nPars, fitter->GetCovarianceMatrix() );  
  TH1F *up=new TH1F(*data);  up->SetLineColor(kBlue); 
  TH1F *dn=new TH1F(*data);  dn->SetLineColor(kBlue);
  for (int ib=1; ib<=data->GetNbinsX(); ib++){
    double x=data->GetBinCenter(ib);
    //wmodel.ParameterGradient(x,pars,grad);
    //double sigma=sigmaFit(grad, *COV, nPars);
    //double sigma=sigmaFit2(grad, *COV, nPars);
    double sigma=OneSD(wmodel,x,pars,*COV,nPars);
    up->SetBinContent(ib,model(x)+sigma);
    dn->SetBinContent(ib,model(x)-sigma);
  }
  up->Draw("same,hist,c");
  dn->Draw("same,hist,c");
  c1->Print("fit.png");
  //return;
  TCanvas *c2=new TCanvas("c2","bands",dh/2+80,60,dh/2,dh/3);

  TH1F *uprat=new TH1F(*up); uprat->Divide(&model);
  TH1F *dnrat=new TH1F(*dn); dnrat->Divide(&model);
  uprat->SetLineWidth(2);  
  dnrat->SetLineWidth(2);  

  // data/fit
  TH1F *data2=(TH1F*)data->Clone("ratio");
  data2->SetTitle("Ratio;x;data/fit");
  data2->Divide(&model);
  data2->SetMaximum(data2->GetMaximum()*1.1);
  data2->Draw();
  
  uprat->Draw("hist,c,same");
  dnrat->Draw("hist,c,same");

  c2->Print("bands.png");

}


