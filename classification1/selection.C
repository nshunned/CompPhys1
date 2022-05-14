// Apply Selection cus to S/B sample 
// Plot variable for surviving events
// Give signal significance for events passing selection

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

bool SB1=true;

// use this function to define events passing cuts
// a simple example of one cut is given
bool pass(Float_t *v){
  Float_t v0=v[0] , v1=v[1], v2=v[2], v3=v[3], v4=v[4];
  if (SB1) 
    return v0>64 && v1<118.8 && v2>0.2 && v3>72 && v4>-5;
  else 
    return v0>68 && v1>12 && v2>0.2 && v3>76 && v4>-4;
								 

  Float_t var5=v[4];
  if (var5>-8) return true;
  return false;
}



void makeHists(TH1F **hs, TH1F **hb, TTree *s, TTree *b, TString var);


// fill histograms with events passing our selection cuts

void selection(double signalWgt=1.0){
  if (signalWgt<0.9) SB1=false;
  gStyle->SetOptStat(false);
  // access the signal and backgound samples
  TFile *tf=new TFile("simulation.root");
  TTree *ts=(TTree*)tf->Get("TreeS");
  TTree *tb=(TTree*)tf->Get("TreeB");
  // for saving selected distros.
  TFile *tfo=new TFile("postSelection.root","recreate"); 

  // 1. make histograms for all five variables
  const int nvar=5;
  TH1F *hs[nvar], *hb[nvar];
  TString vars[nvar]={"var1","var2","var3","var4","var5"};
  for (int i=0; i<nvar; i++) makeHists(&hs[i],&hb[i],ts,tb,vars[i]);

  // 2. loop over trees and fill histograms w/ event passing cuts
  Float_t svar[nvar], bvar[nvar];
  TBranch *sb, *bb;
  for (int i=0; i<nvar; i++){  // fetch the 5 variables into an array
    sb=ts->GetBranch(vars[i]);
    sb->SetAddress(&svar[i]);
    bb=tb->GetBranch(vars[i]);
    bb->SetAddress(&bvar[i]);
  }
 
  // fill signal histograms
  for (int n=0; n<ts->GetEntries(); n++){ // loop over tree and read values
    ts->GetEntry(n);
    if (pass(svar)) {
      for (int i=0; i<nvar; i++) hs[i]->Fill(svar[i],signalWgt);
    }
  }
  // fill bkg histograms
  for (int n=0; n<tb->GetEntries(); n++){
    tb->GetEntry(n);
    if (pass(bvar)) {
      for (int i=0; i<nvar; i++) hb[i]->Fill(bvar[i]);
    }
  }


  // 3. calculate signal significance
  cout << "Signal Significance" << endl;
  double ns=ts->GetEntries()*signalWgt;
  double nb=tb->GetEntries();
  cout << "Before cuts: " << ns/TMath::Sqrt(ns+nb) << endl;
  double sPass=hs[0]->GetEntries()*signalWgt;
  double bPass=hb[0]->GetEntries();  
  cout << "After  cuts: " << sPass/TMath::Sqrt(sPass+bPass) << endl;
  cout << "Events in/out" << endl;
  cout << "Signal:\t\t" << ns << " " << sPass << endl;
  cout << "Background:\t" << nb << " " << bPass << endl;
  cout << "Results saved in: " << tf->GetName() << endl;
  tfo->Write();
  //  tfo->Close();

  // 4. plot results
  TCanvas *canSB=new TCanvas("cSB","Signal and background",1450,300);
  canSB->Divide(nvar,1);
  for (int i=0; i<nvar; i++){
    canSB->cd(i+1); 
    hs[i]->Scale(1/hs[i]->Integral());
    hb[i]->Scale(1/hb[i]->Integral());

    if (hb[i]->GetMaximum()>hs[i]->GetMaximum()){
      hb[i]->Draw();
      hs[i]->Draw("same"); 
    }
    else{
      hs[i]->Draw();
      hb[i]->Draw("same"); 
    }
   }

}


// Find the min/max values of some variable is present in two trees
void minmax(double &min, double &max, TTree *s, TTree *b, TString var){
  min=TMath::Min(s->GetMinimum(var),b->GetMinimum(var));
  max=TMath::Max(s->GetMaximum(var),b->GetMaximum(var));
}

// Create histograms with reasonable x-axis ranges
void makeHists(TH1F **hs, TH1F **hb, TTree *s, TTree *b, TString var){
  double min,max;
  minmax(min,max,s,b,var);

  *hs = new TH1F("hs_"+var,"Signal: "+var,100,min,max);
  (*hs)->SetLineColor(kRed);

  *hb = new TH1F("hb_"+var,"Background: "+var,100,min,max);
  (*hb)->SetLineColor(kBlue);
}
