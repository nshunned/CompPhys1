// Apply Selection cuts to S/B sample 
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
#include "cuts.h"  // contains the function to define events passing cuts
using namespace std;


void makeHists(TH1F **hs, TH1F **hb, TTree *s, TTree *b, TString var);


// fill histograms with events passing our selection cuts

void selection(double signalWgt=1.0, bool passAll=false) {
  gStyle->SetOptStat(false);
  // access the signal and backgound samples
  TString datFile="simulation.root";
  TFile *tf=new TFile(datFile);
  TTree *ts=(TTree*)tf->Get("TreeS");
  TTree *tb=(TTree*)tf->Get("TreeB");
  // for saving selected distros.
  TFile *tfo=new TFile(datFile.ReplaceAll(".root","_sel.root"),"recreate"); 

  // 1. make histograms for all five variables and 2D plots
  const int nvar=5;
  TH1F *hs[nvar], *hb[nvar];
  TString vars[nvar]={"var1","var2","var3","var4","var5"};
  for (int i=0; i<nvar; i++) makeHists(&hs[i],&hb[i],ts,tb,vars[i]);

  TH2F *hs2d[nvar*(nvar-1)/2];
  TH2F *hb2d[nvar*(nvar-1)/2];
  int count=0;
  for (int i=0; i<nvar; i++) {
    double y1=hs[i]->GetXaxis()->GetXmin();
    double y2=hs[i]->GetXaxis()->GetXmax();
    for (int j=i+1; j<nvar; j++) {
      double x1=hs[j]->GetXaxis()->GetXmin();
      double x2=hs[j]->GetXaxis()->GetXmax();
      TString name="S__"+vars[i]+"_v_"+vars[j]; 
      TString title=name+";"+vars[j]+";"+vars[i];
      hs2d[count]=new TH2F(name,title,30,x1,x2,30,y1,y2);
      name="B__"+vars[i]+"_v_"+vars[j];
      title=name+";"+vars[j]+";"+vars[i];
      hb2d[count]=new TH2F(name,title,30,x1,x2,30,y1,y2);
      count++;
    }
  }


  // 2. loop over trees and fill histograms w/ event passing cuts
  Float_t svar[nvar], bvar[nvar];
  TBranch *sb, *bb;
  for (int i=0; i<nvar; i++) {  // fetch the 5 variables into an array
    sb=ts->GetBranch(vars[i]);
    sb->SetAddress(&svar[i]);
    bb=tb->GetBranch(vars[i]);
    bb->SetAddress(&bvar[i]);
  }
  // fill signal histograms
  for (int n=0; n<ts->GetEntries(); n++) { // loop over tree and read values
    ts->GetEntry(n);
    if (passAll || pass(svar)) {
      for (int i=0; i<nvar; i++) {
        hs[i]->Fill(svar[i],signalWgt);
      }
      count=0;
      for (int i=0; i<nvar; i++) {
        for (int j=i+1; j<nvar; j++) { 
          hs2d[count]->Fill(svar[j],svar[i]); 
          count++;
        }
      }
    }
  }
  // fill bkg histograms
  for (int n=0; n<tb->GetEntries(); n++) {
    tb->GetEntry(n);
    if (passAll || pass(bvar)) {
      for (int i=0; i<nvar; i++) hb[i]->Fill(bvar[i]);
      count=0;
      for (int i=0; i<nvar; i++) {
        for (int j=i+1; j<nvar; j++) { 
          hb2d[count]->Fill(bvar[j],bvar[i]); 
          count++;
        }
      }
    }
  }

  // 3. plot results
  TCanvas *canSB=new TCanvas("cSB2","Signal and background",1450,300);
  canSB->Divide(nvar,1);
  for (int i=0; i<nvar; i++) {
    canSB->cd(i+1); 
    hs[i]->SetMinimum(0);
    hs[i]->SetMaximum(1.05*TMath::Max(hs[i]->GetMaximum(),hb[i]->GetMaximum()));
    hs[i]->Draw("hist"); 
    hb[i]->Draw("hist,same");
   }

  TCanvas *can2D=new TCanvas("c2D2","2D Plots",1450,500);
  can2D->Divide(5,2);
  count=0;
  for (int i=0; i<nvar; i++) {
    for (int j=i+1; j<nvar; j++) {
      can2D->cd(count+1);
      hs2d[count]->Draw();
      hb2d[count]->Draw("same,cont2");
      count++;
    }
  }
  

  // 4. calculate signal significance
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

  cout << "Saving histgrams in file: " << tfo->GetName() << endl;
  tfo->Write();
  //  tfo->Close();
}


// Find the min/max values of some variable is present in two trees
void minmax(double &min, double &max, TTree *s, TTree *b, TString var) {
  min=TMath::Min(s->GetMinimum(var),b->GetMinimum(var));
  max=TMath::Max(s->GetMaximum(var),b->GetMaximum(var));
}

// Create histograms with reasonable x-axis ranges
void makeHists(TH1F **hs, TH1F **hb, TTree *s, TTree *b, TString var) {
  double min,max;
  minmax(min,max,s,b,var);

  *hs = new TH1F("hs_"+var,"Signal: "+var,100,min,max);
  (*hs)->SetLineColor(kRed);
  (*hs)->Sumw2();

  *hb = new TH1F("hb_"+var,"Background: "+var,100,min,max);
  (*hb)->SetLineColor(kBlue);
  (*hb)->Sumw2();
}

