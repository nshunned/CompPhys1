// Apply Selection cuts to data sample 
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


void makeHists(TH1F **hd, TTree *d, TString var, TFile *tf);


// fill histograms with events passing our selection cuts
// we pass the MC file b/c the data should be binned to match the
// signal and background templates that we have already defined
void selectD(TString datFile, TString MCfile, bool passAll=false){
  gStyle->SetOptStat(false);
  // access the signal and backgound samples
  TFile *tfSim=new TFile(MCfile);
  TFile *tf=new TFile(datFile);
  TTree *td=(TTree*)tf->Get("TreeD");
  // for saving selected distros.
  TFile *tfo=new TFile(datFile.ReplaceAll(".root","_sel.root"),"recreate"); 

  // 1. make histograms for all five variables and 2D plots
  const int nvar=5;
  TH1F *hs[nvar];
  TString vars[nvar]={"var1","var2","var3","var4","var5"};
  for (int i=0; i<nvar; i++) makeHists(&hs[i],td,vars[i],tfSim);

  TH2F *hs2d[nvar*(nvar-1)/2];
  int count=0;
  for (int i=0; i<nvar; i++){
    double y1=hs[i]->GetXaxis()->GetXmin();
    double y2=hs[i]->GetXaxis()->GetXmax();
    for (int j=i+1; j<nvar; j++){
      double x1=hs[j]->GetXaxis()->GetXmin();
      double x2=hs[j]->GetXaxis()->GetXmax();
      TString name="D__"+vars[i]+"_v_"+vars[j]; 
      TString title=name+";"+vars[j]+";"+vars[i];
      hs2d[count]=new TH2F(name,title,30,x1,x2,30,y1,y2);
      count++;
    }
  }


  // 2. loop over trees and fill histograms w/ event passing cuts
  Float_t dvar[nvar];
  TBranch *db;
  for (int i=0; i<nvar; i++){  // fetch the 5 variables into an array
    db=td->GetBranch(vars[i]);
    db->SetAddress(&dvar[i]);
  }
 
  // fill data histograms
  for (int n=0; n<td->GetEntries(); n++){ // loop over tree and read values
    td->GetEntry(n);
    if (passAll || pass(dvar)) {
      for (int i=0; i<nvar; i++) hs[i]->Fill(dvar[i]);
      count=0;
      for (int i=0; i<nvar; i++){
        for (int j=i+1; j<nvar; j++){ 
          hs2d[count]->Fill(dvar[j],dvar[i]); 
          count++;
        }
      }
    }
  }

  // 3. plot results
  TCanvas *canD=new TCanvas("cD2","Data Distributions",1450,300);
  canD->Divide(nvar,1);
  for (int i=0; i<nvar; i++){
    canD->cd(i+1); 
    hs[i]->Draw(); 
   }

  TCanvas *can2DD=new TCanvas("c2DD2","2D Data Plots",1450,500);
  can2DD->Divide(5,2);
  count=0;
  for (int i=0; i<nvar; i++){
    for (int j=i+1; j<nvar; j++){
      can2DD->cd(count+1);
      hs2d[count]->Draw();
      count++;
    }
  }
  
  cout << "Saving histograms in file: " << tfo->GetName() << endl;
  tfo->Write();
  //  tfo->Close();
}




// Create histograms with reasonable x-axis ranges
void makeHists(TH1F **hd, TTree *d, TString var, TFile *tfSim){
  TH1F *hsim=(TH1F*)tfSim->Get("hb_"+var);
  int nbins=hsim->GetNbinsX();
  double min=hsim->GetBinLowEdge(1);
  double max=hsim->GetBinLowEdge(nbins)+hsim->GetBinWidth(nbins);
  *hd = new TH1F("hd_"+var,"Data: "+var,nbins,min,max);
  (*hd)->SetLineColor(kBlack);
  (*hd)->Sumw2();
}

