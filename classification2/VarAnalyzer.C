// Example of a generic analysis tool
// This example is written to explore trees describing signal and background
// classes of events.  A real world example might automatically detect variable
// names and number of dimensions, but this is written for our limited example
// of 5 variables and the case of one signal model and one background model.

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

// Find the min/max values of some variable is present in two trees
void minmax(double &min, double &max, TTree *s, TTree *b, TString var){
  min=TMath::Min(s->GetMinimum(var),b->GetMinimum(var));
  max=TMath::Max(s->GetMaximum(var),b->GetMaximum(var));
}

// Fill signal and background histgrams for some variable
// An optional weight may be aplied to the signal
// Signal uses red line, background uses blue line
void filler(TH1F **hs, TH1F **hb, TTree *s, TTree *b, TString var, 
	    double sigWgt=1.0){
  double min,max;
  minmax(min,max,s,b,var);

  *hs = new TH1F("hs_"+var,"Signal: "+var,100,min,max);
  (*hs)->SetLineColor(kRed);

  *hb = new TH1F("hb_"+var,"Background: "+var,100,min,max);
  (*hb)->SetLineColor(kBlue);

  s->Draw(var+">>+"+(*hs)->GetName());
  b->Draw(var+">>+"+(*hb)->GetName());

  // signal scaling
  (*hs)->Scale(sigWgt);

  // make sure both historams are fully visible when plotted together
  double ymax=TMath::Max((*hs)->GetMaximum(),(*hb)->GetMaximum());
  (*hs)->SetMaximum(ymax*1.05);
  (*hb)->SetMaximum(ymax*1.05);
}

void filler2D(TH2F **hs, TH2F **hb, TTree *s, TTree *b, 
	      TString var1, TString var2);

double rho(TTree *t, TString var1, TString var2);

void cutScan(TH1F *s, TH1F *b, TH1F **cutL, TH1F **cutH, TString var);


// Perform an intitial scan of our variables
// Optionally apply a weight to reduce the signal fraction,
// thus affecting results for the signal significance in a cut scan

void VarAnalyzer(double signalWgt=1.0){
  gStyle->SetOptStat(false);
  // access the signal and backgound samples
  TString datFile="simulation.root";
  TFile *tf=new TFile(datFile);
  TTree *ts=(TTree*)tf->Get("TreeS");
  TTree *tb=(TTree*)tf->Get("TreeB");
  TFile *tfo=new TFile(datFile.ReplaceAll(".root","_vars.root"),"recreate"); 

  // 1. fill histograms for all five variables
  const int nvar=5;
  TH1F *hs[nvar], *hb[nvar];
  TString vars[nvar]={"var1","var2","var3","var4","var5"};
  for (int i=0; i<nvar; i++) filler(&hs[i],&hb[i],ts,tb,vars[i],signalWgt);

  // 2. find optimal cuts for singal significance
  TH1F *hsigL[nvar], *hsigH[nvar]; // histogram for low/high cuts on samples
  for (int i=0; i<nvar; i++) cutScan(hs[i],hb[i],&hsigL[i],&hsigH[i],vars[i]);

  // 3. fill histograms of 2D distributions
  TH2F *hs2d[nvar*(nvar-1)/2];
  TH2F *hb2d[nvar*(nvar-1)/2];
  int count=0;
  for (int i=0; i<nvar; i++){
    for (int j=i+1; j<nvar; j++){
      filler2D(&hs2d[count],&hb2d[count],ts,tb,vars[i],vars[j]);
      count++;
    }
  }
  
  // 4. calulate correlations
  TH2F *sRho=new TH2F("sRho","Signal correlations;Var A;Var B",
		      5,0.5,5.5,5,0.5,5.5);
  TH2F *bRho=new TH2F("bRho","Background correlations;Var A;Var B",
		      5,0.5,5.5,5,0.5,5.5);


  
  for (int i=0; i<nvar; i++){
    for (int j=0; j<nvar; j++){
      sRho->Fill(i+1,j+1,rho(ts,vars[i],vars[j]));
      bRho->Fill(i+1,j+1,rho(tb,vars[i],vars[j]));
    }
  }

  // 5. plot S/B distributions and result of cut scans
  TCanvas *canSB=new TCanvas("cSB","Signal and background",1450,500);
  canSB->Divide(nvar,2);
  for (int i=0; i<nvar; i++){
    canSB->cd(i+1); 
    hs[i]->Draw(); 
    hb[i]->Draw("same");
  }
  for (int i=0; i<nvar; i++){
    canSB->cd(i+6); 
    hsigL[i]->Draw(); 
    hsigH[i]->Draw("same");
  }
  canSB->Print("canSB.png");

  // 6. plot 2D results
  TCanvas *can2D=new TCanvas("c2D","2D Plots",1450,500);
  can2D->Divide(5,2);
  count=0;
  for (int i=0; i<nvar; i++){
    for (int j=i+1; j<nvar; j++){
      can2D->cd(count+1);
      hs2d[count]->Draw();
      hb2d[count]->Draw("same,cont2");
      count++;
    }
  }
  can2D->Print("can2D.png");
  
  // 7. plot correlation results
  TCanvas *canRho=new TCanvas("cRho","Var. Correlations",600,300);
  sRho->GetZaxis()->SetRangeUser(-1, 1);
  bRho->GetZaxis()->SetRangeUser(-1, 1);
  sRho->GetXaxis()->SetNdivisions(nvar);
  sRho->GetYaxis()->SetNdivisions(nvar);
  bRho->GetXaxis()->SetNdivisions(nvar);
  bRho->GetYaxis()->SetNdivisions(nvar);
 
  canRho->Divide(2,1);
  canRho->cd(1)->SetRightMargin(0.2); // to make labels more visible
  sRho->Draw("colz");
  canRho->cd(2)->SetRightMargin(0.2);
  bRho->Draw("colz");

  canRho->Print("canRho.png");

  cout << "Saving histgrams in file: " << tfo->GetName() << endl;
  tfo->Write(); // save histograms  
}

// Scan cut locations for events above low cut (solid line)
// and below high cut (dashed line)
// Return histograms showing signal significance for values of the cuts
void cutScan(TH1F *s, TH1F *b, TH1F **cutL, TH1F **cutH, TString var){
  *cutL=new TH1F("CL"+var,"C_L:"+var,s->GetNbinsX(),
		s->GetXaxis()->GetXmin(),s->GetXaxis()->GetXmax());
  *cutH=new TH1F("CH"+var,"C_H:"+var,s->GetNbinsX(),
		s->GetXaxis()->GetXmin(),s->GetXaxis()->GetXmax());
  for (int i=1; i<=s->GetNbinsX(); i++){ // scan possible cut values
    double ns=0, nb=0;
    ns=s->Integral(i+1,s->GetNbinsX());  // count events above cut
    nb=b->Integral(i+1,s->GetNbinsX());
    double sig=ns/TMath::Max(1.0,TMath::Sqrt(ns+nb));
    (*cutL)->SetBinContent(i,sig);
    ns=s->Integral()-ns;
    nb=b->Integral()-nb;
    sig=ns/TMath::Max(1.0,TMath::Sqrt(ns+nb));
    (*cutH)->SetBinContent(i,sig);
    (*cutH)->SetLineStyle(2);
  }
}


// make scatter plots of two variables
void filler2D(TH2F **hs, TH2F **hb, TTree *s, TTree *b, 
	      TString var1, TString var2){
  double xmin,xmax, ymin, ymax;
  minmax(xmin,xmax,s,b,var2);
  minmax(ymin,ymax,s,b,var1);

  TString name="S__"+var1+"_v_"+var2; 
  TString title=name+";"+var2+";"+var1;
  *hs = new TH2F(name,title,30,xmin,xmax,30,ymin,ymax);
  name="B__"+var2+"_v_"+var1;
  *hb = new TH2F(name,title,30,xmin,xmax,30,ymin,ymax);

  s->Draw(var1+":"+var2+">>+"+(*hs)->GetName());
  b->Draw(var1+":"+var2+">>+"+(*hb)->GetName());
}


double rho(TTree *t, TString var1, TString var2){
  if (var1==var2) return 1;
  // associate variables in the tree with local variables
  t->ResetBranchAddresses();  // clear any memory associations
  TBranch *b1=t->GetBranch(var1);
  TBranch *b2=t->GetBranch(var2);
  Float_t v1,v2;
  b1->SetAddress(&v1);
  b2->SetAddress(&v2);

  double sumx=0, sumy=0, sumxy=0, sumx2=0, sumy2=0;
  int n=t->GetEntries();
  for (int i=0; i<n; i++){ // loop over tree a read values
    t->GetEntry(i);
    sumx+=v1;
    sumy+=v2;
    sumxy+=v1*v2;
    sumx2+=v1*v1;
    sumy2+=v2*v2;
  }
  return (n*sumxy-sumx*sumy)/
    ( TMath::Sqrt((n-1)*sumx2-sumx*sumx) *
      TMath::Sqrt((n-1)*sumy2-sumy*sumy) );
}
