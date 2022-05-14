#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

using namespace std;


double scan(TH1F *s, TH1F *b){
  int nbins=s->GetNbinsX();
  double sig=0;
  int bestcut=0;
  int lowcut=1; // lowside cut
  for (int i=1; i<nbins; i++){
    double lo_s=s->Integral(1,i);
    double hi_s=s->Integral(i+1,nbins);
    double lo_b=b->Integral(1,i);
    double hi_b=b->Integral(i+1,nbins);
    double sig1=lo_s/sqrt(lo_s+lo_b);
    double sig2=hi_s/sqrt(hi_s+hi_b);
    if (sig1>sig) {sig=sig1; bestcut=i; lowcut=0;}  // try high side cut
    if (sig2>sig) {sig=sig2; bestcut=i; lowcut=1;}  // try low side
    //   cout << i << " " << sig << " " << sig1 << " " << sig2 << endl;
  }
  if (lowcut) cout << "Max significance cut for low  cut on ";
  else cout << "Max significance cut for high cut on ";
	     
  cout << s->GetName() << " = " 
       << sig << "  at bin " << bestcut << " : " << s->GetBinLowEdge(bestcut+1) << endl;
  return sig;
}


void scanner(float sigwgt=1.0){
  TFile *in=new TFile("simulation.root");
  TTree *s=(TTree*)in->Get("TreeS");
  TTree *b=(TTree*)in->Get("TreeB");

  s->SetWeight(sigwgt);
  // make S,B histograms
  s->Draw("var1>>hs1(100,0,400)");
  s->Draw("var2>>hs2(100,0,120)");
  s->Draw("var3>>hs3(100,0,1)");
  s->Draw("var4>>hs4(100,0,400)");
  s->Draw("var5>>hs5(100,-30,20)");
  b->Draw("var1>>hb1(100,0,400)");
  b->Draw("var2>>hb2(100,0,120)");
  b->Draw("var3>>hb3(100,0,1)");
  b->Draw("var4>>hb4(100,0,400)");
  b->Draw("var5>>hb5(100,-30,20)");

  double ns=( (TH1F*)gDirectory->Get("hs1") )->Integral();
  double nb=( (TH1F*)gDirectory->Get("hb1") )->Integral();
  cout << "Initial S/sqrt(S+B): " << ns/sqrt(ns+nb) << endl;
  
  // scan histos
  scan((TH1F*)gDirectory->Get("hs1"),(TH1F*)gDirectory->Get("hb1"));
  scan((TH1F*)gDirectory->Get("hs2"),(TH1F*)gDirectory->Get("hb2"));
  scan((TH1F*)gDirectory->Get("hs3"),(TH1F*)gDirectory->Get("hb3"));
  scan((TH1F*)gDirectory->Get("hs4"),(TH1F*)gDirectory->Get("hb4"));
  scan((TH1F*)gDirectory->Get("hs5"),(TH1F*)gDirectory->Get("hb5"));
}
