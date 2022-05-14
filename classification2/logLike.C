using namespace std;

void makeHists(TH1F **hs, TH1F **hb, TTree *s, TTree *b, TString var);
void makeHists2(TH2F **hs2, TH2F **hb2, TTree *s, TTree *b, TString varX, TString varY);

void logLike(TString datFile, double signalWgt=1.0) {
  gStyle->SetOptStat(false);
  // access the signal and backgound samples
  TFile *fSim=new TFile("simulation.root");
  TTree *ts=(TTree*)fSim->Get("TreeS");
  TTree *tb=(TTree*)fSim->Get("TreeB");
  // access the data
  TFile *fDat=new TFile(datFile);
  TTree *td=(TTree*)fDat->Get("TreeD");
  // for saving selected distros.
  TFile *fOut=new TFile("logLike.root","recreate");

  // 1. make histograms for all five variables and 2D plots
  const int nvar=5;
  TH1F *hs[nvar], *hb[nvar];
  TH2F *hs2, *hb2;
  TString vars[nvar]={"var1","var2","var3","var4","var5"};
  for (int i=0; i<nvar; i++) {
    makeHists(&hs[i],&hb[i],ts,tb,vars[i]);
  }
  int X=1, Y=2;
  makeHists2(&hs2,&hb2,ts,tb,TString("var")+X,TString("var")+Y);

  // 2. loop over trees and fill histograms w/ event passing cuts
  Float_t svar[nvar], bvar[nvar], dvar[nvar];
  TBranch *sb, *bb, *db;
  for (int i=0; i<nvar; i++) {  // fetch the 5 variables into an array
    sb=ts->GetBranch(vars[i]);
    sb->SetAddress(&svar[i]);
    bb=tb->GetBranch(vars[i]);
    bb->SetAddress(&bvar[i]);
    db=td->GetBranch(vars[i]);
    db->SetAddress(&dvar[i]);
  }
  // fill signal histograms
  for (int n=0; n<ts->GetEntries(); n++) { // loop over tree and read values
    ts->GetEntry(n);
    for (int i=0; i<nvar; i++) {
      hs[i]->Fill(svar[i],signalWgt);
    }
    hs2->Fill(svar[0],svar[1],signalWgt);
  }
  // fill bkg histograms
  for (int n=0; n<tb->GetEntries(); n++) {
    tb->GetEntry(n);
    for (int i=0; i<nvar; i++) {
      hb[i]->Fill(bvar[i]);
    }
    hb2->Fill(bvar[0],bvar[1]);
  }

  // 3. normalize signal and background
  for (int i=0; i<nvar; i++) {
    hs[i]->Scale(1./hs[i]->Integral());
    hb[i]->Scale(1./hb[i]->Integral());
  }
  hs2->Scale(1./hs2->Integral());
  hb2->Scale(1./hs2->Integral());

  // 4. make histograms for LLs and LLR
  TH1F *hlls[3], *hllr[3], *hll2[3];
  TString lls[3]={"sll_s","bll_s","dll_s"};
  TString llr[3]={"sll_r","bll_r","dll_r"};
  TString ll2[3]={"sll_2","bll_2","dll_2"};
  for (int i=0; i<3; i++) {
    hlls[i]=new TH1F(lls[i],"Singal LL",440,-210,10);
    hllr[i]=new TH1F(llr[i],"S to B LLR",440,-210,10);
    hll2[i]=new TH1F(ll2[i],"S to B LLR2",440,-210,10);
    hlls[i]->SetLineWidth(2);
    hllr[i]->SetLineWidth(2);
    hll2[i]->SetLineWidth(2);
  }
  hlls[0]->SetLineColor(kRed); hlls[1]->SetLineColor(kBlue); hlls[2]->SetLineColor(kBlack);
  hllr[0]->SetLineColor(kRed); hllr[1]->SetLineColor(kBlue); hllr[2]->SetLineColor(kBlack);
  hll2[0]->SetLineColor(kRed); hll2[1]->SetLineColor(kBlue); hll2[2]->SetLineColor(kBlack);

  // 5. fill histograms for LLs and LLR
  // signal events compared to signal and background
  for (int n=0; n<ts->GetEntries(); n++) {
    ts->GetEntry(n);
    // LL
    double LLs=0, LL2=0;
    for (int i=0; i<nvar; i++) {
      double probSS=hs[i]->GetBinContent(hs[i]->FindBin(svar[i]));
      LLs+=(probSS==0 ? (-50) : log(probSS));
      if ((i!=X) && (i!=Y))
        LL2+=(probSS==0 ? (-50) : log(probSS));
    }
    double probSS=hs2->GetBinContent(hs2->FindBin(svar[X],svar[Y]));
    LL2+=(probSS==0 ? (-50) : log(probSS));
    // LLR
    double LLr=LLs;
    for (int i=0; i<nvar; i++) {
      double probSB=hb[i]->GetBinContent(hb[i]->FindBin(svar[i]));
      LLr-=(probSB==0 ? (-50) : log(probSB));
      if ((i!=X) && (i!=Y))
        LL2-=(probSB==0 ? (-50) : log(probSB));
    }
    double probSB=hb2->GetBinContent(hb2->FindBin(svar[X],svar[Y]));
    LL2-=(probSB==0 ? (-50) : log(probSB));
    // fill
    hlls[0]->Fill(LLs);
    hllr[0]->Fill(LLr);
    hll2[0]->Fill(LL2);
  }
  // background events compared to signal and background
  for (int n=0; n<tb->GetEntries(); n++) {
    tb->GetEntry(n);
    // LL
    double LLb=0, LL2=0;
    for (int i=0; i<nvar; i++) {
      double probBS=hs[i]->GetBinContent(hs[i]->FindBin(bvar[i]));
      LLb+=(probBS==0 ? (-50) : log(probBS));
      if ((i!=X) && (i!=Y))
        LL2+=(probBS==0 ? (-50) : log(probBS));
    }
    double probBS=hs2->GetBinContent(hs2->FindBin(bvar[X],bvar[Y]));
    LL2+=(probBS==0 ? (-50) : log(probBS));
    // LLR
    double LLr=LLb;
    for (int i=0; i<nvar; i++) {
      double probBB=hb[i]->GetBinContent(hb[i]->FindBin(bvar[i]));
      LLr-=(probBB==0 ? (-50) : log(probBB));
      if ((i!=X) && (i!=Y))
        LL2-=(probBB==0 ? (-50) : log(probBB));
    }
    double probBB=hb2->GetBinContent(hb2->FindBin(bvar[X],bvar[Y]));
    LL2-=(probBB==0 ? (-50) : log(probBB));
    // fill
    hlls[1]->Fill(LLb);
    hllr[1]->Fill(LLr);
    hll2[1]->Fill(LL2);
  }
  // data events compared to signal and background
  for (int n=0; n<td->GetEntries(); n++) {
    td->GetEntry(n);
    // LL
    double LLd=0, LL2=0;
    for (int i=0; i<nvar; i++) {
      double probDS=hs[i]->GetBinContent(hs[i]->FindBin(dvar[i]));
      LLd+=(probDS==0 ? (-50) : log(probDS));
      if ((i!=X) && (i!=Y))
        LL2+=(probDS==0 ? (-50) : log(probDS));
    }
    double probDS=hs2->GetBinContent(hs2->FindBin(dvar[X],dvar[Y]));
    LL2+=(probDS==0 ? (-50) : log(probDS));
    // LLR
    double LLr=LLd;
    for (int i=0; i<nvar; i++) {
      double probDB=hb[i]->GetBinContent(hb[i]->FindBin(dvar[i]));
      LLr-=(probDB==0 ? (-50) : log(probDB));
      if ((i!=X) && (i!=Y))
        LL2-=(probDB==0 ? (-50) : log(probDB));
    }
    double probDB=hb2->GetBinContent(hb2->FindBin(dvar[X],dvar[Y]));
    LL2-=(probDB==0 ? (-50) : log(probDB));
    // fill
    hlls[2]->Fill(LLd);
    hllr[2]->Fill(LLr);
    hll2[2]->Fill(LL2);
  }

  // 6. plot results
  TCanvas *canLLs=new TCanvas("canLLs","Signal Log Likelihood",1000,600);
  hlls[0]->Draw("hist");
  hlls[1]->Draw("hist,same");
  hlls[2]->Draw("hist,same");
  TLegend* leg1=new TLegend(0.125,0.7,0.425,0.9);
  leg1->AddEntry(hlls[0],"Signal");
  leg1->AddEntry(hlls[1],"Background");
  leg1->AddEntry(hlls[2],"Data");
  leg1->Draw();
  
  TCanvas *canLLr=new TCanvas("canLLr","Log Likelihood Ratio",1000,600);
  hllr[0]->Draw("hist");
  hllr[1]->Draw("hist,same");
  hllr[2]->Draw("hist,same");
  TLegend* leg2=new TLegend(0.125,0.7,0.425,0.9);
  leg2->AddEntry(hllr[0],"Signal");
  leg2->AddEntry(hllr[1],"Background");
  leg2->AddEntry(hllr[2],"Data");
  leg2->Draw();

  TCanvas *canLL2=new TCanvas("canLL2","Log Likelihood Ratio 2D",1000,600);
  hll2[0]->Draw("hist");
  hll2[1]->Draw("hist,same");
  hll2[2]->Draw("hist,same");
  TLegend* leg3=new TLegend(0.125,0.7,0.425,0.9);
  leg3->AddEntry(hll2[0],"Signal");
  leg3->AddEntry(hll2[1],"Background");
  leg3->AddEntry(hll2[2],"Data");
  leg3->Draw();

  // 7. save to file
  for (int i=0; i<3; i++) {
    hlls[i]->Write();
    hllr[i]->Write();
    hll2[i]->Write();
  }
  //fOut->Close();
}

// Find the min/max values of some variable is present in two trees
void minmax(double &min, double &max, TTree *s, TTree *b, TString var) {
  min=TMath::Min(s->GetMinimum(var),b->GetMinimum(var));
  max=TMath::Max(s->GetMaximum(var),b->GetMaximum(var));
}

// Create histograms with reasonable x-axis ranges
void makeHists(TH1F **hs, TH1F **hb, TTree *s, TTree *b,TString var) {
  double min,max;
  minmax(min,max,s,b,var);

  *hs=new TH1F("hs_"+var,"Signal: "+var,100,min,max);
  (*hs)->SetLineColor(kRed);
  (*hs)->Sumw2();

  *hb=new TH1F("hb_"+var,"Background: "+var,100,min,max);
  (*hb)->SetLineColor(kBlue);
  (*hb)->Sumw2();
}

void makeHists2(TH2F **hs2, TH2F **hb2, TTree *s, TTree *b, TString varX, TString varY) {
  double minX,maxX;
  minmax(minX,maxX,s,b,varX);
  double minY,maxY;
  minmax(minY,maxY,s,b,varY);

  *hs2=new TH2F("hs2_"+varX+"_vs_"+varY,"Signal: "+varX+", "+varY,100,minX,maxX,100,minY,maxY);
  (*hs2)->SetLineColor(kRed);
  (*hs2)->Sumw2();

  *hb2=new TH2F("hb2_"+varX+"_vs_"+varY,"Signal: "+varX+", "+varY,100,minX,maxX,100,minY,maxY);
  (*hb2)->SetLineColor(kBlue);
  (*hb2)->Sumw2();
}