#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include "TLatex.h"

void setAxes(TH1F* h) {
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);

  h->GetYaxis()->SetTitle("Count");
}

// fit1a.C
void fit1a(int entries=1000, bool save=false) {
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram", 100, 0, 100);
  TRandom2 *generator = new TRandom2(0);  // parameter == seed, 0->use clock
  
  int ntrials=1000;
  TF1 *fitFunc;
  TH1F *fitChi2Hist = new TH1F("fitChi2Hist", "Fit #chi^{2} Distribution", 50,0,120);
  TH1F *fitMeanHist = new TH1F("fitMeanHist", "Fit Mean Distribution", 40,48,52);
  TH1F *fitProbHist = new TH1F("fitProbHist", "Fit Probability Distribution", 50,0,1.5);
  TGraph *reduced = new TGraph();
  reduced->SetName("reduced");
  reduced->SetTitle("Fit Probability vs. Reduced #chi^{2}");

  for (int i=0; i<ntrials; i++) {
    randomHist1->Reset(); // reset histogram bin content to 0

    for (int i=0; i<entries; i++) {
      randomHist1->Fill(generator->Gaus(50,10)); // params: mean, sigma
    }
    randomHist1->Fit("gaus","q");
    fitFunc = randomHist1->GetFunction("gaus");

    fitChi2Hist->Fill(fitFunc->GetChisquare());
    fitMeanHist->Fill(fitFunc->GetParameter(1));
    fitProbHist->Fill(fitFunc->GetProb());
    reduced->SetPoint(i,fitFunc->GetChisquare()/fitFunc->GetNDF(),fitFunc->GetProb());
  }

  // set up canvas
  TCanvas *canvas = new TCanvas("fi1a", "Non-Linear Fit", 1400, 700);
  canvas->Divide(2,2);
  //gStyle->SetOptStat("M");
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.25);
  gStyle->SetTitleSize(0.075,"t");
  gPad->SetLeftMargin(0.125);
  gPad->SetBottomMargin(0.125);

  // make the plots prettier
  setAxes(fitChi2Hist);
  setAxes(fitMeanHist);
  setAxes(fitProbHist);
  fitChi2Hist->GetXaxis()->SetTitle("#chi^{2}");
  fitMeanHist->GetXaxis()->SetTitle("Fit Mean");
  fitProbHist->GetXaxis()->SetTitle("Fit Probability");

  reduced->GetXaxis()->SetLabelSize(0.05);
  reduced->GetYaxis()->SetLabelSize(0.05);
  reduced->GetXaxis()->SetTitleSize(0.05);
  reduced->GetYaxis()->SetTitleSize(0.05);
  reduced->GetXaxis()->SetTitle("Reduced #chi^{2}");
  reduced->GetYaxis()->SetTitle("Fit Probability");
  
  // draw plots
  canvas->cd(1);
  fitChi2Hist->Draw("");
  canvas->cd(2);
  fitMeanHist->Draw("");
  canvas->cd(3);
  fitProbHist->Draw("");
  canvas->cd(4);
  reduced->Draw("AP");
  
  // prints 
  canvas->Print("fits1.pdf","pdf");

  TFile *tf=0;
  if (save) {
    tf = new TFile("histo.root","recreate");
    tf->Write();
    tf->Close();
  }
}
