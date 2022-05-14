#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include "TLatex.h"

using namespace TMath;

void setCanvas() {
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.25);
  gStyle->SetTitleSize(0.075,"t");
  gPad->SetLeftMargin(0.125);
  gPad->SetBottomMargin(0.125);
}

void setAxes(TH1F* h) {
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
}

double nLL0(TH1F* hist, TF1* fit) {
  double nll=0;
    for (int i=0; i<hist->GetNbinsX(); i++) {
      double expected = fit->Eval(hist->GetBinCenter(i));
      double nObs = hist->GetBinContent(i);
      double prob = PoissonI(nObs,expected);
      nll += Log(prob);
    }
  return nll *= -1;
}

double nLL(TH1F *hist, TF1 *fit) {
  double nll=0;
  TRandom2* r = new TRandom2();
  for (int i=0; i<hist->GetNbinsX(); i++) {
    double expected = fit->Eval(hist->GetBinCenter(i));
    int nObs = r->Poisson(expected);
    double prob = PoissonI(nObs,expected);
    nll += Log(prob);
  }
  return -nll;
}

void partI(int entries) {
  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram", 100, 0, 100);
  TRandom2 *generator = new TRandom2(0);  // parameter == seed, 0->use clock
  
  int ntrials=1000;
  TF1 *fitFuncChi2;
  TF1 *fitFuncLog;
  TH1F *fitMeanChi2Hist = new TH1F("fitMeanChi2Hist", "#chi^{2} Fit Mean;Fit Mean;Count", 40,48,52);
  TH1F *fitMeanLogHist = new TH1F("fitMeanLogHist", "Log Likelihood Fit Mean;Fit Mean;Count", 40,48,52);
  TH1F *fitSigmaChi2Hist = new TH1F("fitSigmaChi2Hist", "#chi^{2} Fit Sigma;Fit Sigma;Count", 60,8,12);
  TH1F *fitSigmaLogHist = new TH1F("fitSigmaLogHist", "Log Likelihood Fit Sigma;Fit Sigma;Count", 60,8,12);

  for (int i=0; i<ntrials; i++) {
    // reset histogram bin content to 0
    randomHist1->Reset();
    // fill histogram
    for (int i=0; i<entries; i++) {
      randomHist1->Fill(generator->Gaus(50,10)); // params: mean, sigma
    }

    // fit using Chi2
    randomHist1->Fit("gaus","0q");
    fitFuncChi2 = randomHist1->GetFunction("gaus");
    // remove previous fit & refit using LL
    randomHist1->GetListOfFunctions()->Remove(fitFuncChi2);
    randomHist1->Fit("gaus","0qL");
    fitFuncLog = randomHist1->GetFunction("gaus");

    // fill parameter histograms
    fitMeanChi2Hist->Fill(fitFuncChi2->GetParameter(1));
    fitMeanLogHist->Fill(fitFuncLog->GetParameter(1));
    fitSigmaChi2Hist->Fill(fitFuncChi2->GetParameter(2));
    fitSigmaLogHist->Fill(fitFuncLog->GetParameter(2));
  }

  // set up canvas
  TCanvas *canvas = new TCanvas("PartI", "Non-Linear Fit", 1400,700);
  setCanvas();
  canvas->Divide(2,2);

  // make the plots prettier
  setAxes(fitMeanChi2Hist);
  setAxes(fitMeanLogHist);
  setAxes(fitSigmaChi2Hist);
  setAxes(fitSigmaLogHist);

  // draw plots
  canvas->cd(1); fitMeanChi2Hist->Draw("");
  canvas->cd(2); fitMeanLogHist->Draw("");
  canvas->cd(3); fitSigmaChi2Hist->Draw("");
  canvas->cd(4); fitSigmaLogHist->Draw("");
  
  // print canvas
  canvas->Print("fit1b.pdf","pdf");
  // canvas->Close();
}

void partII() {
  TCanvas *canvas = new TCanvas("PartII", "Non-Linear Fit", 700,400);
  setCanvas();

  TRandom2 *generator = new TRandom2(0);  // parameter == seed, 0->use clock

  // extract histogram from root file
  TFile *file = new TFile("histo25.root");
  TH1F *randomHist = new TH1F(*(TH1F*)file->Get("randomHist1"));
  
  // fit & extract fit, NLL
  randomHist->Fit("gaus","0qL");
  TF1 *fitFunc = randomHist->GetFunction("gaus");
  double nll=nLL0(randomHist,fitFunc);

  // pseudo-experiments & histogram
  int nTrials = 1000;
  int nEntries = randomHist->GetEntries();
  TH1F *nLLHist = new TH1F("nLLHist", "Negative Log Likelihood;NLL;Count", 80,nll-15,nll+25);

  double mean = fitFunc->GetParameter(1);
  double sigma = fitFunc->GetParameter(2);
  for (int i=0; i<nTrials; i++) {
      // reset histogram bin content to 0
      randomHist->Reset();
      // fill histogram
      for (int i=0; i<nEntries; i++) {
        randomHist->Fill(generator->Gaus(mean,sigma)); // params: mean, sigma
      }

      // fit pseudo-experiment histogram using LL
      randomHist->Fit("gaus","0qL");
      fitFunc = randomHist->GetFunction("gaus");

      // fill NLL histogram
      nLLHist->Fill(nLL(randomHist,fitFunc));
  }
  setAxes(nLLHist);
  nLLHist->Draw();

  TLine *nLLLine = new TLine(nll,0,nll,nLLHist->GetBinContent(nLLHist->GetXaxis()->FindBin(nll)));
  nLLLine->SetLineWidth(2);
  nLLLine->SetLineColor(2);
  nLLLine->Draw();

  canvas->Print("fit1b_nll.pdf","pdf");
  //canvas->Close();
  // file->Close();
}

void partIII() {
  // NLL & chi2 contours
  TGraph *nLLContour = new TGraph();
  TGraph *chi2Contour = new TGraph();

  // extracting histograms
  TFile *file = new TFile("histo25.root");
  TH1F *randomHist1 = new TH1F(*(TH1F*)file->Get("randomHist1"));
  file->Close();
  file = new TFile("histo1k.root");
  TH1F *randomHist2 = new TH1F(*(TH1F*)file->Get("randomHist1"));
  file->Close();

  // fit & extract fit, mean
  randomHist1->Fit("gaus","0qL");
  randomHist2->Fit("gaus","0q");
  TF1 *fitFunc1 = randomHist1->GetFunction("gaus");
  TF1 *fitFunc2 = randomHist2->GetFunction("gaus");
  double mean1 = fitFunc1->GetParameter(1);
  double mean2 = fitFunc2->GetParameter(1);
  
  int nSteps=500;
  for (int i=-nSteps; i<nSteps; i++) {
    // vary fit mean
    double newMean1 = mean1+i*0.1;
    double newMean2 = mean2+i*0.1;
    fitFunc1->FixParameter(1,newMean1);
    fitFunc2->FixParameter(1,newMean2);

    // remove fit & refit
    randomHist1->GetListOfFunctions()->Remove(fitFunc1);
    randomHist2->GetListOfFunctions()->Remove(fitFunc2);
    randomHist1->Fit(fitFunc1,"0qBL");
    randomHist2->Fit(fitFunc2,"0qB");

    fitFunc1 = randomHist1->GetFunction("gaus");
    fitFunc2 = randomHist2->GetFunction("gaus");
    nLLContour->SetPoint(i+nSteps, newMean1, nLL0(randomHist1,fitFunc1));
    chi2Contour->SetPoint(i+nSteps, newMean2, fitFunc2->GetChisquare());
  }

  nLLContour->SetTitle("NLL Contour");
  nLLContour->GetXaxis()->SetLabelSize(0.05);
  nLLContour->GetYaxis()->SetLabelSize(0.05);
  nLLContour->GetXaxis()->SetTitleSize(0.05);
  nLLContour->GetYaxis()->SetTitleSize(0.05);
  nLLContour->GetXaxis()->SetTitle("Fit Mean");
  nLLContour->GetYaxis()->SetTitle("Fit NLL");
  nLLContour->GetYaxis()->SetTitleOffset(0.9);
  
  chi2Contour->SetTitle("#chi^{2} Contour");
  chi2Contour->GetXaxis()->SetLabelSize(0.05);
  chi2Contour->GetYaxis()->SetLabelSize(0.05);
  chi2Contour->GetXaxis()->SetTitleSize(0.05);
  chi2Contour->GetYaxis()->SetTitleSize(0.05);
  chi2Contour->GetXaxis()->SetTitle("Fit Mean");
  chi2Contour->GetYaxis()->SetTitle("Fit #chi^{2}");
  chi2Contour->GetYaxis()->SetTitleOffset(0.8);

  // set up new canvas
  TCanvas *canvas = new TCanvas("PartIII", "Non-Linear Fit", 1400,400);
  setCanvas();
  canvas->Divide(2,1);

  // draw contours
  canvas->cd(1); nLLContour->Draw("");
  canvas->cd(2); chi2Contour->Draw("");
  
  canvas->Print("fit1b_contours.pdf","pdf");
  //canvas->Close();
}

// fit1a.C
void fit1b(int entries=1000, bool save=false) {
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  partI(entries);
  partII();
  partIII();
  
  TFile *tf=0;
  if (save) {
    tf = new TFile("histo.root","recreate");
    tf->Write();
    tf->Close();
  }
}
