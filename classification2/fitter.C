#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFractionFitter.h>

#include <iostream>
using namespace std;

void fitter(TH1F* mbkg, TH1F *msig, TH1F *data){
   gStyle->SetOptStat(kFALSE); 
   TCanvas *tc=new TCanvas("tcfit");
   
   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   mc->Add(msig);
   mc->Add(mbkg);
   data->SetLineWidth(2);
   data->SetLineColor(kBlack);

   TFractionFitter* fit = new TFractionFitter(data, mc, "V"); // initialize
   fit->Constrain(1,0.0,1.0);             // constrain bkg fraction to be between 0 and 1
   //fit->Constrain(0,0.0,1.0);
   //fit->SetRangeX(1,15);                // use only the first 15 bins in the fit
   Int_t status = fit->Fit();               // perform the fit
   cout << "fit status: " << status << endl;
   
   TH1F* result = (TH1F*) fit->GetPlot();
   if (status == 0) {                       // check on fit status
     result->SetLineWidth(3);
     data->Draw("Ep");
     result->Draw("same");
   }

   //return;

   Double_t params[2], errors[2];
   Double_t ndata=data->Integral();
   fit->GetResult(0,params[0],errors[0]);
   fit->GetResult(1,params[1],errors[1]);
   cout <<  "signal/background fraction: " << params[0] << "  / " << params[1] << endl;

   TH1F *hSfit=(TH1F*)fit->GetMCPrediction(0);
   TH1F *hBfit=(TH1F*)fit->GetMCPrediction(1);
   // fitted signal
   hSfit->Scale( 1./hSfit->Integral() * params[0]*ndata );
   hSfit->SetLineColor(kRed);
   hSfit->SetLineWidth(2);
   //hSfit->SetFillColor(kRed);
   // fitted background 
   hBfit->Scale( 1./hBfit->Integral() * params[1]*ndata );
   hBfit->SetLineColor(kBlue);
   hBfit->SetLineWidth(2);
   //hBfit->SetFillColor(kBlue);

   double mu= hSfit->Integral(); // # of signal events
   double sigma = hSfit->Integral()*errors[0]/params[0]; // fit error on # signal events
   cout << "Number of data events: " << data->Integral() << endl;
   cout << "Background Events: " << hBfit->Integral() << endl;
   cout << "Signal Events: " << mu << " +- "  
	<<  sigma << "  mu/sigma= "  << mu/sigma << endl;

   hBfit->Draw("same,hist");
   hSfit->Draw("same,hist");

   auto tl=new TLegend(0.6,0.7,0.9,0.9);
   tl->AddEntry(fit->GetPlot(),"S+B fit","l");
   tl->AddEntry(hBfit,"Bkg fit","l");
   tl->AddEntry(hSfit,"Sig fit","l");
   tl->Draw();

   
   /*
   THStack *hs=new THStack("hs","stacked MC histograms");
   hs->Add(h0);
   hs->Add(h1);

   data->SetMinimum(0.1);
   data->Draw("Ep");
   hs->Draw("same");
   result->Draw("same");
   data->Draw("Ep,same");
   */

 }
