#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TRandom.h"
#include "TNamed.h"
#include "TAxis.h"
#include "math.h"

void fifth(){
    // file to write out our plots
    TFile* file = new TFile("sixplots.root","recreate");

    // canvas to draw our plots
    TCanvas* canvas = new TCanvas("Six Plots", "Six Plots", 1350, 700);
    canvas->Divide(3,2);
    gStyle->SetOptStat("0");
    gStyle->SetTitleSize(0.075,"t");
    gPad->SetLeftMargin(0.125);
    gPad->SetBottomMargin(0.125);

    // a random generator
    TRandom* generator=new TRandom(0);  // parameter == seed, 0->use clock

    // plot #1
    canvas->cd(1);
    TH1F* one = new TH1F("one", "Normally Distributed",400,0,400);
    one->GetXaxis()->SetLabelSize(0.04);
    one->GetYaxis()->SetLabelSize(0.04);
    one->GetYaxis()->SetTitle("Count");
    one->GetYaxis()->SetTitleSize(0.04);

    for (int i=0; i<10000; i++)
        one->Fill(generator->Gaus(200,15));
    one->Draw();
    one->Write();

    // plot #2
    canvas->cd(2);
    TH1F* two = new TH1F("two", "Uniformly Distributed",400,0,400);
    two->GetXaxis()->SetLabelSize(0.04);
    two->GetYaxis()->SetLabelSize(0.04);
    two->GetYaxis()->SetTitle("Count");
    two->GetYaxis()->SetTitleSize(0.04);

    for (int i=0; i<10000; i++)
        two->Fill(generator->Integer(400));
    two->Draw();
    two->Write();

    // plot #3
    canvas->cd(3);
    TH1F* three = new TH1F(*one);
    three->SetNameTitle("three", "Uniform + Normal");
    three->GetXaxis()->SetLabelSize(0.04);
    three->GetYaxis()->SetLabelSize(0.04);
    three->GetYaxis()->SetTitle("Count");
    three->GetYaxis()->SetTitleSize(0.04);

    three->Add(two);
    three->Draw();
    three->Write();

    // plot #4
    canvas->cd(4);
    TH1F* four = new TH1F("four", "Negatively Natural", 100,0,10);
    four->GetXaxis()->SetLabelSize(0.04);
    four->GetYaxis()->SetLabelSize(0.04);
    four->GetYaxis()->SetTitle("Weighted Count");
    four->GetYaxis()->SetTitleSize(0.04);

    for(int i=0; i<four->GetNbinsX(); i++){
        double center = four->GetXaxis()->GetBinCenter(i);
        four->Fill(center, exp(-center));
    }
    four->Draw();
    four->Write();

    // plot #5
    canvas->cd(5);
    TH2F* five = new TH2F("five", "2D Normal Distribution", 200,-100,100,200,-100,100);
    for(int i=0; i<1000; i++)
        five->Fill(generator->Gaus(0,10),generator->Gaus(1,20));
    five->Draw();
    five->Write();

    // plot #6
    canvas->cd(6);
    double max[] = {84, 90, 87, 94, 79, 85, 86, 88, 85, 88, 97, 97, 76, 79, 86, 93, 76, 78, 74, 82, 91, 93, 93, 82, 83, 92, 83, 90, 91, 76};
    double avg[] = {76.2, 75.8, 76.4, 79.6, 75.3, 74.2, 73.8, 74.1, 75, 75, 78.3, 81.5, 70.8, 71.6, 75.4, 77.6, 71.9, 67.7, 63.4, 65.1, 74.4, 78.9, 80.8, 73.2, 68.9, 73.1, 71.4, 73.1, 73.2, 73.3};
    double min[] = {70, 70, 69, 67, 72, 69, 63, 60, 67, 67, 68, 68, 67, 67, 69, 63, 67, 60, 55, 50, 58, 66, 68, 62, 54, 59, 61, 65, 67, 71};

    TGraph* maximum = new TGraph();
    TGraph* average = new TGraph();
    TGraph* minimum = new TGraph();

    maximum->SetName("Maximum");
    average->SetName("Average");
    minimum->SetName("Minimum");
    
    for(int i=1; i<=30; i++){
        maximum->SetPoint(i-1, i, max[i-1]);
        average->SetPoint(i-1, i, avg[i-1]);
        minimum->SetPoint(i-1, i, min[i-1]);
    }

    maximum->SetTitle("September");
    maximum->GetYaxis()->SetRangeUser(40, 100);
    maximum->GetYaxis()->SetTitle("Temperature (Deg F)");
    maximum->GetXaxis()->SetTitle("Day");
    maximum->GetXaxis()->SetTitleSize(0.04);
    maximum->GetYaxis()->SetTitleSize(0.04);

    maximum->SetLineColor(2);
    average->SetLineColor(3);
    minimum->SetLineColor(4);
    maximum->SetLineWidth(2);
    average->SetLineWidth(2);
    minimum->SetLineWidth(2);

    maximum->Draw("");
    average->Draw("SAME");
    minimum->Draw("SAME");
    maximum->Write();
    average->Write();
    minimum->Write();

    TLegend* legend = new TLegend(0.5, 0.32, 0.16, 0.15, "");
    legend->SetTextSize(0.05);
    legend->AddEntry(maximum, "Maximum");
    legend->AddEntry(average, "Average");
    legend->AddEntry(minimum, "Minimum");
    legend->Draw();
  
    canvas->Print("sixplots.pdf","pdf");
}