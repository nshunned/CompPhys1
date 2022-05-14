void plot_inputs(){
  gStyle->SetOptStat(0);

  bool nostack=0;

  TFile *tf=new TFile("mydata.root");
  
  TH1F *back1=(TH1F*)tf->Get("back1");   // model for background 1
  TH1F *back2=(TH1F*)tf->Get("back2");   // model for background 2
  TH1F *signal=(TH1F*)tf->Get("signal"); // model for signal
  TH1F *data=(TH1F*)tf->Get("data");     // some data (combo of 3 above distros)

  back1->SetLineColor(kGreen);
  back2->SetLineColor(kBlue);
  signal->SetLineColor(kRed);
  back1->SetLineWidth(3);
  back2->SetLineWidth(3);
  signal->SetLineWidth(3);
  data->SetLineWidth(3);

  if (nostack==false){
    back1->SetFillColor(kGreen);
    back2->SetFillColor(kBlue);
    signal->SetFillColor(kRed);
  }

  // here we plot the histograms read from the file
  TCanvas *c1=new TCanvas("c1","Inputs",5,20,600,500);
  c1->Divide(2,2);
  c1->cd(1);
  back1->Draw();
  c1->cd(2);
  back2->Draw();
  c1->cd(3);
  signal->Draw();
  c1->cd(4);
  data->Draw("e");

  // this in an example of how to stack the histograms
  // your job is to use TFractionFitter to find the fractions
  // of each component, necessary to reproduce the data
  TCanvas *c2=new TCanvas("c2","stack",605,20,600,500);
  THStack *hs=new THStack("hs","stacked MC histograms");
  hs->Add(back1);
  hs->Add(back2);
  hs->Add(signal);
  if (nostack)
    hs->Draw("nostack");
  else
    hs->Draw();

  TLegend *leg = new TLegend(0.1,0.7,0.45,0.9);
  leg->AddEntry(signal,"Signal");
  leg->AddEntry(back1,"Background 1");
  leg->AddEntry(back2,"Background 2");

  leg->Draw();
}
