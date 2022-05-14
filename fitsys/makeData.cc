// generate data with systematics
// "True" (parent) distribution is a linear function
// First sample the parent pdf
// Then add a systematic that affects the shape
// Finally add a systematic that shifts the overall normaliztion

void makeData(){
  gStyle->SetOptStat(0);  // turn off histogram stat box
  TFile * tf=new TFile("mydata.root","recreate");

  TCanvas *c1=new TCanvas("ca","My data",800,800);
  c1->Divide(2,2);

  c1->cd(1);
  int nevents=10000;
  //TF1 *fun1=new TF1("f1","[0]*(2*x+5)/750.",0,25);   // parent pdf, p0=66 2/3, p1=26 2/3 
  TF1 *fun1=new TF1("f1","[0]*(0.75*x+3.75)/750.",0,25); // p0=50, p1=10
  //TF1 *fun1=new TF1("f1","[0]*(1*x+2.5)/750.",0,25); // p0=50, p1=10
  fun1->SetParameter(0,1);
  TH1F *h1=new TH1F ("h1","Sample, no systematics;x;# events",25,0,25);
  h1->FillRandom("f1",nevents);
  h1->Draw("e");
  fun1->SetLineColor(kGreen);
  fun1->SetParameter(0,nevents);
  fun1->Draw("same");

  c1->cd(2);
  TF1 *fun2=new TF1("f2","6*x*x/TMath::Power(1.4,x-0.5)",0,25);  // syst 1
  fun2->SetTitle("Systematic 1;x;#Delta events");
  fun2->SetLineColor(1);
  fun2->Draw();

  TH1F *h2=new TH1F(*h1);  // copy of histogram w/ systematics to be added
  h2->SetTitle("Sample w/ systematics");
  h2->SetName("Data");
  TRandom *tr = new TRandom();  // get a random # generator
  tr->SetSeed();
  double sigma=6;    // scale for systematic 1
  double shift=tr->Gaus(0,1);
  h2->Add(fun2,shift);  // shift by # S.D.
  h2->SetLineColor(2);
  
  // systematic 2 is an offset
  c1->cd(3);
  shift=tr->Gaus(0,1);  // # S.D. for offset shift
  TF1 *fun3=new TF1("f3","20",0,25);
  fun3->SetTitle("Systematic 2;x;#Delta events");
  h2->Add(fun3,shift);
  fun3->Draw();

  c1->cd(4);
  h2->Draw("e");
  fun1->Draw("same");
    

  // make error band around true pdf
  TH1F *hup=new TH1F(*h1); hup->Reset();
  TH1F *hdn=new TH1F(*h1); hdn->Reset();

  // uncorrelated systeatics add in quadrature
  for (int i=1; i<=hup->GetNbinsX();i++){
    double x=hup->GetBinCenter(i);
    double y=fun1->Eval(x);
    double s=TMath::Sqrt(fun2->Eval(x)*fun2->Eval(x)+fun3->Eval(x)*fun3->Eval(x));
    hup->SetBinContent(i,y+s);
    hdn->SetBinContent(i,y-s);
  }
  
  hup->Draw("same");
  hdn->Draw("same");

  c1->Print("ca.png");
  // write the data and the histograms showing 1 S.D. systematics
  TH1F *hs1=new TH1F(*h1);  hs1->Reset(); hs1->SetName("syst1");
  hs1->Add(fun2);
  TH1F *hs2=new TH1F(*h2);  hs2->Reset(); hs2->SetName("syst2");
  hs2->Add(fun3);

  h2->Write();
  hs1->Write();
  hs2->Write();
  // tf->Close();  // keep hte file open in order to see the plots
}

