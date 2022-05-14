THStack *Stack(TH1F *h1, TH1F *h2);


// this is basically an ad hoc definition of the signal region
// defined by looking at the 2D S and B distrinbutions
float select(Float_t *var){
  bool pass=true;
  pass=var[0]>50 && var[0]<250;
  pass=pass && var[1]>0.2;
  pass=pass && var[1]<1.2*var[3]-25 && var[1]>1.3*var[3]-112.5;
  pass=pass && var[3]>60 and var[4]>-7.5;

  if (pass) return 1.0;
  return 0;
}


// In this example we assume a S/B of 0.02 (default), with a total of 10200
// events and calculate the signal significance before and after applying
// some 2D selections

void select2d(float ns=200., float nb=10000.){
  //gStyle->SetOptStat(0);
  TFile *in=new TFile("simulation.root");
  TTree *s=(TTree*)in->Get("TreeS");
  TTree *b=(TTree*)in->Get("TreeB"); 
  float snorm=1.0*ns/s->GetEntries();
  float bnorm=1.0*nb/b->GetEntries();
  cout << "Initial signal significance " << ns/TMath::Sqrt(ns+nb) << endl;
  Float_t bvar[5];
  Float_t svar[5];
  // the static  TString::Format method return a TString following
  // a printf style format
  for (int i=1; i<=5; i++){
    b->SetBranchAddress(TString::Format("var%d",i),&bvar[i-1]);
    s->SetBranchAddress(TString::Format("var%d",i),&svar[i-1]);
  }

  // here we make some histograms to study distributions
  // before and after the cuts are made
  TH1F *var1s=new TH1F("var1s","var1s",100,0,400);
  TH1F *var2s=new TH1F("var2s","var2s",100,0,120);
  TH1F *var3s=new TH1F("var3s","var3s",100,0,1);
  TH1F *var4s=new TH1F("var4s","var4s",100,0,400);
  TH1F *var5s=new TH1F("var5s","var5s",100,-30,20);
  TH1F *var1b=new TH1F("var1b","var1b",100,0,400);
  TH1F *var2b=new TH1F("var2b","var2b",100,0,120);
  TH1F *var3b=new TH1F("var3b","var3b",100,0,1);
  TH1F *var4b=new TH1F("var4b","var4b",100,0,400);
  TH1F *var5b=new TH1F("var5b","var5b",100,-30,20);
  TH1F *var1cs=new TH1F("var1cs","var1cs",100,0,400);
  TH1F *var2cs=new TH1F("var2cs","var2cs",100,0,120);
  TH1F *var3cs=new TH1F("var3cs","var3cs",100,0,1);
  TH1F *var4cs=new TH1F("var4cs","var4cs",100,0,400);
  TH1F *var5cs=new TH1F("var5cs","var5cs",100,-30,20);
  TH1F *var1cb=new TH1F("var1cb","var1cb",100,0,400);
  TH1F *var2cb=new TH1F("var2cb","var2cb",100,0,120);
  TH1F *var3cb=new TH1F("var3cb","var3cb",100,0,1);
  TH1F *var4cb=new TH1F("var4cb","var4cb",100,0,400);
  TH1F *var5cb=new TH1F("var5cb","var5cb",100,-30,20);

  // below we'll access the histograms by name to help with automation
  TH1F* h;
  for (int n=0;n<s->GetEntries();n++){ // loop over signal model
    s->GetEntry(n);
    for (int i=1; i<=5; i++){
      h=(TH1F*)gDirectory->Get(TString::Format("var%ds",i));
      h->Fill(svar[i-1],snorm);  // fill before cuts
      h->SetLineColor(kRed);
      h->SetFillColor(kRed);
      h->SetFillStyle(1002);
      h=(TH1F*)gDirectory->Get(TString::Format("var%dcs",i));
      h->Fill(svar[i-1],snorm*select(svar));  // fill after cuts
      h->SetLineColor(kRed);
      h->SetFillColor(kRed);
      h->SetFillStyle(1002);
    }   
  }
  for (int n=0;n<b->GetEntries();n++){ // loop over bkg model
    b->GetEntry(n);
    for (int i=1; i<=5; i++){
      h=(TH1F*)gDirectory->Get(TString::Format("var%db",i));
      h->Fill(bvar[i-1],bnorm);  // fill before cuts
      h->SetLineColor(kBlue);
      h->SetFillColor(kBlue);
      h->SetFillStyle(1002);
      h=(TH1F*)gDirectory->Get(TString::Format("var%dcb",i));
      h->Fill(bvar[i-1],bnorm*select(bvar));  // fill after cuts
      h->SetLineColor(kBlue);
      h->SetFillColor(kBlue);
      h->SetFillStyle(1002);
    }
  }

  // plot stacks of signal and bkg 
  TCanvas *csb2=new TCanvas("csb2","csb2",1200,600);
  TH1F *hs, *hb;
  csb2->Divide(5,2);
  for (int i=1;i<=10;i++){
    csb2->cd(i);
    if (i<6) {
      hs=(TH1F*)gDirectory->Get(TString::Format("var%ds",i));
      hb=(TH1F*)gDirectory->Get(TString::Format("var%db",i));
      Stack(hb,hs)->Draw("hist");
    }
    else  {
      hs=(TH1F*)gDirectory->Get(TString::Format("var%dcs",i-5));
      hb=(TH1F*)gDirectory->Get(TString::Format("var%dcb",i-5));
      Stack(hb,hs)->Draw("hist");
    }
  }
  float nscut=var1cs->Integral();
  float nbcut=var1cb->Integral();

  cout << "Signal significance after cuts " <<
    nscut/TMath::Sqrt(nscut+nbcut) << endl;

  cout << "stacked variable distributions: pre-cuts (top row), post-cuts (bottom row)" << endl;
}


THStack *Stack(TH1F *h1, TH1F *h2){
    THStack *stack=new THStack();
    stack->Add(h1);
    stack->Add(h2);
    return stack;
} 
