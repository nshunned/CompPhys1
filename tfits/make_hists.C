// make_hists.C

#include "TROOT.h"   // needed for gROOT class
#include "TH1F.h"    // defines TH1F class
#include "TF1.h"
#include "TFile.h"   // defines TFILE class
#include "TRandom2.h" // ROOT random number generator
#include <iostream>


void make_hists(){
  gROOT->Reset();  // useful to reset ROOT to a cleaner state
  TFile *file = new TFile("mydata.root","recreate");
  TRandom *generator=new TRandom2(0);  // parameter == seed, 0->use clock
  TH1F *back1 = new TH1F("back1", "Bkg Histogram1", 50, 0, 50);
  TH1F *back2 = new TH1F("back2", "Bkg Histogram2", 50, 0, 50);
  TH1F *signal=new TH1F("signal", "Signal Histogram", 50, 0, 50);
  TH1F *data=new TH1F("data", "Data Histogram;x;# of Events", 50, 0, 50);

  const Int_t NMC=10000;
  const Int_t NDATA=1000;

  TF1 *f1 = new TF1("f1","pow((x/50.),3)",0,50);
  TF1 *f2 = new TF1("f2","1/sqrt(x+1)",0,50);


  for (int i=0;i<NMC;i++) {
    signal->Fill(generator->Landau(30,3));
    back1->Fill(f1->GetRandom()); 
    back2->Fill(f2->GetRandom()); 
  }

  Float_t b1=100-generator->Rndm()*10;
  Float_t b2=100-generator->Rndm()*10;

  for (int i=0;i<NDATA;i++) {
    data->Fill(generator->Landau(30,3));
    if (generator->Rndm()*100<b1)
      data->Fill(f1->GetRandom()); 
    if (generator->Rndm()*100<b2)
      data->Fill(f2->GetRandom()); 
  }

  file->Write();  // writes all objects to the file
  file->Close();  // comment out to prevent histogram display from being erased
  std::cout << "Created mydata.root" << endl;
}
