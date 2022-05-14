// forth.C
// Here we specify a named function for ROOT.
// The function is meant to be compiled into a shared 
// library that will be available to the ROOT framework
// Usage in ROOT: .L forth.C+
//                randomHist2() or randomHist2(int)
// Note: you must reload/compile the function after making changes!
//

// include files are now needed b/c this code will be compiled w/ the 
// standard c++ compiler

#include "TROOT.h"   // needed for gROOT class
#include "TH1F.h"    // defines TH1F class
#include "TFile.h"   // defines TFILE class
#include "TRandom.h" // ROOT random number generator
#include <iostream>
using namespace std;

// Here we explicitly access the random numebr generator class in ROOT
// and use Fill to add each entry to the histogram

void randomHist4(int entries=1000){
  gROOT->Reset();  // useful to reset ROOT to a cleaner state
  TFile *file = new TFile("randomHist4.root","recreate");
  TRandom *generator=new TRandom(0);  // parameter == seed, 0->use clock
  TH1F *randomHist4 = new TH1F("randomHist4", "Random Histogram", 100, 0, 100);
  for (int i=0 ; i<entries ; i++){
    randomHist4->Fill(generator->Gaus(50,10)); // params: mean, sigma
  }
  // randomHist4->Draw(); // closing the file will blank out the histogram display
  file->Write();  // writes all objects to the file
  cout << "Wrote file: " << file->GetName() << endl;
  file->Close();
}
