// third.C
// Here we specify a named function for ROOT.
// The function is meant to be compiled into a shared 
// library that will be available to the ROOT framework
// Usage in ROOT: .L third.C+
//                randomHist2() or randomHist2(int)
// Note: you must reload/compile the function after making changes!
//

// include files are now needed b/c this code will be compiled w/ the 
// standard c++ compiler

#include "TROOT.h"   // needed for gROOT class
#include "TH1F.h"    // defines TH1F class

void randomHist2(int entries=1000){
  gROOT->Reset();  // useful to reset ROOT to a cleaner state [optional]
  TH1F *my_hist = new TH1F("my_hist", "My First Histogram", 100, -5, 5);
  my_hist->FillRandom("gaus", entries);
  my_hist->Draw();
}
