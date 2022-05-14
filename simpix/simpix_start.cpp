// simple example of using ROOT libraries in a C++ program with graphics
// and use of TASImage class

#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TRandom1.h"

#include "assert.h"

#include <ctime>
#include <cstdio>
#include <iostream>

using namespace std;
using namespace TMath;

TRandom *gen = new TRandom();

Double_t GetDistance(UInt_t pix1, UInt_t pix2) {
  // rgb for pix1
  UInt_t r1 = pix1 & 0xffff0000;
  UInt_t g1 = pix1 & 0xff00ff00;
  UInt_t b1 = pix1 & 0xff0000ff;
  // rgb for pix2
  UInt_t r2 = pix2 & 0xffff0000;
  UInt_t g2 = pix2 & 0xff00ff00;
  UInt_t b2 = pix2 & 0xff0000ff;
  return sqrt( (Double_t) ((r1-r2)*(r1-r2) + (g1-g2)*(g1-g2) + (b1-b2)*(b1-b2)) );
}

Double_t GetTotalDistance(Long_t numPix, UInt_t* tgtPix, UInt_t* outPix) {
  Double_t distance=0;
  for (Long_t i=0; i<numPix; i++) {
    distance += GetDistance(tgtPix[i],outPix[i]);
  }
  return distance;
}

Double_t GetDL(Long_t numPix, UInt_t* tgtPix, UInt_t* outPix, Long_t i, Long_t j) {
  Double_t di = GetDistance(tgtPix[i],outPix[i]) + GetDistance(tgtPix[j],outPix[j]);
  Double_t df = GetDistance(tgtPix[i],outPix[j]) + GetDistance(tgtPix[j],outPix[i]);
  return df-di;
}

Double_t GetT(Long_t numPix, UInt_t* tgtPix, UInt_t* outPix) {
  Double_t dL=0, dLMax=0;
  for (Long_t i=0; i<numPix; i++) {
    dL = GetDL(numPix,tgtPix,outPix,(Long_t)(gen->Rndm()*numPix),(Long_t)(gen->Rndm()*numPix));
    if (dL > dLMax) dLMax = dL;
  }
  return dLMax *= 100;
}

UInt_t Update(Long_t numPix, UInt_t* tgtPix, UInt_t* outPix, Double_t T) {
  Long_t i = (Long_t)(gen->Rndm()*numPix);
  Long_t j = (Long_t)(gen->Rndm()*numPix);
  Double_t dL = GetDL(numPix,tgtPix,outPix,i,j);
  if (dL <= 0 || gen->Rndm() < exp(-dL/T)) {
    UInt_t temp = outPix[i];
    outPix[i] = outPix[j];
    outPix[j] = temp;
    return 1;
  }
  return 0;
}

int main(int argc, char **argv){
  if (argc<3) {
    cout << "Usage: simapix_start image1 image2 <output=out.png>" << endl;
    return 0; 
  }
  TString fsrc = argv[1];
  TString ftgt = argv[2];
  TString fout;
  fout = (argc>3 ? argv[3] : "out.png");
  cout << "Reading images: source = " << fsrc << " target = " << ftgt << endl;
  cout << "Output = " << fout << endl;

  TApplication theApp("App", &argc, argv);

  // create image objects
  TASImage* src = new TASImage(fsrc.Data());
  TASImage* tgt = new TASImage(ftgt.Data());
  TASImage* out = new TASImage(*src); // start with copy of source

  // Test image geometry, exit if they are not the same dimensions
  assert ( src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight() );
  cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;
  Long_t numPix = src->GetWidth()*src->GetHeight();

  // *** The work happens here
  // access the pixels for the output image 
  // each pixel is a 32-bit word, 1 byte each for (alpha,red,green,blue)
  // don't touch alpha (bits 31:28)
  // UInt_t* srcPix = src->GetArgbArray();
  UInt_t* tgtPix = tgt->GetArgbArray();
  UInt_t* outPix = out->GetArgbArray();  

  // initial distance
  const Double_t L0 = GetTotalDistance(numPix,tgtPix,outPix);
  cout << "The initial color distance is: " << L0 << endl;
  // initial temperature
  Double_t T = GetT(numPix,tgtPix,outPix);
  cout << "The initial temperature is: " << T << endl;

  // Metropolis
  UInt_t success=0;
  // Double_t l,lMin=L0;
  // start
  // double sec;
  clock_t timer = clock();
  for (int i=0; i<100000; i++) {
    for (int j=0; j<100*numPix; j++) {
      success += Update(numPix,tgtPix,outPix,T);
      if (success==100) break;
    }
    success=0;
    T *= 0.9;
    
    // l = GetTotalDistance(numPix,tgtPix,outPix);
    // if (l < lMin) {
      // sec = (double)(clock()-timer)/CLOCKS_PER_SEC;
      // lMin = l;
      // printf("The color distance after optimizaton is: %-7f (i = %-6i T = %-7.3f time (s) = %-6.2f)\n", l,i,T,sec);
    // }
  }
  cout << (double)(clock()-timer)/CLOCKS_PER_SEC << " Fin." << endl;

  
  // // flip the image
  // for (int i=0;i< numPix/2; i++){
  //   unsigned pxl=outPix[i];
  //   outPix[i]=outPix[numPix-i-1];
  //   outPix[numPix-i-1]=pxl;
  // }

  // *************************


  // print the results
  TCanvas *c1 = new TCanvas("c1", "images", 640, 480);
  c1->Divide(2,2);

  c1->cd(1);
  c1->Draw();
  src->Draw("X");
  c1->cd(2);
  tgt->Draw("X");
  c1->cd(3);
  out->Draw("X");
  c1->Print("collage.png");
  
  // save the new image
  out->WriteImage(fout.Data());

  // coment out the lines for running in batch mode
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();

}
