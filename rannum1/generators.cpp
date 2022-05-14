// This code may be compiled to make a stand alone exe
// or it can be run from the ROOT command line as:
// .L generators.cpp+
// GeneratorsTest()

#include "TApplication.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom2.h"
#include <cstdlib>
#include <iostream>

using namespace std;

// simple example to use different types of random number generators
// to generate numbers in range [0:1)
void GeneratorsTest(){
  TRandom tr(0);    // seed = 0 is default
  TRandom2 tr2(0);  // and will generate a seed based on the time
  int seed=tr.GetSeed();
  srand(seed);      // grab a seed and us it to set other generators
  srand48(seed);
  cout << "rand()    ";
  for (int i=0;i<5;i++) cout << 1.0* rand() / RAND_MAX << " ";
  cout << endl;
  cout << "drand48() ";
  for (int i=0;i<5;i++) cout << drand48() << " ";
  cout << endl;
  cout << "TRandom   ";
  for (int i=0;i<5;i++) cout << tr.Uniform() << " ";
  cout << endl;
  cout << "TRandom2  ";
  for (int i=0;i<5;i++) cout << tr2.Uniform() << " ";
  cout << endl;
}

int main(int argc, char **argv) {
  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this can be ignored
  TApplication theApp("App", &argc, argv);
  GeneratorsTest();

  // include these two lines if you want to view graphics in ROOT
  // cout << "To exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program
  // theApp.Run(true);
  return 0;
}
