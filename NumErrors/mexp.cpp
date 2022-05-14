// mexp.cpp

#include "mexp.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

// not a good solution, but try running it for x<1 to begin

const int MAXIT=1000;      // max iterations for series
const double EPS=1e-9;     // stop for corrections smaller than this

#ifdef EXTERNC
extern "C" {
#endif
    
  double fact(double x){
    if (x < 1)
      return 1;
    return x*fact(x-1);
  }

  double mexp(double x) {
    if (x > 1)
      return mexp(x/2)*mexp(x/2);
    return mexp_0(x);
  }

  double mexp_0(double x) {
    double prev = 1.0;
    double value = 1.0;
    int sign;
    double term;

    for (int n = 1; n < MAXIT; n++) {
      n%2 ? sign = -1 : sign = 1;
      term = sign*prev*x/n;
      prev = fabs(term);
      value += term;
      
      if (fabs(term/value) < EPS || !isfinite(value))
	break;
    }
    return value;
  }

  void usage(char* name) {
    printf("Usage %s <positive value>\n",name);
  }
  
#ifdef EXTERNC
}
#endif

int main(int argc, char *argv[]) {
  if (argc < 2) {
    usage(argv[0]);
    return -1;
  }
  
  double x = atof(argv[1]);
  if (x < 0) {
    usage(argv[0]);
    return -1;
  }

  double val = mexp(x);
  double relerr = fabs(val-exp(-x))/exp(-x);

  printf("exp(-%lf) = %lf, rel error = %lg\n",x,val,relerr);

  return 0;
}
