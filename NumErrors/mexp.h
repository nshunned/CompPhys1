// mexp.h
#pragma once

#ifdef EXTERNC
extern "C" {
#endif
  double fact(double x);
  double mexp(double x);
  double mexp_0(double x);
  void usage(char* name);
#ifdef EXTERNC
}
#endif
