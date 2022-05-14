/* **********************************************************************
 * Integration example: compareinteg.cpp                                *
 * Compare results in integrators for different #'s of fcn calls        *
 *    			                                                *
 * Output files showing # of divisions and value of integral calculation*
 ************************************************************************
 */

#include <stdio.h>
#include <math.h>
#include "integrators.hpp"     

// note: leave XMIN=0.0 for now (b/c of a bug in the gaus integral code
const int MAXINTERVAL=500;             // max number of intervals
const double XMIN=0.0;                 // ranges of integration 
const double XMAX=1.0;

// the function we want to integrate
// note the use of a function pointer below to pass this function
// to the integration routines
double f (double x) { 
  return (cos(exp(100*x)));            
}
double true_value=0.6321205588285577L;

int main() {
  // make output files for the results from each integrator
  FILE *pTrap = fopen("trap.dat","w"); 
  FILE *pSimp = fopen("simp.dat","w");
  FILE *pGaus = fopen("gaus.dat","w");
  FILE *pGaus2 = fopen("gaus2.dat","w");

  double trap_last=0, simp_last=0, gaus_last=0;
  double gaus2_last=0, gaus2_new=0;
  GaussInt gaussint2;

  double trap_new=0, simp_new=0, gaus_new=0;
  double trap_d=0, simp_d=0, gaus_d=0, gaus2_d=0;  // difference between iterations 
  double trap_e=0, simp_e=0, gaus_e=0, gaus2_e=0;  // actual error
  
  int nInvl=2;   // start with 2 intervals
  while (nInvl<=MAXINTERVAL){
    
    trap_new = trapez(f, nInvl, XMIN, XMAX);
    trap_d=fabs(trap_new-trap_last);
    trap_e=fabs(trap_new-true_value);    
    trap_last=trap_new;
    if (nInvl>2)
      fprintf(pTrap, "%i\t%.15lf\t%.15lf\t%.15lf\n", nInvl, trap_d, trap_e, trap_new);
    
    if (nInvl%2) {// simpson's rule needs odd # of intervals
      simp_new = simpson(f, nInvl, XMIN, XMAX);
      simp_d=fabs(simp_new-simp_last);
      simp_e=fabs(simp_new-true_value);
      simp_last=simp_new;
      if (nInvl>3)
	fprintf(pSimp, "%i\t%.15lf\t%.15lf\t%.15lf\n", nInvl, simp_d, simp_e, simp_new);
    }
    if (nInvl<50) {
      gaus_new = gaussint(f, nInvl, XMIN, XMAX);
      gaussint2.Init(nInvl);
      gaus2_new=gaussint2.Integ(f,XMIN,XMAX);
      gaus_d=fabs(gaus_new-gaus_last);
      gaus_e=fabs(gaus_new-true_value);
      gaus2_d=fabs(gaus2_new-gaus2_last);
      gaus2_e=fabs(gaus2_new-true_value);
      gaus_last=gaus_new;
      gaus2_last=gaus2_new;
      if (nInvl>2){
	fprintf(pGaus, "%i\t%.15lf\t%.15lf\t%.15lf\n", nInvl, gaus_d, gaus_e, gaus_new);
	fprintf(pGaus2, "%i\t%.15lf\t%.15lf\t%.15lf\n", nInvl, gaus2_d, gaus2_e, gaus2_new);
      }
    }

    nInvl++;
  }

  printf("data stored in compareinteg.dat\n");
  // n.b. "true" is only true w/in our limited precision
  printf("result, estimated and true error:\n");
  printf("trap rule: %.15lf\t%.15lf\t%.15lf\n",trap_new,trap_d,trap_e);
  printf("simp rule: %.15lf\t%.15lf\t%.15lf\n",simp_new,simp_d,simp_e);
  printf("gaus rule: %.15lf\t%.15lf\t%.15lf\n",gaus_new,gaus_d,gaus_e);
  printf("gaus2 rule: %.15lf\t%.15lf\t%.15lf\n",gaus2_new,gaus2_d,gaus2_e);
  fclose(pGaus2);
  fclose(pTrap);
  fclose(pSimp);
  fclose(pGaus);
  return 0;
}

