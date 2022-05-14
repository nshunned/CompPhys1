/* 
******************************************************************************
* onespin.c: Uses Metropolis algorithm to generate thermal ensemble for a    *
*	     single spin.                			             *
*  								             *
******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
 
// update state by flipping spin, accept new state if E_new<E_old
// or rand[0,1)< P(new)/P(old)
void update(int *s, double h) {
  int spin, newspin;
  double DeltaBetaE;
  spin = *s;
  newspin = -spin;                      // trial spin flip
  DeltaBetaE = -(newspin-spin)*h;       // beta*(E(new) - E(old))
  // if the new state is at lower energy, accept it
  // if not, accept it with a probability decreasing w/ exp(-DeltaBetaE)
  if ( DeltaBetaE <= 0 || drand48() < exp(-DeltaBetaE) ) *s = newspin; 
}

int main(int argc, char *argv[]) {
  int n, nsweep;
  int spin;
  double h;
  int nplus, nminus;

  printf("Program generates a thermal ensemble for states with one spin.\n\n");

  srand48((long)(time(0)));

  if (argc==3){
    nsweep=atoi(argv[1]);
    h=atof(argv[2]);
  }
  else{
    printf("Enter total number of spin configurations (sweeps) generated:\n");
    scanf("%d", &nsweep); 
    
    printf("Enter value of magnetic field parameter h=H/(k_bT):\n");
    scanf("%lf", &h);
  }

  // our initial state
  spin = 1;

  // Metropolis update loop
  nplus = nminus = 0;
  for(n=0; n<nsweep; n++) {
    update(&spin,h);           // try a Metropolis update
    if (spin == 1) nplus++;                        
    if (spin == -1) nminus++;
  }

  printf("Calculations\n");
  printf("P(+ state) = %lf\n", (double)nplus/(nplus+nminus));
  printf("P(- state) = %lf\n", (double)nminus/(nplus+nminus));

  // Write <magnetization> 
  printf("<sigma> = %lf\n", (double)(nplus-nminus)/(nplus+nminus));
  
  /* Write theoretical prediction for comparison  */
  printf("\nTheory prediction: tanh(h) = %lf\n", tanh(h));   
  
  return(0);
}
