/* 
******************************************************************************
* ising1d.c                                                                  *
* =========                                                                  *
*    Uses Metropolis algorithm to generate thermal ensemble for the          *
*    one-dimensional ising chain with free boundary conditions               *
*    in a magnetic field h.                                                  *
*								             *
*    The state of the lattice is stored as spins +-1 in elements 1     *
*    through N of an array of length N+2.  The free boundary conditions      *
*    are handled by fixing elements 0 and N+1 of the array to be zero.       *
*  								             *
******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/* update_spin (s, env)
** -----------
** Do a metropolis update on a spin *s whose "environment" is env.
** The environment should be set beforehand so that the dependence
** of the total beta*energy on the value of the selected spin is
** given by (selected spin)*env 
*/


void update_spin(int *s, float env) {
  int spin = *s;
  int newspin = ( drand48() < 0.5 ? 1 : -1 );
  double DeltaBetaE = -(newspin-spin)*env;       /* beta*(E(new) - E(old)) */
  if ( DeltaBetaE <= 0 || drand48() < exp(-DeltaBetaE) ) *s = newspin; 
}

/* sweep (spin, N, beta, h)
** -----
** Sweep once through all the sites of the lattice spin of length N,
** trying an update at each site with inverse temperature beta and external
** magnetic field parameter h.
*/
void sweep(int *spin, int N, float beta, float h) {
  for(int ns=1; ns<=N; ns++) {
    update_spin(&(spin[ns]), beta*(spin[ns-1] + spin[ns+1]) + h);
  }
}

/* InitializeAllSpinUp (spin, N)
** =============================
**   Initialize all the N spins of the lattice "spin" to be +1.  Also
** initialize the two fake "boundary" spins to zero to implement free
** boundary conditions.
*/

void InitializeAllSpinUp(int *spin, int N) {
  for(int ns=1; ns<=N; ns++) spin[ns] = 1;
  spin[0] = spin[N+1] = 0;             // no interaction w/ endpoints
}


int main(){
  int ntherm=10000;           // initial sweeps to "thermalize" the system
  const int VisualDisplay=1;  // 1 or 0 to turn on/off display of chains 

  int nsweep, N, nmag=0, ntotal=0;
  int *spin;
  float beta, h;
  double expected;


  printf("Program generates a thermal ensemble for chain of N spins\n");
  printf("with free boundary conditions.\n\n");

  srand48((long)time(0));               // seed the number generator

  printf("Enter number of spins in chain:\n");
  scanf("%d", &N); 

  printf("Enter total number of configurations (sweeps) generated:\n");
  scanf("%d", &nsweep); 

  printf("Enter value of magnetic field parameter h:\n");
  scanf("%f", &h);

  printf("Enter temperature parameter beta (= 1/kT):\n");
  scanf("%f", &beta);

  // make room for N spins and two fake boundary spins
  spin = (int*)malloc((N+2)*sizeof(int)); 
  InitializeAllSpinUp(spin, N);  // cold start

  // sweep ntherm times to thermalize system
  for(int i=0; i<ntherm; i++) sweep(spin,N,beta,h);

  // now sweep through lattice nsweep times
  // we will average the magnitization of the system
  // over nsweeps to get the average "thermalized" value
  nmag = ntotal = 0;
  for(int n=0; n<nsweep; n++) {
    sweep(spin,N,beta,h);
    for(int ns=1; ns<=N; ns++) {
      nmag += spin[ns];
      ntotal++;
      if (VisualDisplay) {
        if(spin[ns]==1) printf("+");
        if(spin[ns]==-1) printf("-");
      }
    }
    if (VisualDisplay) printf("\n");
   }

   printf("Magnetization of thermalized system: <s> = %f\n", (float)nmag/ntotal);
   expected = exp(beta)*sinh(h) / 
     sqrt( exp(2*beta)*sinh(h)*sinh(h) + exp(-2*beta) );
   printf("Expected result for infinite length chain over long time: %f\n",
	  expected);

   return(0);
}
