/* 
******************************************************************************
* ising2d_vs_T.c                                                             *
* =========                                                                  *
*    This program plots the result of the average spin mangentization        *
* vs. temperature T for the 2D Ising model with free boundary conditions.    *
*								             *
* Storage: The state of the lattice is stored as spins +-1 in elements       *
*    [1..N][1..N] of a 2-d array of size (N+2) by (N+2).  The free boundary  *
*    conditions are handled by fixing elements [0][j], [N+1][j], [i][0],     *
*    and [i][N+1] of the array to be zero.                                   *
*                                                                            *
* Compile:  cc -o ising2d_vs_T ising2d_vs_T.c -lm -lcurses                   *
*  								             *
******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>    /* for sleep system call */
#include <time.h>


/* if you don't have drand48 uncomment the following two lines */
/*    #define drand48 1.0/RAND_MAX*rand
      #define srand48 srand                */


/*
** PARAMETERS
** ----------
**   NX, NY:        The lattice is NX by NY sites.
**   ntherm:        Number of "thermalization" sweeps to do before starting.
**   VisualDisplay: 1 or 0 to turn on/off display of spins.
**   SleepTime:     seconds to pause between displays of each sweep.
*/

#define NX 64
#define NY 64

int ntherm = 1000;
const int VisualDisplay = 1;
const int SleepTime = 300000;  // in microseconds

/*
**  THE LATTICE
*/

int spin[NX+2][NY+2];       /* put fake spins around edge of lattice */


/* update_spin (s, env)
** -----------
**   Do a metropolis update on a spin *s whose "environment" in env.
** The environment should be set beforehand so that the dependence
** of the total beta*energy on the value of the selected spin is
** given by
**               (selected spin)*env .
**
** sweep (spin, N, beta, h)
** -----
**   Sweep once through all the sites of the lattice spin of length N,
** trying an update at each site with inverse temperature beta and external
** magnetic field parameter h.
*/

void update_spin(int *s, float env)
{
  int spin, newspin;
  double DeltaBetaE;

  spin = *s;
  newspin = ( drand48() < 0.5 ? 1 : -1 );
  DeltaBetaE = -(newspin-spin)*env;       /* beta*(E(new) - E(old)) */
  if ( DeltaBetaE <= 0 || drand48() < exp(-DeltaBetaE) ) *s = newspin; 
}

void sweep(float beta, float h)
{
  int nx,ny;
  float environment;

  for(nx=1; nx<=NX; nx++) 
    for(ny=1; ny<=NY; ny++) {
      environment = 
        beta*(spin[nx][ny-1] + spin[nx][ny+1] + spin[nx-1][ny] + spin[nx+1][ny])
        + h;
      update_spin(&(spin[nx][ny]), environment);
  }
}

/* InitializeHot (spin)
** ====================
**   Initialize all the NX by NY spins of the lattice randomly.  Also
** initialize the two fake "boundary" spins to zero to implement free
** boundary conditions.
*/

void InitializeHot()
{
  int nx, ny;

  for(nx=0; nx<=NX+1; nx++) for(ny=0; ny<=NY+1; ny++) {
    if (nx==0 || nx==NX+1 || ny==0 || ny==NY+1) spin[nx][ny]=0;
    else spin[nx][ny] = (drand48() < 0.5 ? 1 : -1);
  }
}

/* Magnetization
** =============
**   Return the volume average of the spin for the current spin
** configuration.
*/

float Magnetization() {
  int nx,ny,nmag;

  nmag = 0;
  for(nx=1; nx<=NX; nx++) for(ny=1; ny<=NY; ny++) nmag += spin[nx][ny];
  return( (float)nmag/(NX*NY));
}


/* DisplayLattice
** ==============
**   Print out the 2D lattice on the screen.
**   If SleepTime is > 0, then the screen is cleared before each display,
** and the program pauses for SleepTime useconds after each display.
*/

void DisplayLattice(float T) {
  int nx, ny;

  if (SleepTime>0)  puts( "\033[2J" );
  for (nx=1; nx<=NX; nx++) {
    for (ny=1; ny<=NY; ny++) {
      if (spin[nx][ny]== 1) printf("X");
      if (spin[nx][ny]==-1) printf("-");
    }
    printf("\n");
  }
  printf("T = %f:   magnetization <sigma> = %f\n", T, Magnetization());
  if (SleepTime>0) usleep(SleepTime);
  else printf("\n");
}




int main() {
  int n, nsweep, nx, ny, itemp, ntemp;
  long ntotal, nmag, nenergy, nenergy2;
  float beta, h, Tmax, T;
  FILE *output;
  const char *OutputFileName = "ising2d_vs_T.dat";

  output = fopen(OutputFileName,"w");

  printf("Program calculate <sigma> vs. T for a 2D Ising model of");
  printf("%ix%i spins with free boundary conditions.\n\n",NX,NY);

  srand48((long)time(0));		      // seed the generator 

  
  printf("Enter # sweeps per temperature sample:\n");
  scanf("%d", &nsweep); 

  printf("Enter value of magnetic field parameter h:\n");
  scanf("%f", &h);

  printf("Enter starting value (maximum) of temperature T (=1/beta):\n");
  scanf("%f", &Tmax);

  printf("Enter # temperatures to simulate:\n");
  scanf("%i", &ntemp);

  InitializeHot();

  /* now do ntemp temperatures between Tmax and 0 */
  for (itemp=ntemp; itemp>0; itemp--) {
    T = (Tmax*itemp)/ntemp;
    beta = 1/T;

    /* sweep ntherm times to thermalize system at new temperature */
    for(n=0; n<ntherm; n++) sweep(beta,h);

    /* then sweep through lattice nsweep times for that temperature */
    nmag=ntotal=nenergy=nenergy2=0;
    for(n=0; n<nsweep; n++) {
      sweep(beta,h);
      for(nx=1; nx<=NX; nx++) for(ny=1; ny<=NY; ny++) {
	      nmag += spin[nx][ny];
        nenergy += (spin[nx][ny]*spin[nx+1][ny] + spin[nx][ny]*spin[nx][ny+1]);
	      ntotal++;
      }
      nenergy2+=nenergy*nenergy;
    }
    nenergy *= -1;
    float avg_nenergy = (float)nenergy/nsweep;
    float avg_nenergy2 = (float)nenergy2/nsweep;
    float capacity = (avg_nenergy2 - avg_nenergy*avg_nenergy)/(NX*NX*NY*NY)*beta*beta;
    fprintf(output, "%f %f %f %f\n", T, (float)nmag/ntotal, avg_nenergy, capacity);
    if (VisualDisplay) DisplayLattice(T);
  }

  printf("Output file is %s\n",OutputFileName);

  return(0);
}
