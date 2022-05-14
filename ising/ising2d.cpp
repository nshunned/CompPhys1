/* 
******************************************************************************
* ising2d.c                                                                  *
* =========                                                                  *
*    Uses Metropolis algorithm to generate thermal ensemble for the          *
*    two-dimensional ising lattice with free boundary conditions             *
*    in a magnetic field h.                                                  *
*								             *
* Storage: The state of the lattice is stored as spins +-1 in elements       *
*    [1..N][1..N] of a 2-d array of size (N+2) by (N+2).  The free boundary  *
*    conditions are handled by fixing elements [0][j], [N+1][j], [i][0],     *
*    and [i][N+1] of the array to be zero.                                   *
*                                                                            *
* Compile:  cc -o ising2d ising2d.c -lm -lcurses                             *
*  								             *
******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>      /* for sleep system call */
#include <time.h>

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

int ntherm = 500;
const int VisualDisplay = 1;
const int SleepTime = 100000;   // in microseconds

/*
**  THE LATTICE
*/

int spin[NX+2][NY+2];       // put fake spins around edge of lattice 


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

  for(nx=1; nx<=NX; nx++) for(ny=1; ny<=NY; ny++) {
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

void InitializeHot() {
  printf("Initializing system\n");
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
**   and the program pauses for SleepTime useconds after each display.
*/

void DisplayLattice(int sweep_number) {
  int nx, ny;

  if (SleepTime>0)  puts( "\033[2J" );
  for (nx=1; nx<=NX; nx++) {
    for (ny=1; ny<=NY; ny++) {
      if (spin[nx][ny]== 1) printf("X");
      if (spin[nx][ny]==-1) printf("-");
    }
    printf("\n");
  }
  if (sweep_number>=0)
    printf("sweep %i:   magnetization <sigma> = %f\n", 
	   sweep_number, Magnetization());
  if (SleepTime>0) usleep(SleepTime);
  else printf("\n");
}


/* main
** ====
*/

int main()
{
  int n, nsweep, nx, ny;
  long ntotal, nmag;
  float beta, h;

  printf("Program generates a thermal ensemble a 2D Ising model of");
  printf("%ix%i spins with free boundary conditions.\n\n",
	 NX,NY);

  srand48((long)time(0));               // seed the number generator

  
  printf("Enter total number of configurations generated:\n");
  scanf("%d", &nsweep); 

  printf("Enter value of magnetic field parameter h:\n");
  scanf("%f", &h);

  printf("Enter temperature parameter beta (= 1/kT):\n");
  scanf("%f", &beta);

  InitializeHot();

  /* sweep ntherm times to thermalize system */
  printf("Thermalizing system, %d sweeps\n",ntherm);
  for(n=0; n<ntherm; n++) {
    sweep(beta,h);
    if (VisualDisplay &&  n%(ntherm/10)==0 ) {
      DisplayLattice(-1);
      printf("Thermalization sweep %d\n",n);
      usleep(1000000);
    }
  }

  /* now sweep through lattice nsweep times */
  nmag=ntotal=0;
  for(n=0; n<nsweep; n++) {
    if (VisualDisplay) DisplayLattice(n);
    sweep(beta,h);
    for(nx=1; nx<=NX; nx++) for(ny=1; ny<=NY; ny++) {
      nmag += spin[nx][ny];
      ntotal++;
    }
  }

  if (VisualDisplay) DisplayLattice(nsweep);
  printf("Average Magnetization: <s> = %f\n", (float)nmag/ntotal);

  return(0);
}
