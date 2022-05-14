/* 
******************************************************************************
*  drand68_plot.c:  Makes a 2D plot of random numbers for visual inspection  *
*                   using drand48().                                         *
*  								             *
******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "TRandom.h"
#include "TRandom2.h"

int main(int argc, char **argv) {
  double x1, x2, y1, y2, z1, z2;
  int nvals=100000;
  long seed=0;
  FILE *outfile, *outfile2, *outfile3;
  int c;
  while ((c = getopt (argc, argv, "hs:n:")) != -1){
    switch (c){
    case 's':
      seed=atoi(optarg);
      break;
    case 'n':
      nvals=atoi(optarg);
      break;
    case 'h':
      printf("usage: drand48 -s [seed(0)] -n [number of pairs(1000)] \n");
      return(0);
    }
  }
  printf("seed=%d\n",seed);
    
  outfile = fopen("rand_100k.dat","w");
  outfile2 = fopen("TRandom1_100k.dat","w");
  outfile3 = fopen("TRandom2_100k.dat","w");

  TRandom* r = new TRandom();
  TRandom1* r1 = new TRandom1();

  for(int n=0; n<nvals; n++){
    x1 = rand();
    x2 = rand();
    fprintf(outfile, "%f\t%f\n", x1, x2);
    /*
    y1 = r->Rndm();
    y2 = r->Rndm();
    fprintf(outfile1, "%f\t%f\n", y1, y2);

    z1 = r1->Rndm();
    z2 = r1->Rndm();
    fprintf(outfile2, "%f\t%f\n", z1, z2);
    */
  }

  fclose(outfile);
  fclose(outfile2);
  fclose(outfile3);

  return (0);
}
