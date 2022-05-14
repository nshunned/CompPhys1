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

// Generate pairs of random numners using drand48
// This example also shows the usage of getopt to do argument parsing
// The default behaviour is:
// to seed drand48 based on the clock:  -s 0
// generate 1000 pairs: -n 1000
// and create output files drand48_pairs.(dat/bin)

int main(int argc, char **argv) {
  double x1, x2;
  int nvals=100000;
  long seed=0;
  FILE *outfile, *outfile2;
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
    
  outfile = fopen("drand48_pairs.dat","w");
  outfile2 = fopen("drand48_pairs.bin","wb");

  printf("Program generates 2-dimensional distribution of ");
  printf("random numbers for plotting.\n\n");

  if (seed==0) srand48(time(0));
  else srand48(seed);                  // seed the number generator
  //printf("seed=%d\n",seed);

  for(int n=0; n<nvals; n++){
    x1 = drand48();
    x2 = drand48();
    fprintf(outfile, "%f\t%f\n", x1, x2);
    fwrite((void *)&x1,sizeof(double),1,outfile2);
    fwrite((void *)&x2,sizeof(double),1,outfile2);
  }
  printf("Output file for plotting is drand48_pairs.dat\n");

  fclose(outfile);
  fclose(outfile2);

  return (0);
}
