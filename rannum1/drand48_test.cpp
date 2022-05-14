/* 
************************************************************************
*  drand48_test.c: Program tests drand48() for nearest-neighbor        *
*  correlations.                                                       *
*								       *
************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
 
/* if you don't have drand48 uncomment the following two lines */
/*    #define drand48 1.0/RAND_MAX*rand
      #define srand48 srand                */

int main(int argc, char **argv) {
   double x1, x2;
   double sum, error;
   double nsigma;
   int nvals=100000;
   long seed=0;

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

   if (seed==0) srand48(time(0));
   else srand48(seed);                  // seed the number generator
   
   x1 = drand48();
   sum = 0;
   for(int n=0; n<nvals; n++){
     x2 = drand48();
     sum += x1*x2;
     x1 = x2; 
   }
   sum /= nvals;           // < x_i * x_i+1 >
   error = sum/0.25-1;    // relative error: < x_i * x_i+1 > / ( <x_i><x_i+1> ) - 1

   // using estimate sigma_mu = stdDeV/sqrt(N)
   double stdDev=sqrt(1./12.);
   double sigma=stdDev/sqrt(nvals/2.);  // expected error on mean
   nsigma = (sum-0.25)/sigma;

   printf("Generated %d pairs\n",nvals);
   printf("<x_i x_i+1> = %lf\t relative error = %lf\n", sum, error);
   printf("Result is roughly %f std deviations from the correct answer \n",nsigma);

   return(0);
}
