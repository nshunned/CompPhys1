/* 
************************************************************************
*  badrandom_test.c:                                                   * 
*  Program applies nearest neigbor test using badrandom		       *
*								       *
************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
long M;
double a, c;

double badrandom(double x)
{
   double rin, xout;
   rin = x*M;
   xout =  (a*rin + c)/M;
   xout -= (int) (xout);
   return xout;
}

int main()
{
   double x1, x2;
   double sum, error;
   double nsigma;
   int nmax;
   double xstart=0.5;

   printf("Program tests badrandom number generator for nearest-neighbor ");
   printf("correlations.\n\n");

   printf("Enter a, c, and M parameters for linear congruent random ");
   printf("number generator:\n");
   scanf("%lf %lf %li", &a, &c, &M);

   printf("Enter total number of randoms pairs:\n");
   scanf("%d", &nmax);

   x1 = xstart;
   sum = 0;
   for(int n=0; n<nmax; n++) {
     x2 = badrandom(x1);
     sum += x1*x2;
     x1 = x2; 
   }
   sum /= nmax;
   error = sum/0.25-1;

   // using estimate sigma_mu = stdDeV/sqrt(N)
   double stdDev=sqrt(1./12.);
   double sigma=stdDev/sqrt(nmax/2.);  // expected error on mean
   nsigma = (sum-0.25)/sigma;

   printf("<x_i x_i+1> = %lf\t relative error = %lf\n", sum, error);
   printf("Result is roughly %f std deviations from correct answer \n",nsigma);

   return(0);
}
