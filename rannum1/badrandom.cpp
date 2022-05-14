/* 
**************************************************************************
*  badrandom_plot.c:                                                     *
*        Makes a 2D plot of random numbers for visual inspection         *
*        using a linear congruence random number generator.              *
*                                                                        *
* I use the name "badrandom" below because it could be a bad random      *
* number generator depending on your choices of a, c, and M (and         *
* even for good choices could be bad for very serious applications       *
* needing lots of random numbers.                                        *
*                                                                        *
*                                                                        *
**************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
long M;
double a, c;

/* I use the name "badrandom" below because it could be a bad random
** number generator depending on your choices of a, c, and M (and
** even for good choices could be bad for very serious applications
** needing lots of random numbers.
*/

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
   double xstart=0.;
   double x1, x2;
   int n, nmax;
   FILE *outfile, *outfile2;

   outfile = fopen("badrandom_pairs.dat","w");
   outfile2 = fopen("badrandom_pairs.bin","wb");

   printf("Program generates 2-dimensional distribution of random numbers ");
   printf("for plotting.\n");

   printf("Enter a, c, and M parameters for linear congruent random ");
   printf("number generator:\n");
   scanf("%lf %lf %li", &a, &c, &M);

   printf("Enter total number of random pairs:\n");
   scanf("%d", &nmax);
   x1 = xstart;

   for(n=0; n<nmax; n++)
   {
     x2 = badrandom(x1);
     x1 = badrandom(x2);;
     fprintf(outfile, "%lf\t%lf\n", x1, x2);
     fwrite((void *)&x1,sizeof(double),1,outfile2);
     fwrite((void *)&x2,sizeof(double),1,outfile2);
   }
   printf("Output for plotting in file badrandom_pairs.dat\n");

   fclose(outfile);
   fclose(outfile2);

   return(0);
}

