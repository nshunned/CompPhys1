// n_volume.cpp : Calculate volume of an n-sphere
#include <stdio.h>  // Standard I/O headers, also in C++ can use <cstdio>
#include <math.h>   // std math library headers - define PI, " " <cmath>
 
int main() {
  double n, V;  // Double precision variables
  printf ("Enter the dimension of the sphere \n");  // Request input
  scanf("%lf", &n);  // Read from standard input
  V = pow(M_PI, n/2)/tgamma(n/2+1);
  // Calculate volume
  printf("%s :: dimension n = %lf, volume V = %lf\n", __FILE__, n, V);  // Print results

  /*
  for(int n = 0; n <= 20; n++) {
    double V = pow(M_PI, n/2)/tgamma(n/2+1);
    printf("n = %lf, V = %lf\n", (double)n, V);
  }
  */
} // End main
