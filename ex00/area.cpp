// area.cpp : Calculate area of a circle 
#include <stdio.h>  // Standard I/O headers, also in C++ can use <cstdio>
#include <math.h>   // std math library headers - define PI, " " <cmath>
 
int main() {
  double rad, A;  // Double precision variables
  printf ("Enter the radius of a circle \n");  // Request input
  scanf("%lf", &rad);  // Read from standard input
  A = rad * rad * M_PI;
  // Calculate area
  printf("%s :: radius r= %lf, area A = %lf\n", __FILE__, rad, A);  // Print results
} // End main
