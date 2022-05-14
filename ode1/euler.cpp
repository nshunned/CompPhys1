#include <iostream>
#include <getopt.h>
#include <cmath>
using namespace std;


void usage(char **argv){
  fprintf(stderr, "\nUsage: %s [options]\n",argv[0]);
  fprintf(stderr, " -v vinit: initial speed [default 10 m/s]\n");
  fprintf(stderr, " -a theta0: initial angle above horizontal [45deg]\n");
  fprintf(stderr, " -t step: time step for approximation [.01s]\n");
  fprintf(stderr, " -h: print this message\n");
}

int main(int argc, char **argv) {
  double vinit = 10;  // m/s
  double theta0 = 45;  // deg
  double dt = 0.01;   // time step [s]
  int opt;
  while ((opt = getopt(argc, argv, "v:a:t:h")) != -1) {
    switch (opt) {
    case 'v':
      vinit = atof(optarg);
      break;
    case 'a':
      theta0 = atof(optarg);
      break;
    case 't':
      dt = atof(optarg);
      break;
    case 'h':
      usage(argv);
      return 0;
      break;
    default: /* '?' */
      ;
    }
  }

  printf("Simulating projectile motion with params:\n");
  printf("(vinit,theta0,dt)=(%7.2lf,%7.2lf,%7.2f)\n",vinit,theta0,dt);
  
  // fill in the blanks

  // given the boundary conditions at t=0, namely:
  // x(0) = y(0) = 0
  // vx(0) = vinit * cos(theta0) , vy(0) = vinit * sin(theta0)
  // define a method for stepping forward in time to calculate the
  // new position of the projectile in steps of dt
  // your simulation should stop at the last time step before y goes negative

  
  // write out a data file containing at least the following information
  // t   x(t)   y(t)   vx(t)   vy(t)  

  FILE *out;
  out = fopen("euler.dat", "w");
  fprintf(out, "t         x(t)      y(t)      vx(t)     vy(t)      K(t)      U(t)      E(t)\n");

  double g = -9.81;

  double t = 0; double x = 0; double y = 0;
  double vx = vinit*cos(theta0/360*2*M_PI);
  double vy = vinit*sin(theta0/360*2*M_PI);
  double K = 0.5*vinit*vinit;
  double U = g*y;
  double E = K+U;

  while(y >= 0) {
    fprintf(out, "%-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f\n", t, x, y, vx, vy, K, U, E);
    
    t += dt;
    x = x + vx*dt;
    y = y + vy*dt;
    vy = vy + g*dt;
    K = 0.5*(vx*vx + vy*vy);
    U = abs(g*y);
    E = K + U;
  }

  printf("Data stored in euler.dat\n");
  fclose(out);

  // writing out a data file containing the vollowing information which depends on varying dt
  // dt   x_f(dt)   y_f(dt)   x_f_err   y_f_err
  
  out = fopen("landing.dat", "w");
  fprintf(out, "dt        x_f(dt)   y_f(dt)   x_f_err\n");

  double x_f_act = vinit*vinit*sin(2*theta0/360*2*M_PI)/abs(g);
  for(int i = 0; i < 15; i++) {
    t = 0; x = 0; y = 0;
    vx = vinit*cos(theta0/360*2*M_PI);
    vy = vinit*sin(theta0/360*2*M_PI);

    double x_f, y_f;
    while(y >= 0) {
      t += dt;
      x_f = x;
      y_f = y;
      x = x + vx*dt;
      y = y + vy*dt;
      vy = vy + g*dt;
    }
    fprintf(out, "%-9.6f %-9.6f %-9.6f %-9.6f\n", dt, x_f, y_f, abs(x_f-x_f_act));
    dt /= 2;
  }
  
  printf("Data stored in landing.dat\n");
  fclose(out);

  return 0;
}
