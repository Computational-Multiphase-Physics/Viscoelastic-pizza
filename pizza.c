/** Title: Encapsulation
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Jan 07, 2021
*/

// 1 is Si Pool, 2 is Water Drop and 3 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "log-conform-Elastic_v9.h"
// #include "tension.h"

#define tsnap (1e-1)
#define tsnap2 (1e-4)
#define R0 (1e0)
#define h0 (1e0)

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)

double tmax, Ldomain, Omega, El, Re, We, MuAirPizza, RhoR;
int MAXlevel;
scalar Gpd[];

// u.n[top] = neumann(0.);
// p[top] = dirichlet(0.);

int main(int argc, char const *argv[]) {
  
  MAXlevel = 10;
  tmax = 4e0; //atof(argv[3]);
  Ldomain = 6.0; //atof(argv[4]);

  Omega = 1e0; // it is the repeating variable

  El = 1e-1; // G/\rho\Omega^2R^2
  Re = 1e1; // \rho\OmegaR^2/\mu
  We = 1e30; // \rho\Omega^2R^3/\gamma
  MuAirPizza = 1e-3; // viscosity ratio \mu_{air}/\mu_{pizza dough}
  RhoR = 1e-2; // density ratio \rho_{air}/\rho_{pizza dough}


  L0=Ldomain;
  X0=0.0; Y0=0.;
  init_grid (1 << (6));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.000; mu1 = 1e0/Re;
  rho2 = RhoR; mu2 = MuAirPizza/Re;
  // polymers
  Gp = Gpd;

  // f.sigma = 1e0/We;

  fprintf(ferr, "Level %d tmax %g, El %g, Re %g, De Infty, We Infty\n", MAXlevel, tmax, El, Re);
  run();
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y){
    double ff = (f[] + f[0,-1])/2.;
    if (y > 0.)
      av.y[] += ff*fm.y[]*sq(Omega);
  }
}

event properties (i++) {
  foreach () {
   Gpd[] = El*clamp(f[], 0., 1.); // (f[] > 1.-1e-6 ? El: 0.); //// this is an artificial patch for now. The code has issues with VE terms in the interfacial cells!
  //  lambdad[] = De*clamp(f[], 0., 1.); //f[] > 1.-1e-6 ? De: 0.); //De*clamp(f[], 0., 1.);
  }
}

event init(t = 0){
  if(!restore (file = "dump")){
    /**
    We can now initialize the volume fractions in the domain. */
    refine(sq(x)/h0 + sq(y) < 1.05*R0 && sq(x)/h0 + sq(y) > 0.95*R0 && level<MAXlevel);
    fraction(f, R0-(sq(x)/h0 + sq(y)));
    // refine(x<(h0/2.0)+0.025 && y < 1e0+(h0/2.0)+0.025 && level<MAXlevel);
    // fraction(f, y < 1e0+(h0/2.0) ? sq(h0/2.0)-(sq(x)+sq(y-(h0/2.0)-1e0)) : (h0/2.0)-x);
  }
  // mask (x > 10.0 ? right : none);
}

event adapt(i++) {
  adapt_wavelet ((scalar *){f, u.x, u.y},
    (double[]){fErr, VelErr, VelErr},
    MAXlevel, MAXlevel-4);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.y[])+sq(u.x[]))*f[];
  }

  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf(fp, "Level %d tmax %g, El %g, Re %g, De Infty, We Infty\n", MAXlevel, tmax, El, Re);
    fprintf (fp, "i dt t ke\n");
  } else {
    fp = fopen ("log", "a");
  }

  fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
  fclose(fp);
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  assert(ke > -1e-10);
  // dump (file = "dump");
  if (ke < 1e-6 && i > 100){
    fprintf(ferr, "kinetic energy too small now! Stopping!\n");
    fp = fopen ("log", "a");
    fprintf(fp, "kinetic energy too small now! Stopping!\n");
    fclose(fp);
    return 1;
  }

}
