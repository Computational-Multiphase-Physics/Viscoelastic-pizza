/** Title: Elastic Pizza
# Version: 2.0
# Updated: Nov 03, 2024

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# Main feature: 
- 2D+axi viscoelastic scalar implementation of the pizza problem.
- There is no surface tension (infinite Weber number).
*/

// 1 is Pizza and 2 is air
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseVE.h"
#include "navier-stokes/conserving.h"
#include "log-conform-viscoelastic-scalar-2D.h"
#include "tension.h" // uncomment to make Weber number finite

#define tsnap (1e-2)
#define tsnap2 (1e-4)

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

double tmax, Ldomain, Pi, tEtas, tLam, Ec, MuAirPizza, RhoR, H0;
int MAXlevel;

// top is outflow
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
// right is outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

// bottom forcing axi!
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));

/*
The charateristic scales: $R_0$, $G$, $\rho$: initial radius of the blob, elastic (pizza) modulus, and (pizza) density.

The characteristic frequency is $\Omega_c = \sqrt{G/(\rho R_0^2)}$: frequency required to significantly stretch the blob.

For Oldroyd-B, the parameters are: 
- Dimensionless retardation time: $t_{\eta s} = (\eta_s/G)*\Omega_c = \eta_s/\sqrt{\rho G R_0^2}$. $t_{\eta s}$ is the time it takes to show an elastic response
  - For a purely elastic solid, $t_{\eta s}$ is 0
  - $t_{\eta s} \to \infty$ is liquid. 

- Dimensionless relaxation time: $t_\lambda = \lambda * \Omega_c = \lambda * \sqrt{G/(\rho R_0^2)}$. $t_\lambda$ is the time it takes to show a viscous response
  - For a viscous liquid, $t_\lambda$ is 0 and for solids, $t_\lambda \to \infty$.
  - A purely elastic solid has $t_{\eta s} = 0$ and $t_\lambda \to \infty$.

### Additionally, there are two other dimensionless numbers:

#### Pizza number: $\Pi = \rho\Omega^2R_0^2/G$.

This is the ratio of inertial to elastic forces. This is the square of the ratio of two characteristic frequencies: $(\Omega/\Omega_c)^2$.<br>
$\Pi \gg 1$ is needed to show significant deformation of the blob. This is similar to the Bond number.

#### Elasto-capillary number: $Ec = \gamma/(G R_0)$.

This compares the elastic to capillary forces. To start with, we can assume $G \gg \gamma/R$.
- Another way to think about it is that $Ec = \Omega_\gamma^2/\Omega_c^2 = \frac{\gamma/(\rho R_0^3)}{G/(\rho R_0^2)} = \gamma/(G R_0)$.

## Oldroyd-B model implementation
For details of the constitutive model [Oldroyd-B model](https://en.wikipedia.org/wiki/Oldroyd-B_model), see: [https://github.com/comphy-lab/Viscoelastic3D](https://github.com/comphy-lab/Viscoelastic3D).
*/

char comm[80], restartFile[80], logFile[80];
int main(int argc, char const *argv[]) {
  
  MAXlevel = 9;
  tmax = 1e1; //atof(argv[3]);
  Ldomain = 4.0; //atof(argv[4]);

  /*
  Aspect ratio: H0/R. This is the initial height of the blob that will be stretched.
  */
  H0 = 1e0; // atof(argv[5]);

  /*
  Pizza number: \rho\Omega^2R^2/G. 
  */
  Pi = 4e0; // atof(argv[6]); 

  /* 
  Dimensionless retardation time: tEtas = (\eta_s/G)*\omega_c. For purely elastic solid, we need to keep tEas \to 0.  
  */
  tEtas = 5e-1; // atof(argv[7]); 
  /*
  Dimensionless relaxation time: tLam = \lambda * \omega_c. For purely elastic solid, we need to keep tLam \to \infty.
  */
  tLam = 1e30; // atof(argv[8]);   
  /*
  Elasto-capillary number: \gamma/(G R_0). For now, this value is only a placeholder. We will assume G \gg \gamma/R.
  */
  Ec = 1e-2; // atof(argv[9]); 

  MuAirPizza = 1e-2; // atof(argv[10]); // viscosity ratio \mu_{air}/\mu_{pizza dough}. This value needs to be close to 0. 
  RhoR = 1e-3; // atof(argv[11]); // density ratio \rho_{air}/\rho_{pizza dough}. This value needs to be very close to 0.


  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (6));

  /*
  In the folder called "intermediate", we will store the snapshot files as the simulation progresses. See event: writingFiles.
  */
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  /*
  We save the restart file every t = tsnap such that if the simulation is interrupted, we can recover from the last saved state.
  */
  sprintf (restartFile, "restart");
  /*
  We also save the log file every few time steps.
  */
  sprintf (logFile, "logFile.dat");


  /*
  Pizza */
  rho1 = 1.000; mu1 = tEtas;
  G1 = 1e0; 
  lambda1 = tLam; // infinite relaxation time: phase 1 is Kelvin--Voigt solid.
  /*
  Air */
  rho2 = RhoR; mu2 = MuAirPizza*mu1;
  G2 = 0.; lambda2 = 0.; // air is modeled as purely Newtonian. 
  /*
  Surface tension
  */
  f.sigma = Ec;

  // CFL and tolerance
  CFL = 1e-2;
  TOLERANCE = 1e-4;

  fprintf(ferr, "Level %d tmax %g, Pi %g, tEtas %g, tLam %g, Ec %g\n", MAXlevel, tmax, Pi, tEtas, tLam, Ec);
  run();
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y){
    if (y > 1e-20){
      double ff = (f[] + f[0,-1])/2.;
      av.y[] += ff*fm.y[]*Pi*(1.-RhoR)/rho(ff);
    }
  }
}

event init(t = 0){
  if(!restore (file = restartFile)){
    /**
    We can now initialize the volume fractions in the domain. 
    For sqrt(sq(x/H0) + sq(y)) < 1, f = 1.
    For sqrt(sq(x/H0) + sq(y)) > 1, f = 0.
    */
    // refine the interface region
    refine(sq(x/H0) + sq(y) < sq(1.05) && sq(x/H0) + sq(y) > sq(0.95) && level<MAXlevel);
    // define the volume fraction: VoF color function is continuous across the interface.
    fraction(f, 1e0-sqrt(sq(x/H0) + sq(y)));
  }
}

event adapt(i++) {
  /*
  We adapt the grid to resolve the interface and the flow. For infinite Weber number, these should be fine. 
  For finite Weber number, perhaps adaptation based on curvature is also needed. 
  */
 scalar KAPPA[];
 curvature (f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, A11, A12, A22, KAPPA},
    (double[]){fErr, VelErr, VelErr, AErr, AErr, AErr, KErr},
    MAXlevel, 4);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = restartFile);
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 0.5*(2*pi*y)*rho(f[])*(sq(u.y[])+sq(u.x[]))*sq(Delta);
  }

  scalar pos[];
  position (f, pos, {0,1,0});
  double rmax = statsf(pos).max;

  if (pid() == 0){
    static FILE * fp;
    if (i == 0) {
      fprintf (ferr, "i dt t ke rmax\n");
      fp = fopen (logFile, "w");
      fprintf(fp, "Level %d tmax %g, Pi %g, tEtas %g, tLam %g, Ec %g\n", MAXlevel, tmax, Pi, tEtas, tLam, Ec);
      fprintf (fp, "i dt t ke rmax\n");
    } else {
      fp = fopen (logFile, "a");
    }

    fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, rmax);
    fclose(fp);
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, rmax);
    
    if (ke < -1e-10){
      dump(file="FailedDumpKeNegative");
      fprintf(ferr, "Something very bad is happening! Stopping!\n");
      return 1;
    }
    if (ke < 1e-10 && i > 100){
      dump(file="StoppedDump");
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
    if (ke > 1e4){
      dump(file="FailedDump");
      fprintf(ferr, "kinetic energy blew up! Stopping!\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
    if (rmax > 0.95*Ldomain){
      dump(file="FailedDumpRmax");
      fprintf(ferr, "Rtip has approached very close to the domain boundary! Stopping!\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "Rtip has approached very close to the domain boundary! Stopping!\n");
      fclose(fp);
      return 1;
    }
  }

}

event end (t = tmax + tsnap) {
  dump(file="FinalDump");

  static FILE * fp;
  fp = fopen (logFile, "a");
  fprintf(fp, "Level %d tmax %g, Pi %g, tEtas %g, tLam %g, Ec %g\n", MAXlevel, tmax, Pi, tEtas, tLam, Ec);
  fprintf(fp, "End of simulation.\n");
  fclose(fp);
}
