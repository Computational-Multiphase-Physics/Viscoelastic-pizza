# Viscoelastic-pizza

## Background:

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


## How to run the code:

One core:

```shell
qcc -O2 -Wall -disable-dimensions pizza.c -o pizza -lm 
./pizza
```

OpenMP parallelization:

```shell 
qcc -O2 -Wall -disable-dimensions pizza.c -o pizza -lm -fopenmp 
export OMP_NUM_THREADS=16
./pizza
```
Here, change 16 to the number of available threads. 

OpenMPI parallelization:

```shell
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions pizza.c -o pizza -lm
mpirun -np 16 ./pizza
```
Here, change 16 to the number of available cores. 


**Note:** In the [pizza.c](pizza.c) file, you can uncomment and edit ``argv" parts to pass parameters from terminal. 

# Postprocess: 

## Facets_Xjet.py
- This is a python script to create interfacial facets in raster (png) format.
- It takes the output from `./getFacet2D` and plots the facets using matplotlib.
- It also calculates rMax and vMax using the output from `./getRmaxNV` and plots (rMax) as a blue marker.

## VideoAxi.py
- This is a python script to create a video from the facets data.
- It takes the output from `./getFacet2D` and plots the facets in a video.
- It also gets the velocity gradient tensor norm and the elastic energy (normalized by G) and plots using imshow in matplotlib.

## gettingFacetsInPDF.ipynb
- This is a notebook to create interfacial facets in vector format. 
- It takes the output from `./getFacet2D` and plots the facets in a PDF.
- It also calculates rMax and vMax using the output from `./getRmaxNV`. You can also look at the logFile.dat that is created while running the simulation to get rMax(t) which will have a much better time frequency as compared to the output from `./getRmaxNV` (which is restricted to number of snapshots we same instead of saving at every time step as in logFile.dat).