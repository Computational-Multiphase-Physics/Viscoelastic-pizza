#!/bin/bash

cd 1000/
qcc -O2 -Wall -fopenmp -disable-dimensions pizza.c -o pizza -lm
export OMP_NUM_THREADS=8
./pizza 0.0
cd ../

cd 1001/
qcc -O2 -Wall -fopenmp -disable-dimensions pizza.c -o pizza -lm
export OMP_NUM_THREADS=8
./pizza 1e-2
cd ../

cd 1002/
qcc -O2 -Wall -fopenmp -disable-dimensions pizza.c -o pizza -lm
export OMP_NUM_THREADS=8
./pizza 1e-1
cd ../


cd 1003/
qcc -O2 -Wall -fopenmp -disable-dimensions pizza.c -o pizza -lm
export OMP_NUM_THREADS=8
./pizza 1e0
cd ../
