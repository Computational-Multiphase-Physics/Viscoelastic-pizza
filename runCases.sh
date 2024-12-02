#!/bin/bash

qcc -O2 -Wall -fopenmp -disable-dimensions pizza.c -o pizza -lm
export OMP_NUM_THREADS=8
./pizza 10