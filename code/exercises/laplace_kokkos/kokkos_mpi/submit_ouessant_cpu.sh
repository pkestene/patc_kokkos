#!/bin/bash
#BSUB -x
#BSUB -J laplace2d_kokkos                 # Job name
#BSUB -n 10                               # total number of MPI tasks
#BSUB -o laplace2d_kokkos.%J.out          # stdout filename
#BSUB -q compute                          # queue name
#BSUB -R "affinity[core(1):cpubind=core]" # affinity
#BSUB -W 00:05                            # maximum runtime

module load at/10.0
module load gcc/4.8/ompi/2.1
module load cuda/9.0

OMP_DISPLAY_ENV=true

# for MPI + OpenMP on Power8
EXE_NAME=laplace2d_kokkos.omp
mpirun -n ${LSB_DJOB_NUMPROC} ./$EXE_NAME 

