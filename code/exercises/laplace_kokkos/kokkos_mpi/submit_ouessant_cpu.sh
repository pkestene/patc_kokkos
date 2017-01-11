#!/bin/bash
#BSUB -x
#BSUB -J poisson_kokkos                   # Job name
#BSUB -n 10                               # total number of MPI tasks
#BSUB -o poisson_kokkos.%J.out            # stdout filename
#BSUB -q compute                          # queue name
#BSUB -R "affinity[core(1):cpubind=core]" # affinity
#BSUB -W 00:05                            # maximum runtime

module load gcc/4.8/ompi/1.10 
moudle load cuda/8.0
module load kokkos/openmp_gnu485_dev

OMP_DISPLAY_ENV=true

# for MPI + OpenMP on Power8
EXE_NAME=laplace2d_kokkos.omp
mpirun -n ${LSB_DJOB_NUMPROC} ./$EXE_NAME 

