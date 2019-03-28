#!/bin/bash
#BSUB -x
#BSUB -n 4                                # total number of MPI tasks
#BSUB -J laplace2d_kokkos                 # Job name
#BSUB -o laplace2d_kokkos.%J.out          # stdout filename
#BSUB -q computet1                        # queue name
#BSUB -a p8aff(10,8,1,balance)            # threads/task, SMT mode, nb logical cores per thread
#BSUB -R "affinity[core(1):cpubind=core]" # affinity
#BSUB -W 00:05                            # maximum runtime

module load at/10.0
module load gcc/4.8/ompi/2.1
module load cuda/9.2

OMP_DISPLAY_ENV=true

# for MPI + OpenMP on Power8
EXE_NAME=laplace2d_kokkos_mpi.omp
mpirun --report-bindings -n ${LSB_DJOB_NUMPROC} ./$EXE_NAME 

