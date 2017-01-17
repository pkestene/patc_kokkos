#!/bin/bash
#BSUB -x
#BSUB -J test_mpi_kokkos_openmp                    # Job name
#BSUB -n 4                                         # total number of MPI task
#BSUB -o test_mpi_kokkos_openmp.%J.out             # stdout filename
#BSUB -q compute                                   # queue name
#BSUB -R "affinity[core(10):cpubind=core]"         # affinity
#BSUB -R 'span[ptile=2]'                           # tile : number of MPI task/node
#BSUB -W 00:05


module load gcc/4.8 ompi/1.10

# number of OpenMP thread per MPI task
OMP_NUM_THREADS=20

EXE_NAME=test_mpi_kokkos.omp

# report bindings for cross-checking
mpirun --report-bindings ./$EXE_NAME
