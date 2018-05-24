#!/bin/bash
#BSUB -x
#BSUB -J test_mpi                          # Job name
#BSUB -n 4                                 # total number of MPI task
#BSUB -o test_mpi.%J.out                   # stdout filename
#BSUB -q computet1                         # queue name
#BSUB -a p8aff(1,8,1,balance)              # 1 thread/task, SMT=8, 1
#BSUB -R 'span[ptile=20]'                  # tile : number of MPI task/node
#BSUB -W 00:05                             # maximum runtime

# load some MPI module
module load gcc/4.8/ompi/2.1

# report bindings for cross-checking
mpirun --report-bindings ./test_mpi
