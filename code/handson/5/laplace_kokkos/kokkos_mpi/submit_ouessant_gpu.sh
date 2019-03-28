#!/bin/bash
#BSUB -x
#BSUB -n 4                                         # total number of MPI task
#BSUB -gpu "num=4:mode=exclusive_process:mps=no:j_exclusive=yes"
#BSUB -J laplace_kokkos                            # Job name
#BSUB -o laplace_kokkos.%J.out                     # stdout filename
#BSUB -q computet1                                 # queue name
#BSUB -a p8aff(5,8,1,balance)                      # 5 threads/task, so that only 2 tasks/CPU
#BSUB -R 'span[ptile=4]'                           # tile : number of MPI task/node (1 MPI task <--> 1 GPU)
#BSUB -W 00:05                                     # maximum runtime

module load at/10.0 gcc/4.8/ompi/1.10 cuda/9.2

OMP_DISPLAY_ENV=true
NUMBER_OF_GPUS_PER_NODES=4

# for MPI + CUDA on Pascal P100
# option --ndevices tells Kokkos to activate all GPU on node seen by hwloc
EXE_NAME=laplace2d_kokkos_mpi.cuda
mpirun --report-bindings ./$EXE_NAME --ndevices=$NUMBER_OF_GPUS_PER_NODES
