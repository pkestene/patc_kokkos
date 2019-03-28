#!/bin/bash
#BSUB -x
#BSUB -n 8                                          # number of MPI tasks
#BSUB -gpu "num=4:mode=exclusive_process:mps=no:j_exclusive=yes"
#BSUB -J test_mpi_kokkos_cuda                       # Job name
#BSUB -o test_mpi_kokkos_cuda.%J.out                # stdout filename
#BSUB -q computet1                                  # queue name
#BSUB -a p8aff(5,1,1,balance)                       # 5 threads/task, so that only 2 tasks/CPU, 1 task <-> 1 GPU
#BSUB -R 'span[ptile=4]'                            # tile : number of MPI task/node (1 MPI task <--> 1 GPU)
#BSUB -W 00:05


module load at/10.0 gcc/4.8/ompi/2.1 cuda/9.2

# This variable is normally set by the job scheduler
# As of January, 9th 2017, we enforce its value here
# to make sure all GPU devices can be used by our job.
#CUDA_VISIBLE_DEVICES=0,1,2,3

NUMBER_OF_GPUS_PER_NODES=4

EXE_NAME=test_mpi_kokkos.cuda


# Default behaviour : all mpi tasks are binded to the same GPU
mpirun --report-bindings ./$EXE_NAME

echo " "
echo "##############################################"
echo "##############################################"
echo "##############################################"
echo " "

# Nominal bahavior: each mpi task binded to a different GPU
mpirun --report-bindings ./$EXE_NAME --ndevices=$NUMBER_OF_GPUS_PER_NODES
