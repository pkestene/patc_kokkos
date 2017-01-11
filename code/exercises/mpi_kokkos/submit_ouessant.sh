#!/bin/bash
#BSUB -x
#BSUB -J test_mpi_kokkos
#BSUB -n 4
#BSUB -o test_mpi_kokkos.%J.out
#BSUB -q compute
#BSUB -R "affinity[core(5):cpubind=core]"
#BSUB -R "select[ngpus>0] rusage [ngpus_shared=1]"
#BSUB -W 00:05


module load gcc/4.8/ompi/1.10
module load cuda/8.0

# This variable is normally set by the job scheduler
# As of January, 9th 2017, we enforce its value here
# to make sure all GPU devices can be used by our job.
CUDA_VISIBLE_DEVICES=0,1,2,3

NUMBER_OF_GPUS_PER_NODES=4

EXE_NAME=test_mpi_kokkos.cuda


# all mpi tasks are binded to the same GPU
mpirun -n ${LSB_DJOB_NUMPROC} ./$EXE_NAME

echo " "
echo "##############################################"
echo "##############################################"
echo "##############################################"
echo " "

# each mpi task are binded to a different GPU
mpirun -n ${LSB_DJOB_NUMPROC} ./$EXE_NAME --ndevices=$NUMBER_OF_GPUS_PER_NODES
