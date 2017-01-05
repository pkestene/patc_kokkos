#!/bin/bash
#BSUB -x
#BSUB -J test_kokkos
#BSUB -n 1
#BSUB -o test_kokkos.%J.out
#BSUB -q compute
#BSUB -R "affinity[core(1):cpubind=core:distribute=pack]"
#BSUB -R "span[ptile=10]"
#BSUB -W 00:05

module load gcc/4.8/ompi/1.10
module load cuda/8.0

# here setup env for kokkos

MY_HOME=/pwrhome/rech/mhb/rmhb001/

EXE_NAME=.....

mpirun -rf ${LSB_RANK_HOSTFILE} $EXE_NAME

