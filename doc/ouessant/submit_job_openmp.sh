#!/bin/bash
#BSUB -x
#BSUB -J test_saxpy_omp
#BSUB -n 1
#BSUB -o test_saxpy_omp.%J.out
#BSUB -q computet1
#BSUB -R "affinity[core(1):cpubind=core:distribute=pack]"
#BSUB -a "p8aff(20,8,8,pack)"
#BSUB -R "span[ptile=10]"
#BSUB -W 00:05

# module load gcc/4.8/ompi/2.1
# module load cuda/9.0

# here we ask
# 20 threads per (MPI) task, only 1 task
# SMT is 8
# all 8 logical core can be used
# remember that cpus_per_core must be smaller than SMT

./saxpy.host

