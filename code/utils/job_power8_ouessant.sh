#!/bin/bash
#BSUB -x
#BSUB -J test_kokkos_cpu
#BSUB -n 1
#BSUB -o test_kokkos_cpu.%J.out
#BSUB -q computet1
#BSUB -a "p8aff(160,8,8,pack)"
#BSUB -R "span[ptile=1]"
#BSUB -W 00:05

# here we ask
# 160 threads per (MPI) task, only 1 task
# SMT is 8
# all 8 logical core are used
# remember that cpus_per_core must be smaller than SMT
# span[ptile=1] means our job is requesting 1 task per node

# sometimes, it may be needed to load modules for runtime path
#module load at/10
#module load cuda/9.2
#module load gcc/4.8/ompi/2.1


# how to submit this job ? 
# bsub < job_power8_ouessant.sh

./NAME_OF_EXECUTABLE
