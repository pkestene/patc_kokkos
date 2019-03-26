#!/bin/bash
#BSUB -x
#BSUB -n 1
#BSUB -o test_kokkos_gpu.%J.log
#BSUB -q computet1
#BSUB -a p8aff(1,1,1,balance)
# reserve 1 GPU per node
#BSUB -gpu "num=1:mode=exclusive_process:mps=no:j_exclusive=yes"
#BSUB -J test_kokkos_gpu
# 1 tasks per node
#BSUB -R "span[ptile=1]"
#BSUB -W 01:00

# here we ask
# 1 thread per (MPI) task, only 1 task
# SMT is 1
# all 1 logical core used
# remember that cpus_per_core must be smaller than SMT
# span[ptile=1] means our job is requesting 1 task per node

# how to submit this job ? 
# bsub < job_cuda_ouessant.sh

./NAME_OF_EXECUTABLE
