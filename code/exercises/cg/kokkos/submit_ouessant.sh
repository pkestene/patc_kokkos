#!/bin/bash
#BSUB -x
#BSUB -J test_kokkos
#BSUB -n 1
#BSUB -o test_kokkos.%J.out
#BSUB -q compute
#BSUB -R "affinity[core(1):cpubind=core:distribute=pack]"
#BSUB -R "span[ptile=10]"
#BSUB -R "select[ngpus>0] rusage [ngpus_shared=1]"
#BSUB -W 00:05

module load at/10.0 ompi/2.1
module load cuda/9.0

# This variable is normally set by the job scheduler
# As of January, 9th 2017, we enforce its value here
# to make sure all GPU devices can be used by our job.
CUDA_VISIBLE_DEVICES=0,1,2,3

# init the name of your executable
./cg.cuda



