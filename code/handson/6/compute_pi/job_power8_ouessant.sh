#!/bin/bash
#BSUB -x
#BSUB -J test_compute_pi_openmp_cpu
#BSUB -n 1
#BSUB -o test_compute_pi_openmp_cpu.%J.out
#BSUB -q computet1
#BSUB -a "p8aff(160,8,8,pack)"
#BSUB -R "span[ptile=10]"
#BSUB -W 00:05

# here we ask
# 160 threads per (MPI) task, only 1 task
# SMT is 8
# all 8 logical core are used
# remember that cpus_per_core must be smaller than SMT

# submit job using bsub < job_power8_ouessant.sh

./compute_pi.openmp -niter 100000000
