#!/bin/bash
#BSUB -x
#BSUB -J test_saxpy_omp
#BSUB -n 1
#BSUB -o test_saxpy_omp.%J.out
#BSUB -q computet1
#BSUB -a "p8aff(20,8,8,pack)"
#BSUB -R "span[ptile=1]"
#BSUB -W 00:05

# here we ask
# 20 threads per (MPI) task, only 1 task
# SMT is 8
# all 8 logical core can be used
# remember that cpus_per_core must be smaller than SMT

./saxpy.host -l  10000000 -nrepeat 100
./saxpy.host -l  20000000 -nrepeat 100
./saxpy.host -l  50000000 -nrepeat 100
./saxpy.host -l 100000000 -nrepeat 100
./saxpy.host -l 200000000 -nrepeat 100

