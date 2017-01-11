#!/bin/bash
#BSUB -x
#BSUB -J poisson_kokkos                            # Job name
#BSUB -n 4                                         # total number of MPI task
#BSUB -o poisson_kokkos.%J.out                     # stdout filename
#BSUB -q compute                                   # queue name
#BSUB -R "affinity[core(5):cpubind=core]"          # affinity
#BSUB -R "select[ngpus>0] rusage [ngpus_shared=1]" # select GPU compute mode
#BSUB -W 00:05                                     # maximum runtime

module load gcc/4.8/ompi/1.10 
moudle load cuda/8.0
module load kokkos/cuda80_gnu485_dev_p100

OMP_DISPLAY_ENV=true

# for MPI + CUDA on Pascal P100
# option --ndevices tells Kokkos to activate all GPU on node seen by hwloc
EXE_NAME=laplace2d_kokkos.cuda
mpirun -n ${LSB_DJOB_NUMPROC} ./$EXE_NAME --ndevices=$NUMBER_OF_GPUS_PER_NODES
