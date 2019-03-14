#!/bin/bash
#BSUB -x
#BSUB -n 1
# reserve 1 GPU per node
#BSUB -gpu "num=1:mode=exclusive_process:mps=no:j_exclusive=yes"
#BSUB -J stencil_cuda
#BSUB -o stencil_cuda.%J.log
#BSUB -q computet1
#BSUB -a p8aff(1,8,1,balance)
# 1 tasks per nodes 
#BSUB -R "span[ptile=1]"
#BSUB -W 01:00

# launch this script with bsub < job_cuda_ouessant.sh

./stencil.cuda -b
