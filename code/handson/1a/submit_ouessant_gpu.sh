#!/bin/bash
#BSUB -x
#BSUB -gpu "num=1:mode=exclusive_process:mps=no:j_exclusive=yes"
#BSUB -J test_device_query
#BSUB -o test_device_query.%J.out
#BSUB -n 1
#BSUB -q computet1
#BSUB -W 00:05

module load at/11.0 cuda/9.2

echo "########################################################"
./deviceQuery
echo "########################################################"
env | grep CUDA

