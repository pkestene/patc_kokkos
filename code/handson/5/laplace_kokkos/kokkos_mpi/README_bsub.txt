# to submit a job using only Power8 (OpenMP)
# default is 10 MPI tasks, perform bench to increase
bsub < submit_ouessant_cpu.sh

# to submit a job using GPU's (1 MPI task per GPU)
bsub < submit_oeussant_gpu.sh
