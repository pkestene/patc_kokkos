This example is loosely adapted from Nvidia's OpenAcc tutorial:
https://github.com/NVIDIA-OpenACC-Course/nvidia-advanced-openacc-course-sources

- task1: parallelize computation, but not memory (i.e. every MPI task allocates the entier domain). Parallelization strategy is divide by chunk along y-axis.
Each MPI task only needs to modify loop indexes iy_start, iy_end to compute its assigned chunk.
