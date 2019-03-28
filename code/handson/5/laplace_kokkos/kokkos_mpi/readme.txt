# modules
module load at/10.0 cuda/9.2 gcc/4.8/ompi/2.1

# Don't forget to export KOKKOS_PATH, e.g.
export KOKKOS_PATH=~/kokkos-tutorial/kokkos

# build OpenMP CPU
make KOKKOS_DEVICES=OpenMP

# example run on ouessant (login nodes - CPU - Power8)
export OMP_NUM_THREADS=10
mpirun -np 4 ./laplace2d_kokkos_mpi.omp 2048 100

# build Cuda
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Kepler37
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Maxwell50
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Pascal60

# example run on ouessant (login nodes - GPU - K80)
mpirun -np 4 ./laplace2d_kokkos_mpi.cuda 2048 100 --ndevices=4

# To run on compute node, use submit_ouessant_gpu.sh
