# Don't forget to export KOKKOS_PATH, e.g.
export KOKKOS_PATH=~/kokkos-tutorial/kokkos

# build OpenMP CPU
make KOKKOS_DEVICES=OpenMP

# build Cuda
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Kepler37
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Maxwell50
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Pascal60

# example run on ouessant (login nodes)
mpirun -np 4 ./laplace2d_kokkos_mpi.cuda 2048 100 --ndevices=4
