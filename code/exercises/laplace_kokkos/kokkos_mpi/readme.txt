# build OpenMP CPU

# build Cuda
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Kepler37
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Maxwell50
