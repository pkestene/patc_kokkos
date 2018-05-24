# build OpenMP CPU
make KOKKOS_DEVICES=OpenMP

# build Cuda
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Kepler37
make -j 4 KOKKOS_DEVICES=Cuda KOKKOS_ARCH=Maxwell50
