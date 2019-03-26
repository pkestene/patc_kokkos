Just edit :
- main.cpp
- mandelbrot.h

fill the blanks marked with /* TODO */


# minimal cmake information for Power8
mkdir build_openmp; cd build_openmp
cmake -DKOKKOS_ENABLE_OPENMP=ON -DKOKKOS_ENABLE_HWLOC=ON -DKOKKOS_ARCH=Power8 ..

# minimal cmake info for Kepler K80
mkdir build_k80; cd build_k80
cmake -DKOKKOS_ENABLE_CUDA=ON -DKOKKOS_ENABLE_OPENMP=ON -DKOKKOS_ENABLE_HWLOC=ON -DKOKKOS_ARCH=Kepler37 ..

# you may also turn CMAKE_BUILD_TYPE to ON (set optimization flags)

