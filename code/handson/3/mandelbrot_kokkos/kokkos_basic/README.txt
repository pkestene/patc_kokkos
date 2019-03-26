Just edit :
- main.cpp
- mandelbrot.h

fill the blanks marked with /* TODO */


# minimal cmake information
mkdir build_openmp; cd build_openmp
cmake -DKOKKOS_ENABLE_OPENMP=ON -DKOKKOS_ENABLE_HWLOC=ON -DKOKKOS_ARCH=Power8 ..
