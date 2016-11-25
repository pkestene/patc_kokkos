A minimalistic kokkos example to compute Mandelbrot set and illustrate asynchronous memory copy (i.e. overlap between a computationnal functor and a deep copy operation).

Three versions provided:

* basic : Mandelbrot set is computed in a single Kokkos functor. Can be used
  with either OpenMP or Cuda backend.

* pipeline0 performs computations piece by piece

* pipeline1 (only meaningfull when used with CUDA+OpenMP) performs computations piece by piece, but the loop over the pieces
  is parallelized using an OpenMP Kokkos functor, so that the different pieces
  can be computed in different Cuda streams.


# What is kokkos ?

A modern C++ based programming model for HPC applications designed for portability across multiple hardware architectures (multicore, GPU, KNL, Power8, ...) and also providing as efficient as possible performance.

# Build the basic version

0. Need to have installed [kokkos](https://github.com/kokkos/kokkos)

   * Kokkos backend can be CUDA, OpenMP, ...
   * Compiler can be nvcc_wrapper, g++, xlc++, ...

1. Set env variable KOKKOS_PATH to the root directory where Kokkos is installed

2. cd basic; make

3. run

   ./mandelbrot.omp (or ./mandelbrot.cuda)


With default parameters (image of size 8192x8192), some performance of the basic version:
 * Nvidia K80 : 1.5 seconds 
 * Power8 (g++ 4.8.5)
   * 20  threads : 50.1 seconds
   * 40  threads : 27.7 seconds
   * 60  threads : 16.9 seconds
   * 160 threads : 10.5 seconds
   

NB: version pipeline1 require kokkos to be configured with both CUDA and OPENMP
backends, and lambda function enabled. Example command line configuration:

generate_makefile.bash --with-cuda --arch=${YOUR_CUDA_ARCH} --prefix=${SOMEWHERE} --with-cuda-options=enable_lambda --with-openmp
