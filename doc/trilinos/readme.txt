#################################################
1. configure_tpetra_kokkos_cuda_nvcc_wrapper.sh
#################################################

This script allows you to build current trilinos (git branch master / develop)
https://github.com/trilinos/Trilinos

How ?
In trilinos sources : mkdir build; cd build
./configure_tpetra_kokkos_cuda_nvcc_wrapper.sh

Then: make; make install

#################################################
2. module environment
#################################################

To simplify the use of trilinos, you may use a "module"

Subdirectory module contains an example.
Just modify file modulefiles/trilinos/kokkos_dev by changing "topdir"
which should point to your installed trilinos.

Then you can do: "module use /home/.../modulefiles"
module load trilinos/kokkos_dev

Your environment is then setup for using trilinos.

You may still need to set OMPI_CXX to nvcc_wrapper.

#################################################
3. tpetra_example 
#################################################

This directory contains a minimal tpetra example

#################################################
4. trilinos_kokkos_examples
#################################################

This directory contains examples found on the official Trilinos tutorial
https://github.com/trilinos/Trilinos_tutorial/wiki/TrilinosHandsOnTutorial#examples-illustrating-kokkos-thread-scalable-expressions-to-generate-tpetra-objects

They have just been slightly modify for lambda --> [=] replaced by KOKKOS_LAMBDA


