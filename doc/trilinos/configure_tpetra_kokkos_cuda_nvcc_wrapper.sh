#!/bin/bash
#
# This script uses CUDA, OpenMP, and MPI.
#
# Before invoking this script, set the OMPI_CXX environment variable
# to point to nvcc_wrapper, wherever it happens to live.  (If you use
# an MPI implementation other than OpenMPI, set the corresponding
# environment variable instead.)
#
# You will need to export the following variable pointing to nvcc_wrapper found
# in Trilinos sources.
# export OMPI_CXX=/home/pkestene/install/trilinos/github/trilinos/packages/kokkos/config/nvcc_wrapper
#
#

rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*
EXTRA_ARGS=$@
MPI_PATH="/usr"
CUDA_PATH="/usr/local/cuda-8.0/"

#
# As long as there are any .cu files in Trilinos, we'll need to set
# CUDA_NVCC_FLAGS.  If Trilinos gets rid of all of its .cu files and
# lets nvcc_wrapper handle them as .cpp files, then we won't need to
# set CUDA_NVCC_FLAGS.  As it is, given that we need to set
# CUDA_NVCC_FLAGS, we must make sure that they are the same flags as
# nvcc_wrapper passes to nvcc.
#
#CUDA_NVCC_FLAGS="-gencode;arch=compute_50,code=sm_50;--verbose;-I${MPI_PATH}/include"
#CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-Xcompiler;-Wall,-ansi,-fopenmp"
#CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-O3;-DKOKKOS_USE_CUDA_UVM"
#CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-O3"

#CUDA_NVCC_FLAGS="-arch=sm_50 -Xcompiler -Wall,-ansi,-fopenmp -DKOKKOS_USE_CUDA_UVM -DKOKKOS_CUDA_USE_LAMBDA -expt-extended-lambda"

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH="/home/pkestene/local/trilinos_kokkos_dev" \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D CMAKE_CXX_FLAGS:STRING="-g -Wall -std=c++11 -arch=sm_50 -Xcompiler -Wall,-ansi,-fopenmp,-std=c++11 -expt-extended-lambda -DKOKKOS_USE_CUDA_UVM -DKOKKOS_CUDA_USE_LAMBDA" \
  -D CMAKE_C_FLAGS:STRING="-g -Wall" \
  -D CMAKE_FORTRAN_FLAGS:STRING="" \
  -D CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS="" \
  -D Trilinos_ENABLE_Triutils=OFF \
  -D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=ON \
  -D Trilinos_ENABLE_DEBUG:BOOL=OFF \
  -D Trilinos_ENABLE_CHECKED_STL:BOOL=OFF \
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
  -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
  -D BUILD_SHARED_LIBS:BOOL=OFF \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  \
  \
  -D CMAKE_CXX_COMPILER:FILEPATH="${MPI_PATH}/bin/mpicxx" \
  -D CMAKE_C_COMPILER:FILEPATH="${MPI_PATH}/bin/mpicc" \
  -D MPI_CXX_COMPILER:FILEPATH="${MPI_PATH}/bin/mpicxx" \
  -D MPI_C_COMPILER:FILEPATH="${MPI_PATH}/bin/mpicc" \
  -D CMAKE_Fortran_COMPILER:FILEPATH="${MPI_PATH}/bin/mpif77" \
  -D MPI_EXEC:FILEPATH="${MPI_PATH}/bin/mpirun" \
  -D MPI_EXEC_POST_NUMPROCS_FLAGS:STRING="-bind-to;socket;--map-by;socket;env;CUDA_MANAGED_FORCE_DEVICE_ALLOC=1;CUDA_LAUNCH_BLOCKING=1;OMP_NUM_THREADS=2" \
  \
  \
  -D Trilinos_ENABLE_CXX11:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D Trilinos_ENABLE_OpenMP:BOOL=ON \
  -D Trilinos_ENABLE_ThreadPool:BOOL=ON \
  \
  \
  -D TPL_ENABLE_CUDA:BOOL=ON \
  -D CUDA_TOOLKIT_ROOT_DIR:FILEPATH="${CUDA_PATH}" \
  -D CUDA_PROPAGATE_HOST_FLAGS:BOOL=OFF \
  -D TPL_ENABLE_Thrust:BOOL=OFF \
  -D Thrust_INCLUDE_DIRS:FILEPATH="${CUDA_PATH}/include" \
  -D TPL_ENABLE_CUSPARSE:BOOL=OFF \
  -D TPL_ENABLE_Cusp:BOOL=OFF \
  -D Cusp_INCLUDE_DIRS="/home/crtrott/Software/cusp" \
  -D CUDA_VERBOSE_BUILD:BOOL=OFF \
  -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS} \
  \
  \
  -D TPL_ENABLE_HWLOC=ON \
  -D HWLOC_INCLUDE_DIRS="/usr/include" \
  -D HWLOC_LIBRARY_DIRS="/usr/lib/x86_64-linux-gnu/" \
  -D TPL_ENABLE_BinUtils=OFF \
  -D TPL_ENABLE_BLAS:STRING=ON \
  -D TPL_ENABLE_LAPACK:STRING=ON \
  -D TPL_ENABLE_MKL:STRING=OFF \
  -D TPL_ENABLE_HWLOC:STRING=ON \
  -D TPL_ENABLE_GTEST:STRING=ON \
  -D TPL_ENABLE_SuperLU=ON \
  -D TPL_ENABLE_BLAS=ON \
  -D TPL_ENABLE_LAPACK=ON \
  -D TPL_SuperLU_LIBRARIES="/usr/lib/x86_64-linux-gnu/libsuperlu.a" \
  -D TPL_SuperLU_INCLUDE_DIRS="/usr/include/superlu/" \
  \
  \
  -D Trilinos_Enable_Kokkos:BOOL=ON \
  -D Trilinos_ENABLE_KokkosCore:BOOL=ON \
  -D Trilinos_ENABLE_TeuchosKokkosCompat:BOOL=ON \
  -D Trilinos_ENABLE_KokkosContainers:BOOL=ON \
  -D Trilinos_ENABLE_TpetraKernels:BOOL=ON \
  -D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON \
  -D Trilinos_ENABLE_TeuchosKokkosComm:BOOL=ON \
  -D Trilinos_ENABLE_KokkosExample:BOOL=ON \
  -D Kokkos_ENABLE_EXAMPLES:BOOL=ON \
  -D Kokkos_ENABLE_TESTS:BOOL=OFF \
  -D KokkosClassic_DefaultNode:STRING="Kokkos::Compat::KokkosCudaWrapperNode" \
  -D TpetraClassic_ENABLE_OpenMPNode=OFF \
  -D TpetraClassic_ENABLE_TPINode=OFF \
  -D TpetraClassic_ENABLE_MKL=OFF \
  -D Kokkos_ENABLE_Cuda_UVM=ON \
  -D Kokkos_ENABLE_HWLOC= ON \
  \
  \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Teuchos_ENABLE_COMPLEX:BOOL=OFF \
  \
  \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
  -D Tpetra_ENABLE_KokkosCore=ON \
  -D Tpetra_ENABLE_Kokkos_DistObject=OFF \
  -D Tpetra_ENABLE_Kokkos_Refactor=ON \
  -D Tpetra_ENABLE_TESTS=OFF \
  -D Tpetra_ENABLE_EXAMPLES=ON \
  -D Tpetra_ENABLE_MPI_CUDA_RDMA:BOOL=ON \
  \
  \
  -D Trilinos_ENABLE_Belos=OFF \
  -D Trilinos_ENABLE_Amesos=OFF \
  -D Trilinos_ENABLE_Amesos2=OFF \
  -D Trilinos_ENABLE_Ifpack=OFF \
  -D Trilinos_ENABLE_Ifpack2=OFF \
  -D Trilinos_ENABLE_Epetra=OFF \
  -D Trilinos_ENABLE_EpetraExt=OFF \
  -D Trilinos_ENABLE_Zoltan=OFF \
  -D Trilinos_ENABLE_Zoltan2=OFF \
  -D Trilinos_ENABLE_MueLu=OFF \
  -D Belos_ENABLE_TESTS=ON \
  -D Belos_ENABLE_EXAMPLES=ON \
  -D MueLu_ENABLE_TESTS=OFF \
  -D MueLu_ENABLE_EXAMPLES=OFF \
  -D Ifpack2_ENABLE_TESTS=ON \
  -D Ifpack2_ENABLE_EXAMPLES=ON \
  $EXTRA_ARGS \
..

