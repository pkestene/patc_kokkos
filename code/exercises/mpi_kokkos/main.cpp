#include <cstdio>
#include <iostream>
#include <fstream>

#include <omp.h>
#include <mpi.h>
#include "kokkos_shared.h"

#ifdef CUDA
#include <cuda.h>
#endif

#include <unistd.h>

using namespace std;

int main(int argc, char* argv[]) {

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
    
  /*
   * Initialize kokkos (host + device)
   */
  Kokkos::initialize(argc, argv);

  int rank, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
  
  {
    std::cout << "##########################\n";
    std::cout << "KOKKOS CONFIG             \n";
    std::cout << "##########################\n";
    
    std::ostringstream msg;
    std::cout << "Kokkos configuration" << std::endl;
    if ( Kokkos::hwloc::available() ) {
      msg << "hwloc( NUMA[" << Kokkos::hwloc::get_available_numa_count()
          << "] x CORE["    << Kokkos::hwloc::get_available_cores_per_numa()
          << "] x HT["      << Kokkos::hwloc::get_available_threads_per_core()
          << "] )"
          << std::endl ;
    }
    Kokkos::print_configuration( msg );
#if defined( CUDA )
    Kokkos::Cuda::print_configuration( msg );
#else
    Kokkos::OpenMP::print_configuration( msg );
#endif
    std::cout << msg.str();
    std::cout << "##########################\n";
  }

#ifdef CUDA
  int cudaDeviceId;
  cudaGetDevice(&cudaDeviceId);
  std::cout << "I'm MPI task #" << rank << " pinned to GPU #" << cudaDeviceId << "\n";
#endif  

  Kokkos::finalize();

  MPI_Finalize();
  
  return 0;
}
