#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<sys/time.h>
#include <vector>
#include <fstream>      // std::ofstream

// Include Kokkos Headers
#include<Kokkos_Core.hpp>

// make the compiler ignore an unused variable
#ifndef UNUSED
#define UNUSED(x) ((void)(x))
#endif

#ifdef USE_DOUBLE
using real_t = double;
#else
using real_t = float;
#endif // USE_DOUBLE

#ifdef KOKKOS_ENABLE_CUDA
#include "CudaTimer.h"
using Timer = CudaTimer;
#elif defined(KOKKOS_ENABLE_OPENMP)
#include "OpenMPTimer.h"
using Timer = OpenMPTimer;
#else
#include "SimpleTimer.h"
using Timer = SimpleTimer;
#endif

using Device = Kokkos::DefaultExecutionSpace;
using DataArray = Kokkos::View<real_t***, Device>;

// ===============================================================
// ===============================================================
KOKKOS_INLINE_FUNCTION
void index2coord(int index,
                 int &i, int &j,
                 int Nx, int Ny)
{
  UNUSED(Nx);
#ifdef KOKKOS_ENABLE_CUDA
  j = index / Nx;
  i = index - j*Nx;
#else
  i = index / Ny;
  j = index - i*Ny;
#endif
} // index2coord - 2d

// ===============================================================
// ===============================================================
KOKKOS_INLINE_FUNCTION
void index2coord(int index,
                 int &i, int &j, int &k,
                 int Nx, int Ny, int Nz)
{
  UNUSED(Nx);
  UNUSED(Nz);
#ifdef KOKKOS_ENABLE_CUDA
  int NxNy = Nx*Ny;
  k = index / NxNy;
  j = (index - k*NxNy) / Nx;
  i = index - j*Nx - k*NxNy;
#else
  int NyNz = Ny*Nz;
  i = index / NyNz;
  j = (index - i*NyNz) / Nz;
  k = index - j*Nz - i*NyNz;
#endif
} // index2coord - 3d

// ===============================================================
// ===============================================================
KOKKOS_INLINE_FUNCTION
int coord2index(int i,  int j,  int k,
                int Nx, int Ny, int Nz)
{
  UNUSED(Nx);
  UNUSED(Nz);
#ifdef KOKKOS_ENABLE_CUDA
  return i + Nx*j + Nx*Ny*k; // left layout
#else
  return k + Nz*j + Nz*Ny*i; // right layout
#endif
}

// version 1
#include "test_stencil_3d_flat.cpp"

// version 2
//#include "test_stencil_3d_flat_vector.cpp"

// version 3
//#include "test_stencil_3d_range.cpp"

// version 4
//#inclue "test_stencil_3d_range_vector.cpp"

// version 5
//#include "test_stencil_3d_range_vector2.cpp"

// ===============================================================
// ===============================================================
// ===============================================================
int main(int argc, char* argv[]) {

  // Parameters
  int n       = 256;  // 3d array linear size 
  int nrepeat =  10;  // number of kernel invocations
  int nteams  =  64;  // default number of teams (for TeamPolicy)
  
  // Read command line arguments
  for(int i=0; i<argc; i++) {
    if( strcmp(argv[i], "-n") == 0) {
      n = atoi(argv[++i]);
    } else if( strcmp(argv[i], "-nrepeat") == 0) {
      nrepeat = atoi(argv[++i]);
    } else if( strcmp(argv[i], "-nteams") == 0) {
      nteams = atoi(argv[++i]);
    } else if( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
      printf("STENCIL 3D Options:\n");
      printf("  -n <int>:         3d linear size (default: 256)\n");
      printf("  -nrepeat <int>:   number of integration invocations (default: 10)\n");
      printf("  -help (-h):       print this message\n");
      return EXIT_SUCCESS;
    }
  }
  
  //Initialize Kokkos
  Kokkos::initialize(argc,argv);

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
  std::cout << msg.str();
  std::cout << "##########################\n";
  
  // run test
  std::cout << "========================================\n";
  std::cout << "reference naive test using 1d flat range\n";
  test_stencil_3d_flat(n, nrepeat);
  
  // std::cout << "========================================\n";
  // std::cout << "reference naive test using 2d flat range and vectorization (no views)\n";
  // test_stencil_3d_flat_vector(n, nrepeat,false);
  
  // std::cout << "========================================\n";
  // std::cout << "reference naive test using 2d flat range and vectorization (with views)\n";
  // test_stencil_3d_flat_vector(n, nrepeat,true);
  
  // std::cout << "========================================\n";
  // std::cout << "reference naive test using 3d range\n";
  // test_stencil_3d_range(n, nrepeat);
  
  // std::cout << "========================================\n";
  // std::cout << "reference naive test using 3d range and vectorization\n";
  // test_stencil_3d_range_vector(n, nrepeat);
  
  // std::cout << "========================================\n";
  // std::cout << "reference naive test using 3d range and vectorization with team policy\n";
  // test_stencil_3d_range_vector2(n, nrepeat,nteams);
  
  // Shutdown Kokkos
  Kokkos::finalize();
  
  return EXIT_SUCCESS;
}
