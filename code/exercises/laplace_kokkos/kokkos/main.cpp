#include <cmath>
#include <cstring>
#include <cstdio>
#include <cassert>

#include "common.h"
#include "params.h"
#include "DataContext.h"
#include "DataContextKokkos.h"

#include "OpenMPTimer.h"

#include "laplace2d_serial_kernel.h"
#include "laplace2d_kokkos_kernel.h"

#include "check_results.h"

int main(int argc, char* argv[])
{

  /*
   * Initialize kokkos (host + device)
   */
#ifdef CUDA
  // Initialize Host mirror device
  Kokkos::HostSpace::execution_space::initialize();
  //const unsigned device_count = Kokkos::Cuda::detect_device_count();

  // Use the first device:
  Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
#else // OpenMP CPU
  Kokkos::initialize(argc, argv);
#endif

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
#if defined( CUDA )
    Kokkos::Cuda::print_configuration( msg );
#else
    Kokkos::OpenMP::print_configuration( msg );
#endif
    std::cout << msg.str();
    std::cout << "##########################\n";
  }


  int NX = 1024;
  int NY = 1024;
  int iter_max = 1000;

#ifdef USE_DOUBLE
  real_t tol = 1e-5;
#else
  real_t tol = 1e-4;
#endif
  
  Params params(NX, NY, iter_max, tol);

#ifdef USE_DOUBLE
  printf("DOUBLE PRECISION computations\n");
#else
  printf("SINGLE PRECISION computations\n");
#endif
  
  printf("Jacobi relaxation Calculation: %d x %d mesh\n", NY, NX);
  
  // allocate data context
  DataContext context(params);

  memset(context.A,    0, NY * NX * sizeof(real_t));
  memset(context.Aref, 0, NY * NX * sizeof(real_t));

  real_t *rhs = context.rhs;
  
  // set rhs
  for (int iy = 1; iy < NY-1; iy++) {
    for( int ix = 1; ix < NX-1; ix++ ) {
      const real_t x = -1.0 + (2.0*ix/(NX-1));
      const real_t y = -1.0 + (2.0*iy/(NY-1));
      rhs[iy*NX+ix] = expr(-10.0*(x*x + y*y));
    }
  }

  // serial computation
  printf("Calculate reference solution and time serial execution.\n");
  OpenMPTimer timer;
  timer.start();
  poisson2d_serial( context, params );
  timer.stop();
  real_t runtime_serial = timer.elapsed();

  // kokkos computation
  printf("Parallel execution with kokkos.\n");
  DataContextKokkos context_kokkos(params);

  // initialize context
  Kokkos::parallel_for( NX*NY, KOKKOS_LAMBDA(const int index) {      
      context_kokkos.A   (index) = 0;
      context_kokkos.Anew(index) = 0;
    });

  Kokkos::parallel_for( NX*NY, KOKKOS_LAMBDA(const int index) {    
      int ix,iy;
      index2coord(index,ix,iy,NX,NY);
      const real_t x = -1.0 + (2.0*ix/(NX-1));
      const real_t y = -1.0 + (2.0*iy/(NY-1));
      context_kokkos.rhs(index) = expr(-10.0*(x*x + y*y));
    });
  
  timer.reset();
  timer.start();
  poisson2d_kokkos( context_kokkos, params );
  timer.stop();
  real_t runtime_kokkos = timer.elapsed();

  if ( check_results(context, context_kokkos, params) )
    printf("Serial and Kokkos results match !\n");
  else
    printf("Serial and Kokkos results don't match !\n");    
  
  printf("serial %dx%d: %8.4f secondes\n",NX,NY,runtime_serial);
  printf("kokkos %dx%d: %8.4f secondes\n",NX,NY,runtime_kokkos);

#ifdef CUDA
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
#else
  Kokkos::finalize();
#endif
  
  return 0;

}
