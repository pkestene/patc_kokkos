#include <cmath>
#include <cstring>
#include <cstdio>
#include <cassert>

#include "common.h"
#include "params.h"
#include "DataContext.h"
#include "DataContextKokkos.h"

#include "OpenMPTimer.h"
#ifdef KOKKOS_ENABLE_CUDA
#include "CudaTimer.h"
#endif // KOKKOS_ENABLE_CUDA

#include "laplace2d_serial_kernel.h"
#include "laplace2d_kokkos_kernel.h"

#include "check_results.h"

#ifdef KOKKOS_ENABLE_CUDA
#include <cuda.h>
#endif // KOKKOS_ENABLE_CUDA

#include <mpi.h>

// ========================================================================
// ========================================================================
void test_laplace(int NX, int NY, int iter_max)
{

  int rank, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  /*
   * Inititalize some parameters.
   */

#ifdef USE_DOUBLE
  real_t tol = 1e-5;
#else
  real_t tol = 1e-4;
#endif
  
  real_t runtime_serial = 0;
  real_t runtime_kokkos = 0;

  /*
   * Perform serial reference computation.
   */
  if (rank == 0)
    printf("Jacobi relaxation Calculation: %d x %d mesh (serial on rank 0)\n", NY, NX);

  // make sure we have correct bounds (for MPI version only)
  // ix_start/ix_end/iy_start/iy_end not used for serial computations
  int ix_start = 1;
  int ix_end   = (NX - 1);
  
  // Ensure correctness if NY%nRanks != 0
  int chunk_size = ceil( (1.0*NY)/nRanks );
  
  int iy_start = rank * chunk_size;
  int iy_end   = iy_start + chunk_size;
  
  // Do not process boundaries
  iy_start = std::max( iy_start,      1 );
  iy_end   = std::min( iy_end  , NY - 1 );
  
  Params params(NX, NY, iter_max, tol, ix_start, ix_end, iy_start, iy_end);
  
  // allocate data context for serial computation
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
  OpenMPTimer timer;
  timer.start();
  poisson2d_serial( context, params );
  timer.stop();
  runtime_serial = timer.elapsed();
  

  /*
   * Now do parallel computation with 
   * - MPI for distributed paralellism 
   * - Kokkos for node-level parallelism
   */
  
  //Wait for all processes to ensure correct timing of the parallel version
  MPI_Barrier( MPI_COMM_WORLD );
  if ( rank == 0) printf("Jacobi relaxation Calculation: %d x %d mesh (MPI+Kokkos)\n", NY, NX);

  // kokkos computation context
  DataContextKokkos context_kokkos(params);

  // initialize context
  Kokkos::parallel_for( NX*NY, KOKKOS_LAMBDA(const int index) {      
      int i,j;
      index2coord(index,i,j,NX,NY);
      context_kokkos.A   (i,j) = 0;
      context_kokkos.Anew(i,j) = 0;
    });

  Kokkos::parallel_for( NX*NY, KOKKOS_LAMBDA(const int index) {    
      int i,j;
      index2coord(index,i,j,NX,NY);
      const real_t x = -1.0 + (2.0*i/(NX-1));
      const real_t y = -1.0 + (2.0*j/(NY-1));
      context_kokkos.rhs(i,j) = expr(-10.0*(x*x + y*y));
    });

#ifdef KOKKOS_ENABLE_CUDA
  CudaTimer timer2;
#else
  OpenMPTimer timer2;
#endif
  timer2.start();
  poisson2d_kokkos( context_kokkos, params );
  timer2.stop();
  runtime_kokkos = timer2.elapsed();

  check_results(context, context_kokkos, params);

  if (rank == 0) {
    printf("serial %dx%d: %8.4f secondes\n",NX,NY,runtime_serial);
    printf("kokkos %dx%d: %8.4f secondes\n",NX,NY,runtime_kokkos);
  }

} // test_laplace

// ========================================================================
// ========================================================================
int main(int argc, char* argv[])
{

  /* Init MPI */
  MPI_Init(&argc, &argv);

  /*
   * Initialize kokkos
   */
  // this regular initialization: let Kokkos use hwloc to decide
  // - how many OpenMP threads ? Default is to use all available core on node
  // - how GPU and MPI tasks are mapped
  Kokkos::initialize(argc, argv);
// #ifdef CUDA
//   // regular init
//   Kokkos::initialize(argc, argv);
// #else // OpenMP
//   // for OpenMP we might want to take control here
//   Kokkos::OpenMP::Initialize
// #endif
    
  int rank, nRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
  
  if (rank==0) {
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
  }

#ifdef KOKKOS_ENABLE_CUDA
  /*
   * Just checking, which MPI rank mapped with which GPU ?
   */
  int cudaDeviceId;
  cudaGetDevice(&cudaDeviceId);
  std::cout << "I'm MPI task #" << rank << " pinned to GPU #" << cudaDeviceId << "\n";
#endif // KOKKOS_ENABLE_CUDA

  int NX = 512;
  int NY = 512;
  int iter_max = 1000;

  if (argc > 1) {
    NX = atoi(argv[1]);
    NY = atoi(argv[1]);
  }

  if (argc > 2)
    iter_max = atoi(argv[2]);
  
  test_laplace(NX,NY,iter_max);
  
  Kokkos::finalize();

  MPI_Finalize();
  
  return 0;

}
