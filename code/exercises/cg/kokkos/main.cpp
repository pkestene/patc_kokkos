#include <cstdlib>
#include <cstdio>

#include "kokkos_shared.h"

#include "vector.h"
#include "vector_functions.h"
#include "matrix.h"
#include "matrix_functions.h"

#define N 100
#define MAX_ITERS 100
#define TOL 1e-12

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);
  
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
  
  matrix A(N);
  vector x(A.num_rows),b(A.num_rows);
  vector r(A.num_rows),p(A.num_rows),Ap(A.num_rows);
  
  double one=1.0, zero=0.0;
  double normr, rtrans, oldtrans, p_ap_dot , alpha, beta;
  int iter=0;
  
  // init matrix
  A.init();
  
  //printf("Rows: %d, nnz: %d\n", A.num_rows, A.row_offsets[A.num_rows]);
  printf("Rows: %d, nnz: %d\n", A.num_rows, A.nnz);
  
  x.init(100000.0);
  b.init(1.0);
  
  waxpby(one, x, zero, x, p);
  matvec(A,p,Ap);
  waxpby(one, b, -one, Ap, r);
  
  rtrans = dot(r,r);
  normr = sqrt(rtrans);
  
  double st = omp_get_wtime();
  do {
    if(iter==0) {
      waxpby(one,r,zero,r,p);
    } else {
      oldtrans=rtrans;
      rtrans = dot(r,r);
      beta = rtrans/oldtrans;
      waxpby(one,r,beta,p,p);
    }
    
    normr=sqrt(rtrans);
  
    matvec(A,p,Ap);
    p_ap_dot = dot(Ap,p);

    alpha = rtrans/p_ap_dot;

    waxpby(one,x,alpha,p,x);
    waxpby(one,r,-alpha,Ap,r);

    if(iter%10==0)
      printf("Iteration: %d, Tolerance: %.4e\n", iter, normr);
    iter++;
  } while(iter<MAX_ITERS && normr>TOL);
  double et = omp_get_wtime();
  
  printf("Total Iterations: %d\n", iter);
  printf("Total Time: %lf s\n", (et-st));
  
  Kokkos::finalize();
  
  return 0;
}
