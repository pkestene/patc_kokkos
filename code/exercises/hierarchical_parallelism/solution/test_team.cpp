#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<sys/time.h>
// Include Kokkos Headers
#include<Kokkos_Core.hpp>


#ifdef KOKKOS_ENABLE_CUDA
#include "CudaTimer.h"
#else
#include "OpenMPTimer.h"
#endif

using Device = Kokkos::DefaultExecutionSpace;
using Array1d = Kokkos::View<double*,Device>;

// ===============================================================
// ===============================================================
// ===============================================================
void test_dotprod(int length, int nrepeat, int nteams, int nthreads) {

  // Allocate Views
  Array1d B("B",length);
  Array1d C("C",length);
  
  // Initialize arrays
  Kokkos::parallel_for(length, KOKKOS_LAMBDA (const int& i) {
      //B(i) = cos(2*M_PI*i/100.);
      //C(i) = cos(2*M_PI*i/100.);
      B(i) = 1.0*i;
      C(i) = 1.0;
  });

  // Time computation
#ifdef KOKKOS_ENABLE_CUDA
  CudaTimer timer;
#else
  OpenMPTimer timer;
#endif

  using team_policy = Kokkos::TeamPolicy<>;
  using team_member = typename team_policy::member_type;

  const team_policy policy(nteams, nthreads);

  double sum = 0;
  
  timer.start();
  for(int k = 0; k < nrepeat; k++) {

    sum = 0.0;
    // Do dotprod
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
    Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA (const team_member& thread,
						   double& lsum) {
			      int i = thread.team_rank() + thread.league_rank() * thread.team_size();
			      // sweep array of size "length"
			      while (i<length) {
				lsum += B(i) * C(i);
				i += thread.league_size()*thread.team_size();
			      }
			      
			    },sum);
#endif
    
  }
  timer.stop();

  printf("Expected result = %f, actual result %f | error=%f\n",(length-1.0)*length/2,fabs(sum-(length-1.0)*length/2));
  
  // Print results
  double time_seconds = timer.elapsed();

  printf("#VectorLength  Time(s) TimePerIterations(s) size(MB) BW(GB/s)\n");
  printf("%13i %8lf %20.3e  %3.3f %3.3f\n",length,time_seconds,time_seconds/nrepeat,1.0e-6*length*3*8,1.0e-9*length*3*nrepeat*8/time_seconds);

} // test_dotprod

// ===============================================================
// ===============================================================
// ===============================================================
int main(int argc, char* argv[]) {

  // Parameters
  int length = 10000000; // length of vectors
  int nrepeat = 10;     // number of integration invocations

  int nteams = 4; // default number of (threads) teams
  int nthreads = 2; // default number of threads per team
  
  // Read command line arguments
  for(int i=0; i<argc; i++) {
    if( strcmp(argv[i], "-l") == 0) {
      length = atoi(argv[++i]);
    } else if( strcmp(argv[i], "-nrepeat") == 0) {
      nrepeat = atoi(argv[++i]);
    } else if( strcmp(argv[i], "-nteams") == 0) {
      nteams = atoi(argv[++i]);
    } else if( strcmp(argv[i], "-nthreads") == 0) {
      nthreads = atoi(argv[++i]);
    } else if( (strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
      printf("SAXPY Options:\n");
      printf("  -l <int>:         length of vectors (default: 10000000)\n");
      printf("  -nrepeat <int>:   number of integration invocations (default: 10)\n");
      printf("  -nteams <int>:    number of (thread) teams\n");
      printf("  -nthreads <int>:  number of threads per teams\n");
      printf("  -help (-h):       print this message\n");
    }
  }
  
  //Initialize Kokkos
  Kokkos::initialize(argc,argv);

  // run test
  test_dotprod(length, nrepeat, nteams, nthreads);

  // Shutdown Kokkos
  Kokkos::finalize();
}
