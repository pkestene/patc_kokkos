#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "common.h"
#include "params.h"
#include "DataContext.h"

#include "OpenMPTimer.h"

#include "laplace2d_serial_kernel.h"

// ========================================================================
// ========================================================================
void test_laplace(int NX, int NY, int iter_max)
{

  real_t tol = 1e-5;
  Params params(NX, NY, iter_max, tol);

  // allocate data context
  DataContext context(params);

  real_t *rhs = context.rhs;
  
  // set rhs
  for (int iy = 1; iy < NY-1; iy++) {
    for( int ix = 1; ix < NX-1; ix++ ) {
      const real_t x = -1.0 + (2.0*ix/(NX-1));
      const real_t y = -1.0 + (2.0*iy/(NY-1));
      rhs[iy*NX+ix] = expr(-10.0*(x*x + y*y));
    }
  }

  int ix_start = 1;
  int ix_end   = (NX - 1);
  
  int iy_start = 1;
  int iy_end   = (NY - 1);

  // serial computation
  OpenMPTimer timer;
  timer.start();
  poisson2d_serial( context, params );
  timer.stop();
  real_t runtime_serial = timer.elapsed();

  printf("%dx%d: %8.4f secondes\n",NX,NY,runtime_serial);

} // test_laplace

// ========================================================================
// ========================================================================
int main(int argc, char* argv[])
{

  int NX = 512;
  int NY = 512;
  int iter_max = 1000;

  if (argc > 1) {
    NX = std::atoi(argv[1]);
    NY = std::atoi(argv[1]);
  }

  if (argc > 2)
    iter_max = std::atoi(argv[2]);
  
  test_laplace(NX,NY,iter_max);
  
  return 0;

}
