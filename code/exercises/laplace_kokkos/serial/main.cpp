#include <cmath>
#include <cstring>
#include <cstdio>
#include <cassert>

#include "common.h"
#include "params.h"
#include "DataContext.h"

#include "laplace2d_kernel.h"

int main(int argc, char* argv[])
{

  int NX = 1024;
  int NY = 1024;
  int iter_max = 1000;
  real_t tol = 1e-5;
  Params params(NX, NY, iter_max, tol);

  // allocate data context
  DataContext context(params);

  
  
  return 0;

}
