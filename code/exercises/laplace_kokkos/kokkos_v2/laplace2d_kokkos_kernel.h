#ifndef LAPLACE2D_KOKKOS_KERNEL_H_
#define LAPLACE2D_KOKKOS_KERNEL_H_

#include "params.h"
#include "DataContextKokkos.h"

#include <limits> // for std::numeric_limits
#ifdef __CUDA_ARCH__
#include <math_constants.h>
#endif // __CUDA_ARCH__

/**
 * Kokkos kernel (reduction)
 */
class Poisson2dKokkosKernel {

public:
  Poisson2dKokkosKernel( DataContextKokkos& context_, Params& params_):
    A(context_.A), Anew(context_.Anew), rhs(context_.rhs), params(params_)
  {}

  // this is a reduce (max) functor
  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index, real_t &error) const
  //void operator()(const int& index) const
  {
    int i,j;
    index2coord(index,i,j,params.NX,params.NY);

    const int NX = params.NX;
    const int NY = params.NY;
    
    if ( j >= 1 and j < NY-1 and
	 i >= 1 and i < NX-1 ) {

      Anew(i,j) = -0.25 * (rhs(i,j) - ( A(i+1,j  ) +
					A(i-1,j  ) +
					A(i  ,j+1) +
					A(i  ,j-1) ));
      
      error = fmaxr( error, fabsr( Anew(i,j)-A(i,j) ) );

    }
  }

  // Tell each thread how to initialize its reduction result.
  KOKKOS_INLINE_FUNCTION
  void init (real_t& dst) const
  {
    // The identity under max is -Inf.
    // Kokkos does not come with a portable way to access
    // floating-point Inf and NaN. 
#ifdef __CUDA_ARCH__
    dst = -CUDART_INF;
#else
    dst = std::numeric_limits<real_t>::min();
#endif // __CUDA_ARCH__
  } // init

  // "Join" intermediate results from different threads.
  // This should normally implement the same reduction
  // operation as operator() above. Note that both input
  // arguments MUST be declared volatile.
  KOKKOS_INLINE_FUNCTION
  void join (volatile real_t& dst,
  	     const volatile real_t& src) const
  {
    // check if reduce value (dst) needs an update from src
    if (dst < src) {
      dst = src;
    }
  } // join

  DataArray A, Anew, rhs;
  Params params;
  
}; // end class Poisson2dKokkosKernel

/**
 * Compute a kokkos solution stored in Anew.
 */
void poisson2d_kokkos( DataContextKokkos& context, Params& params )
{
  int iter_max = params.iter_max;
  real_t tol = params.tol;
  
  int NX = params.NX;
  int NY = params.NY;
  
  int iter  = 0;
  real_t error = 1.0;

  // create a Kokkos Functor for Poisson computation, and launch computation
  Poisson2dKokkosKernel functor(context, params);
  
  while ( error > tol && iter < iter_max ) {
    error = 0.0;

    Kokkos::parallel_reduce(NX*NY, functor, error);

    // copy Anew into A
    Kokkos::deep_copy(context.A, context.Anew);
    
    // Ensure periodic boundary conditions
    Kokkos::parallel_for( NX, KOKKOS_LAMBDA(const int index) {    
    	int ix,iy;
    	index2coord(index,ix,iy,NX,NY);

    	if ( ix >= 1 and ix < NX-1 ) {
    	  context.A(ix,   0) = context.A(ix,NY-2);
    	  context.A(ix,NY-1) = context.A(ix,   1);
    	}
      });
    
    Kokkos::parallel_for( NY, KOKKOS_LAMBDA(const int index) {    
    	int ix,iy;
    	index2coord(index,ix,iy,NX,NY);

    	if ( iy >= 1 and iy < NY-1 ) {
    	  context.A(   0,iy) = context.A(NX-2,iy);
    	  context.A(NX-1,iy) = context.A(   1,iy);
    	}
      });

    if ( (iter % 100) == 0) printf("%5d, %0.6f\n", iter, error);
    iter++;

  } // end while

} // poisson2d_kokkos

#endif // LAPLACE2D_KOKKOS_KERNEL_H_
