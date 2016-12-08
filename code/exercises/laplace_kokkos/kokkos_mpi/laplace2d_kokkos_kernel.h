#ifndef LAPLACE2D_KOKKOS_KERNEL_H_
#define LAPLACE2D_KOKKOS_KERNEL_H_

#include "params.h"
#include "DataContextKokkos.h"

#include <mpi.h>

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
  {
    int i,j;
    index2coord(index,i,j,params.NX,params.NY);
    j += params.iy_start;

    const int NX = params.NX;
    const int NY = params.NY;

    const int ix_start = params.ix_start;
    const int ix_end   = params.ix_end;
    const int iy_start = params.iy_start;
    const int iy_end   = params.iy_end;
    
    if ( j >= iy_start and j < iy_end and
	 i >= ix_start and i < ix_end ) {
      
      Anew(index) = -0.25 * (rhs(index) - ( A( coord2index(i+1,j  ,NX,NY) ) +
					    A( coord2index(i-1,j  ,NX,NY) ) +
					    A( coord2index(i  ,j+1,NX,NY) ) +
					    A( coord2index(i  ,j-1,NX,NY) ) ));
      
      error = fmaxr( error, fabsr( Anew(index)-A(index) ) );

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

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  int mpi_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);

  int ix_start = params.ix_start;
  int ix_end   = params.ix_end;
  int iy_start = params.iy_start;
  int iy_end   = params.iy_end;
  
  // Ensure correctness if NY%nRanks != 0
  int NY_chunk = iy_end - iy_start;

  while ( error > tol && iter < iter_max ) {
    error = 0.0;

    // create a Kokkos Functor for Poisson computation, and launch computation
    Poisson2dKokkosKernel functor(context, params);
    Kokkos::parallel_reduce(NX*NY_chunk, functor, error);

    // MPI reduce to compute global error
    real_t globalerror = 0.0;
    MPI_Allreduce( &error, &globalerror, 1, MPI_REAL_TYPE, MPI_MAX, MPI_COMM_WORLD );
    error = globalerror;
    
    // copy Anew into A for current chunk
    //Kokkos::deep_copy(context.A, context.Anew);
    
    Kokkos::parallel_for( NX*NY_chunk, KOKKOS_LAMBDA(const int index) {    
    	int ix,iy;
    	index2coord(index,ix,iy,NX,NY);
	iy += iy_start;
	
    	if ( ix >= ix_start and ix < ix_end and
	     iy >= iy_start and iy < iy_end) {
    	  context.A(coord2index(ix,iy,NX,NY)) = context.Anew(coord2index(ix,iy,NX,NY));
    	}
      });


    // determine with whom we will have to exchange data
    int top    = (mpi_rank == 0)             ? (mpi_ranks-1) : mpi_rank-1;
    int bottom = (mpi_rank == (mpi_ranks-1)) ? 0             : mpi_rank+1;

    // assumes a CUDA-aware MPI implementation 
    {
      //1. Sent row iy_start (first modified row) to top
      //   receive lower boundary (iy_end) from bottom

      // create subviews
      typedef DataArray::size_type size_type;
      typedef Kokkos::pair<size_type, size_type> pair;
      DataArray rowBlock = Kokkos::subview(image, pair(start,end));


      MPI_Sendrecv( &A[iy_start][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, top   , 0, &A[iy_end][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      
      //2. Sent row (iy_end-1) (last modified row) to bottom receive upper boundary (iy_start-1) from top
      MPI_Sendrecv( &A[(iy_end-1)][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, bottom, 0, &A[(iy_start-1)][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, top   , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
    }



    
    Kokkos::parallel_for( NY, KOKKOS_LAMBDA(const int index) {    
	int ix,iy;
	index2coord(index,ix,iy,NX,NY);

	if ( iy >= 1 and iy < NY-1 ) {
	  context.A(coord2index(   0,iy,NX,NY)) = context.A(coord2index(NX-2,iy,NX,NY));
	  context.A(coord2index(NX-1,iy,NX,NY)) = context.A(coord2index(   1,iy,NX,NY));
	}
      });

    if (mpi_rank == 0) if ( (iter % 100) == 0) printf("%5d, %0.6f\n", iter, error);
    iter++;

  } // end while

} // poisson2d_kokkos

#endif // LAPLACE2D_KOKKOS_KERNEL_H_
