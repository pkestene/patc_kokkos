// ===============================================================
// ===============================================================
// ===============================================================
/**
 * version 2:
 * - data is a 3d array
 * - only loops over i,j are parallelized, loop over k is kept inside
 * - optionally, on can use 1D subview to access data, and help the compiler
 *   to recognize a vectorizable loop
 *
 *
 * subview's are important here. 
 * Without 1d views, you'll have about 25% percent perfomance drop here.
 *
 * Use parameter use_1d_views to activate/deactivate the use of 1d views.
 */
double test_stencil_3d_flat_vector(int n, int nrepeat, bool use_1d_views) {

  uint64_t nbCells = n*n*n;

  uint64_t nbIter = /* TODO : how many iterations ? */;
  
  // Allocate Views
  DataArray x("X",n,n,n);
  DataArray y("Y",n,n,n);
  
  // Initialize arrays
  Kokkos::parallel_for(nbCells, KOKKOS_LAMBDA (const int& index) {
      int i,j,k;
      
      index2coord(index,i,j,k,n,n,n);
      
      x(i,j,k) = 1.0*(i+j+k+0.1);
      y(i,j,k) = 3.0*(i+j+k);
    });

  // Time computation
  Timer timer;
  
  timer.start();

  if (use_1d_views) {

    for(int irepeat = 0; irepeat < nrepeat; irepeat++) {
      
      // Do stencil
      Kokkos::parallel_for(nbIter, KOKKOS_LAMBDA (const int& index) {
	  int i,j;
	  
	  // index = j + n * i -- CPU
	  // index = i + n * j -- GPU
          // TODO : convert index to (i,j)

	  // TODO : create 1d views and uncomment the 6 following lines
	  //auto x_i_j = Kokkos::subview( /* TODO */);
	  //auto y_i_j = Kokkos::subview( /* TODO */);
	  
	  //auto x_im1_j = Kokkos::subview( /* TODO */);
	  //auto x_ip1_j = Kokkos::subview( /* TODO */);
	  
	  //auto x_i_jm1 = Kokkos::subview(x, i, j-1, Kokkos::ALL());
	  //auto x_i_jp1 = Kokkos::subview(x, i, j+1, Kokkos::ALL());
	  
	  /*
	   * subview's are important here. If you uncomment the following lines
	   * you'll about 25% percent perfomance drop here.
	   */
	  if (i>0 and i<n-1 and
	      j>0 and j<n-1) {
	    
	    // vectorization loop
#if defined( __INTEL_COMPILER )
#pragma ivdep
#endif	
	    for (int k=1; k<n-1; ++k) {
	      y_i_j(k) = -5*x_i_j(k) +
		( x_im1_j(k) + x_ip1_j(k) +
		  x_i_jm1(k) + x_i_jp1(k) +
		  x_i_j(k-1) + x_i_j(k+1) );
	    }
	  }
	  
	});
      
    } // end repeat

  } else { // don't use 1 d views
    
    for(int irepeat = 0; irepeat < nrepeat; irepeat++) {
      
      // Do stencil
      Kokkos::parallel_for(nbIter, KOKKOS_LAMBDA (const int& index) {
	  int i,j;
	  //index2coord(index,i,j,k,n,n,n);
	  
	  // index = j + n * i -- CPU
	  // index = i + n * j -- GPU
	  index2coord(index,i,j,n,n);
	  
	  if (i>0 and i<n-1 and
	      j>0 and j<n-1) {
	    
	    // vectorization loop
#if defined( __INTEL_COMPILER )
#pragma ivdep
#endif	
	    for (int k=1; k<n-1; ++k) {
	      y(i,j,k) = -5*x(i,j,k) +
	        ( x(i-1,j,k) + x(i+1,j,k) +
		  x(i,j-1,k) + x(i,j+1,k) +
		  x(i,j,k-1) + x(i,j,k+1) );
	    }
	  }
	  
	});
    } // end repeat
  }
  timer.stop();
  
  // Print results
  double time_seconds = timer.elapsed();

  // 6+1 reads + 1 write
  double dataSizeMBytes = 1.0e-6*nbCells*sizeof(real_t)*(7+1);
  double bandwidth = 1.0e-3*dataSizeMBytes*nrepeat/time_seconds;
  
  printf("#nbCells      Time(s)  TimePerIterations(s) size(MB) BW(GB/s)\n");
  printf("%13lu %8lf %20.3e  %6.3f %3.3f\n",
	 nbCells, time_seconds, time_seconds/nrepeat,
	 dataSizeMBytes,bandwidth);

  return bandwidth;
  
} // test_stencil_3d_flat_vector
