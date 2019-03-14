// ===============================================================
// ===============================================================
// ===============================================================
/**
 * version 3 :
 * same as version 1 (naive) but uses a 3d range policy.
 */
double test_stencil_3d_range(int n, int nrepeat) {

  uint64_t nbCells = n*n*n;
  
  // Allocate Views
  DataArray x("X",n,n,n);
  DataArray y("Y",n,n,n);

  // init 3d range policy
  using Range3D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<3> >;

  Range3D range( /* TODO */ );

  
  // Initialize arrays
  Kokkos::parallel_for("init", range,
		       KOKKOS_LAMBDA (const int& i,
				      const int& j,
				      const int& k) {      
			 x(i,j,k) = 1.0*(i+j+k+0.1);
			 y(i,j,k) = 3.0*(i+j+k);
		       });

  // Time computation
  Timer timer;
  
  timer.start();
  for(int irepeat = 0; irepeat < nrepeat; irepeat++) {
    
    // Do stencil
    Kokkos::parallel_for
      ("stencil compute - 3d range", range, KOKKOS_LAMBDA (const int& i,
							   const int& j,
							   const int& k) {
	if (i>0 and i<n-1 and
	    j>0 and j<n-1 and
	    k>0 and k<n-1 )
	  y(i,j,k) = -5*x(i,j,k) +
	    ( x(i-1,j,k) + x(i+1,j,k) +
	      x(i,j-1,k) + x(i,j+1,k) +
	      x(i,j,k-1) + x(i,j,k+1) );
      });
    
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
  
} // test_stencil_3d_range
