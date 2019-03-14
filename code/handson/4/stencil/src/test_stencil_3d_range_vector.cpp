// ===============================================================
// ===============================================================
// ===============================================================
/**
 * version 4 :
 * same as version 3 but uses a 2d Range policy and keep the loop over
 * index k inside kernel.
 */
double test_stencil_3d_range_vector(int n, int nrepeat) {

  uint64_t nbCells = n*n*n;
  
  // Allocate Views
  DataArray x("X",n,n,n);
  DataArray y("Y",n,n,n);

  // init 2d range policy
  using Range2D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> >;
  using Range3D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<3> >;

  Range2D range2d( {{0,0}}, {{n,n}} );
  Range3D range3d( {{0,0,0}}, {{n,n,n}} );

  
  // Initialize arrays
  Kokkos::parallel_for("init", range3d,
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
      ("stencil compute - 3d range vector", range2d, KOKKOS_LAMBDA (const int& i,
								    const int& j) {

	auto x_i_j = Kokkos::subview(x, i, j, Kokkos::ALL());
	auto y_i_j = Kokkos::subview(y, i, j, Kokkos::ALL());

	auto x_im1_j = Kokkos::subview(x, i-1, j, Kokkos::ALL());
	auto x_ip1_j = Kokkos::subview(x, i+1, j, Kokkos::ALL());

	auto x_i_jm1 = Kokkos::subview(x, i, j-1, Kokkos::ALL());
	auto x_i_jp1 = Kokkos::subview(x, i, j+1, Kokkos::ALL());

	if (i>0 and i<n-1 and
	    j>0 and j<n-1) {

	  // vectorization loop
#if defined( __INTEL_COMPILER )
#pragma ivdep
#endif
	  for (int k=1; k<n-1; ++k)
	    y_i_j(k) = -5*x_i_j(k) +
	      ( x_im1_j(k) + x_ip1_j(k) +
		x_i_jm1(k) + x_i_jp1(k) +
		x_i_j(k-1) + x_i_j(k+1) );
	}
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
  
} // test_stencil_3d_range_vector
