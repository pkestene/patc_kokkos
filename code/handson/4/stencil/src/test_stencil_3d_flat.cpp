// ===============================================================
// ===============================================================
// ===============================================================
/**
 * version 1 : naive
 * - data is a 3d array
 * - all 3 loops parallelized with a single parallel_for
 *
 * See the slides.
 *
 * \return effective bandwidth
 */
double test_stencil_3d_flat(int n, int nrepeat) {

  uint64_t nbCells = n*n*n;
  
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
  for(int irepeat = 0; irepeat < nrepeat; irepeat++) {
    
    // Do stencil
    Kokkos::parallel_for(nbCells, KOKKOS_LAMBDA (const int& index) {
	int i,j,k;

        // TODO convert index to (i,j,k)

	if ( /* TODO */ )
	y(i,j,k) = -5*x(i,j,k) +
	  ( x(i-1,j  ,k  ) + x(i+1,j  ,k  ) +
	    x(i  ,j-1,k  ) + x(i  ,j+1,k  ) +
	    x(i  ,j  ,k-1) + x(i  ,j  ,k+1) );
      });
    
  }
  timer.stop();
  
  // Print results
  double time_seconds = timer.elapsed();

  // compute a sort of bandwidth : 6+1 reads + 1 write
  double dataSizeMBytes = 1.0e-6*nbCells*sizeof(real_t)*(7+1);
  double bandwidth = 1.0e-3*dataSizeMBytes*nrepeat/time_seconds;
  
  printf("#nbCells      Time(s)  TimePerIterations(s) size(MB) BW(GB/s)\n");
  printf("%13lu %8lf %20.3e  %6.3f %3.3f\n",
	 nbCells, time_seconds, time_seconds/nrepeat,
	 dataSizeMBytes,bandwidth);

  return bandwidth;
  
} // test_stencil_3d_flat
