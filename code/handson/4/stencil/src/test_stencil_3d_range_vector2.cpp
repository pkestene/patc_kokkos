// ===============================================================
// ===============================================================
// ===============================================================
/**
 * version 5 :
 * same as version 4 but uses Hierarchical parallelism, i.e.
 * - a TeamPolicy        Kokkos::policy for the outer loop 
 * - a ThreadVectorRange Kokkos::policy for the inner loop (to be vectorized)
 *
 */
double test_stencil_3d_range_vector2(int n, int nrepeat, int nbTeams) {

  uint64_t nbCells = n*n*n;
  
  // Allocate Views
  DataArray x("X",n,n,n);
  DataArray y("Y",n,n,n);

  // init 2d range policy
  using Range2D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> >;
  using Range3D = typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<3> >;

  Range2D range2d( {{0,0}}, {{n,n}} );
  Range3D range3d( {{0,0,0}}, {{n,n,n}} );

  
  // Initialize arrays using a 3d range policy
  Kokkos::parallel_for("init", range3d,
		       KOKKOS_LAMBDA (const int& i,
				      const int& j,
				      const int& k) {      
			 x(i,j,k) = 1.0*(i+j+k+0.1);
			 y(i,j,k) = 3.0*(i+j+k);
		       });

  // get prepared for TeamPolicy
  using team_member_t = typename Kokkos::TeamPolicy<>::member_type;

  
  // Time computation
  Timer timer;
  
  timer.start();
  for(int irepeat = 0; irepeat < nrepeat; irepeat++) {
    
    // Do stencil
    Kokkos::parallel_for
      ("stencil compute - 3d range vector",
       Kokkos::TeamPolicy<> (nbTeams, Kokkos::AUTO),
       KOKKOS_LAMBDA (team_member_t team_member) {

	int i = team_member.team_rank();
	int j = team_member.league_rank();

	int ni = (n+team_member.team_size()  -1)/team_member.team_size();
	int nj = (n+team_member.league_size()-1)/team_member.league_size();
	
	// make sure indexes (i,j) cover the full (nx,ny) range
	for (int jj=0, j=team_member.league_rank();
	     jj<nj;
	     ++jj, j+=team_member.league_size()) {
	  for (int ii=0, i=team_member.team_rank();
	       ii<ni;
	       ++ii, i+=team_member.team_size()) {

	    auto x_i_j = Kokkos::subview(x, i, j, Kokkos::ALL());
	    auto y_i_j = Kokkos::subview(y, i, j, Kokkos::ALL());
	    
	    auto x_im1_j = Kokkos::subview(x, i-1, j, Kokkos::ALL());
	    auto x_ip1_j = Kokkos::subview(x, i+1, j, Kokkos::ALL());
	    
	    auto x_i_jm1 = Kokkos::subview(x, i, j-1, Kokkos::ALL());
	    auto x_i_jp1 = Kokkos::subview(x, i, j+1, Kokkos::ALL());
	    
	    if (i>0 and i<n-1 and
		j>0 and j<n-1) {

	      // Kokkos::parallel_for
	      // 	(Kokkos::ThreadVectorRange(team_member, 1, n-1),
	      // 	 KOKKOS_LAMBDA(int& k) {

	      // 	  y_i_j(k) = -5*x_i_j(k) +
	      // 	  ( x_im1_j(k) + x_ip1_j(k) +
	      // 	    x_i_jm1(k) + x_i_jp1(k) +
	      // 	    x_i_j(k-1) + x_i_j(k+1) );
		  
	      // 	});
	      
	      // vectorization loop
#if defined( __INTEL_COMPILER )
#pragma ivdep
#endif
	      for (int k=1; k<n-1; ++k)
	      	y_i_j(k) = -5*x_i_j(k) +
	      	  ( x_im1_j(k) + x_ip1_j(k) +
	      	    x_i_jm1(k) + x_i_jp1(k) +
	      	    x_i_j(k-1) + x_i_j(k+1) );
		// y(i,j,k) = -5*x(i,j,k) +
		//   ( x(i-1,j,k) + x(i+1,j,k) +
		//     x(i,j-1,k) + x(i,j+1,k) +
		//     x(i,j,k-1) + x(i,j,k+1) );
	      
	    }
	  } // end for ii
	} // end for jj
	
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
  
} // test_stencil_3d_range_vector2
