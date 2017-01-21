#pragma once

#include "vector.h"
#include "matrix.h"

// team policy typedef
typedef Kokkos::TeamPolicy<>       team_policy_t;
typedef team_policy_t::member_type team_member_t ;

class MatvecFunctor {
 
public:
  MatvecFunctor(matrix A, vector x, vector y):
    A(A), x(x), y(y)
  {}

  // this is executed by all threads
  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t& thread) const
  {

    double sum=0;

    // get thread's team rank inside league
    const int i = thread.league_rank();
    
    int row_start = A.row_offsets(i);
    int row_end   = A.row_offsets(i+1);

    // split the integer range [row_start, row_end[
    // among all the threads in the team
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(thread,row_start, row_end), [&] (const int& j, double & local_sum) {
	unsigned int Acol = A.cols(j);
	double Acoef = A.coefs(j);
	double xcoef = x.coefs(Acol);
	local_sum += Acoef*xcoef;
      }, sum);

    // only 1st thread write the team's result
    if (thread.team_rank() == 0)
      y.coefs(i) = sum;

  }
  
  matrix A;
  vector x,y;

};

void matvec(matrix A, vector x, vector y) {


  // instantiate the functor
  MatvecFunctor functor(A, x, y);

  // number of threads teams is A.num_rows (the outer loop bound)
  int num_teams = A.num_rows;
  
  // the total number of threads will be
  // num_teams * team_policy_t::team_size_max( functor )
  team_policy_t policy( num_teams, team_policy_t::team_size_max( functor ) );

  // run the parallel loop
  Kokkos::parallel_for(policy, functor);
  
} // matvec
