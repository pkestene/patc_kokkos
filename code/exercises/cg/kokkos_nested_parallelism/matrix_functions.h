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

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t& thread) const
  {
    
    double sum=0;

    // get thread's team rank inside league
    const int i = thread.league_rank();
    
    int row_start = A.row_offsets(i);
    int row_end   = A.row_offsets(i+1);
    for(int j=row_start;j<row_end;j++) {
      unsigned int Acol = A.cols(j);
      double Acoef = A.coefs(j);
      double xcoef = x.coefs(Acol);
      sum += Acoef*xcoef;
    }
    y.coefs(i) = sum;

  }
  
  matrix A;
  vector x,y;

};

void matvec(matrix A, vector x, vector y, int num_teams) {


  // instantiate the functor
  MatvecFunctor functor(A, x, y);
  
  team_policy_t policy( num_teams, team_policy_t::team_size_max( functor ) );
  
  Kokkos::parallel_for(policy, functor);
  
} // matvec
