#pragma once

#include "vector.h"
#include "matrix.h"

void matvec(matrix A, vector x, vector y) {

  Kokkos::parallel_for(A.num_rows, KOKKOS_LAMBDA(const int i) {

      double sum=0;
      int row_start = A.row_offsets(i);
      int row_end   = A.row_offsets(i+1);
      for(int j=row_start;j<row_end;j++) {
	unsigned int Acol = A.cols(j);
	double Acoef = A.coefs(j);
	double xcoef = x.coefs(Acol);
	sum += Acoef*xcoef;
      }
      y.coefs(i) = sum;
      
    });
  
} // matvec
