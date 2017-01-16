#pragma once

#include<cstdlib>
#include "vector.h"

double dot(vector x, vector y) {
  
  double sum=0;
  unsigned int n=x.n;

  Kokkos::parallel_reduce(n, KOKKOS_LAMBDA(const int i, double &local_sum) {
      local_sum += x.coefs(i) * y.coefs(i);
    }, sum);
  
  return sum;
  
} // dot

void waxpby(double alpha,
	    vector x,
	    double beta,
	    vector y,
	    vector w) {
  
  unsigned int n=x.n;

  Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int i) {
      
      w.coefs(i) = alpha * x.coefs(i) + beta * y.coefs(i);
      
    });

} // waxpy

