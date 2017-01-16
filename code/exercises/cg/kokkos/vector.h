#pragma once

#include "kokkos_shared.h"

struct vector {

  unsigned int n;

  // DataArray coef
  typedef Kokkos::View<double*, DEVICE> DataArray;
  
  // host mirror
  typedef DataArray::HostMirror         DataArrayHost;

  DataArray coefs;

  // constructor
  vector(unsigned int _n) :
    n(_n),  coefs("coefs", _n) {}

  // init
  void init(double val) {
    
    // Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int i) { 
    // 	coefs(i) = val;
    //   });

    // alternate version: use deep_copy with a scalar value second as
    // second argument
    Kokkos::deep_copy(coefs, val);
    
  } // init
  
}; // struct vector


