#ifndef DATA_CONTEXT_KOKKOS_H_
#define DATA_CONTEXT_KOKKOS_H_

#include "common.h"
#include "kokkos_shared.h"

struct DataContextKokkos {

  DataArray A;
  DataArray Anew;
  DataArray rhs;

  DataContextKokkos(Params& params) :
    A   ("A",    params.NX, params.NY),
    Anew("Anew", params.NX, params.NY),
    rhs ("rhs",  params.NX, params.NY)
  {};

  ~DataContextKokkos() {}
  
};

#endif // DATA_CONTEXT_KOKKOS_H_
