#ifndef DATA_CONTEXT_H_
#define DATA_CONTEXT_H_

#include "common.h"

struct DataContext {

  real_t* A;
  real_t* ANew;
  real_t* rhs;

  DataContext(Params& params) {
    int NX = params.NX;
    int NY = params.NY;
    
    A    = new real_t[NX*NY];
    ANew = new real_t[NX*NY];
    rhs  = new real_t[NX*NY];
  };

  ~DataContext() {

    delete[] A;
    delete[] ANew;
    delete[] rhs;
  };
  
};

#endif // DATA_CONTEXT_H_
