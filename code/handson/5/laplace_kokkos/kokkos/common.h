#ifndef COMMON_H_
#define COMMON_H_

#ifdef USE_DOUBLE
    typedef double real_t;
    #define fmaxr fmax
    #define fabsr fabs
    #define expr exp
    #define MPI_REAL_TYPE MPI_DOUBLE
#else
    typedef float real_t;
    #define fmaxr fmaxf
    #define fabsr fabsf
    #define expr expf
    #define MPI_REAL_TYPE MPI_FLOAT
#endif

#endif // COMMON_H_
