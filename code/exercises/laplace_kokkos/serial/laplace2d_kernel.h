#ifndef LAPLACE2D_KERNEL_H_
#define LAPLACE2D_KERNEL_H_

#include "params.h"
#include "DataContext.h"

void poisson2d_serial( DataContext& context, Params& params )
{
  int iter_max = params.iter_max;
  real_t tol = params.tol;
  
  int NX = params.NX;
  int NY = params.NY;

  real_t* A    = context.A;
  real_t* ANew = context.ANew;
  real_t* rhs  = context.rhs;
  
  int iter  = 0;
  real error = 1.0;

  while ( error > tol && iter < iter_max )
    {
      error = 0.0;
      
      for( int iy = 1; iy < NY-1; iy++)
        {
	  for( int ix = 1; ix < NX-1; ix++ )
            {
	      Anew[iy*NX+ix] = -0.25 * (rhs[iy*NX+ix] - ( Aref[ iy   *NX+ix+1] +
							  Aref[ iy   *NX+ix-1] +
							  Aref[(iy-1)*NX+ix] +
							  Aref[(iy+1)*NX+ix] ));
	      //error = fmaxr( error, fabsr(Anew[iy*NX+ix]-Aref[iy*NX+ix]));
            }
        }
      
      
#endif // LAPLACE2D_KERNEL_H_
