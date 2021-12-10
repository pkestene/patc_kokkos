#pragma once

#include<cstdlib>

#include "kokkos_shared.h"

struct matrix
{

  unsigned int N;
  unsigned int num_rows;
  unsigned int nnz;
  using Row_Offsets_Type = Kokkos::View<unsigned int *, Device>;
  using Cols_Type        = Kokkos::View<unsigned int *, Device>;
  using Coefs_Type       = Kokkos::View<double *,       Device>;

  Row_Offsets_Type row_offsets;
  Cols_Type        cols;
  Coefs_Type       coefs;

  // constructor
  matrix(int _N) :
    N(_N),
    num_rows( (N+1)*(N+1)*(N+1) ),
    nnz( 27*num_rows ),
    row_offsets("row_offsets", num_rows+1),
    cols("cols", nnz),
    coefs("coefs", nnz)
  {}

  // init
  void init()
  {

    int offsets_[27];
    double coefs_[27];
    int zstride=N*N;
    int ystride=N;

    int i=0;
    for(int z=-1;z<=1;z++)
    {
      for(int y=-1;y<=1;y++)
      {
	for(int x=-1;x<=1;x++)
        {
	  offsets_[i] = zstride*z+ystride*y+x;
	  if(x==0 && y==0 && z==0)
	    coefs_[i]=27;
	  else
	    coefs_[i]=-1;
	  i++;
	}
      }
    }

    // create host mirror
    Row_Offsets_Type::HostMirror h_row_offsets =
      Kokkos::create_mirror_view(row_offsets);

    Cols_Type::HostMirror h_cols =
      Kokkos::create_mirror_view(cols);

    Coefs_Type::HostMirror h_coefs =
      Kokkos::create_mirror_view(coefs);

    int innz=0;
    for(int i=0; i<num_rows; ++i)
    {
	h_row_offsets(i) = innz;
	for(int j=0;j<27;j++)
        {
	  int n = i+offsets_[j];
	  if(n>=0 && n<num_rows)
          {
	    h_cols(innz)  = n;
	    h_coefs(innz) = coefs_[j];
	    innz++;
	  }
	}

	if (i==0)
	  h_row_offsets(num_rows) = nnz;

    }

    // copy result on device
    Kokkos::deep_copy(row_offsets, h_row_offsets);
    Kokkos::deep_copy(cols,        h_cols);
    Kokkos::deep_copy(coefs,       h_coefs);

  } // init

}; // matrix
