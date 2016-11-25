#ifndef MANDELBROT_H_
#define MANDELBROT_H_

#include "constants.h"
#include "kokkos_shared.h"

/**
 * Kokkos kernel functor to compute mandelbrot set.
 */
class MandelbrotFunctor {

public:
  MandelbrotFunctor(DataArray image, Constants constants):
    image(image), 
    xmin(constants.xmin),
    ymin(constants.ymin),
    dx(constants.dx),
    dy(constants.dy),
    WIDTH(constants.WIDTH),
    HEIGHT(constants.HEIGHT),
    MAX_COLOR(constants.MAX_COLOR),
    MAX_ITERS(constants.MAX_ITERS)    
  {}

  KOKKOS_INLINE_FUNCTION
  unsigned char mandelbrot(int Px, int Py) const
  {

    float x0=xmin+Px*dx;
    float y0=ymin+Py*dy;
    float x=0.0;
    float y=0.0;
    int i;
    for(i=0;x*x+y*y<4.0 && i<MAX_ITERS;i++) {
      float xtemp=x*x-y*y+x0;
      y=2*x*y+y0;
      x=xtemp;
    }
    return (float) MAX_COLOR*i/MAX_ITERS;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index) const
  {
    int i,j;
    index2coord(index,i,j,WIDTH,HEIGHT);

    image(index)=mandelbrot(i,j);

  }

  DataArray image;
  float xmin, ymin, dx, dy, WIDTH, HEIGHT, MAX_COLOR, MAX_ITERS;

}; // end class MandelBrotFunctor

/**
 * Kokkos kernel functor to compute mandelbrot set.
 */
class MandelbrotByBlockFunctor {

  public:
  MandelbrotByBlockFunctor(DataArray     image,
			   Constants     constants,
			   int           block_size) :
    image(image),
    constants(constants),
    block_size(block_size)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& block) const
  {
      
    int start = block * block_size;
    int end   = start + block_size;

    // compute Mandelbrot set on block
    MandelbrotFunctor functor(image, constants);
    Kokkos::parallel_for(Kokkos::RangePolicy<DEVICE>(start,end), functor);     

  }

  DataArray     image;
  Constants     constants;
  int           block_size;

}; // class MandelbrotByBlockFunctor

/**
 * Kokkos kernel functor to compute mandelbrot set.
 */
class MandelbrotCopyFunctor {

  typedef DataArray::size_type size_type;
  typedef Kokkos::pair<size_type, size_type> pair;

public:
  MandelbrotCopyFunctor(DataArray     image,
			DataArrayHost imageHost,
			int           block_size) :
    image(image),
    imageHost(imageHost),
    block_size(block_size)
  {}
    
  KOKKOS_INLINE_FUNCTION
  void operator()(const int& block) const
  {

    int start = block * block_size;
    int end   = start + block_size;

    // create subviews
    DataArray     imageBlock     = Kokkos::subview(image,
						      pair(start,end));
    
    DataArrayHost imageBlockHost = Kokkos::subview(imageHost,
						       pair(start,end));

    // copy back results from device to host
    Kokkos::deep_copy(Kokkos::OpenMP(),imageBlockHost,imageBlock);
    
  }
    
  DataArray     image;
  DataArrayHost imageHost;
  int           block_size;

}; // class MandelbrotCopyFunctor

/**
 * Kokkos kernel functor to compute mandelbrot set.
 */
class MandelbrotByBlockAndCopyFunctor {

  typedef DataArray::size_type size_type;
  typedef Kokkos::pair<size_type, size_type> pair;

public:
  MandelbrotByBlockAndCopyFunctor(DataArray     image,
				  DataArrayHost imageHost,
				  Constants     constants,
				  int           block_size) :
    image(image),
    imageHost(imageHost),
    constants(constants),
    block_size(block_size)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& block) const
  {
      
    int start = block * block_size;
    int end   = start + block_size;

    // compute Mandelbrot set on block
    MandelbrotFunctor functor(image, constants);
    Kokkos::parallel_for(Kokkos::RangePolicy<DEVICE>(start,end), functor);     
    
    // create subviews
    DataArray     imageBlock     = Kokkos::subview(image,
						   pair(start,end));
    
    DataArrayHost imageBlockHost = Kokkos::subview(imageHost,
						   pair(start,end));
    
    // copy back results from device to host
    //Kokkos::deep_copy(Kokkos::OpenMP(),imageBlockHost,imageBlock);
    Kokkos::deep_copy(imageBlockHost,imageBlock);

  }

  DataArray     image;
  DataArrayHost imageHost;
  Constants     constants;
  int           block_size;

}; // class MandelbrotByBlockAndCopyFunctor

#endif

