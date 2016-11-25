#include <cstdio>
#include <iostream>
#include <fstream>

#include <omp.h>

#include "constants.h"
#include "kokkos_shared.h"
#include "mandelbrot.h"

#ifdef CUDA
#include "CudaTimer.h"
#else // OpenMP
#include "OpenMPTimer.h"
#endif

#include <unistd.h>

using namespace std;

int main(int argc, char* argv[]) {

  /*
   * Initialize kokkos (host + device)
   */
#ifdef CUDA
  // Initialize Host mirror device
  Kokkos::HostSpace::execution_space::initialize();
  //const unsigned device_count = Kokkos::Cuda::detect_device_count();

  // Use the first device:
  Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0));
#else // OpenMP CPU
  Kokkos::initialize(argc, argv);
#endif

  {
    std::cout << "##########################\n";
    std::cout << "KOKKOS CONFIG             \n";
    std::cout << "##########################\n";
    
    std::ostringstream msg;
    std::cout << "Kokkos configuration" << std::endl;
    if ( Kokkos::hwloc::available() ) {
      msg << "hwloc( NUMA[" << Kokkos::hwloc::get_available_numa_count()
          << "] x CORE["    << Kokkos::hwloc::get_available_cores_per_numa()
          << "] x HT["      << Kokkos::hwloc::get_available_threads_per_core()
          << "] )"
          << std::endl ;
    }
#if defined( CUDA )
    Kokkos::Cuda::print_configuration( msg );
#else
    Kokkos::OpenMP::print_configuration( msg );
#endif
    std::cout << msg.str();
    std::cout << "##########################\n";
  }

#ifdef CUDA
  CudaTimer timer;
#else
  OpenMPTimer timer;
#endif

  Constants constants = Constants();

  // do computation by block
  int num_blocks = 8;
  int block_size = (constants.HEIGHT/num_blocks)*constants.WIDTH;

  // prepare data array for Mandelbrot set computation
  DataArray     image     = DataArray("image", constants.WIDTH*constants.HEIGHT);
  DataArrayHost imageHost = DataArrayHost("imageHost", constants.WIDTH*constants.HEIGHT);
  
  /*
   * Actual computation :
   * - loop on block in OpenMP
   * - each OpenMP threads launches a compution on GPU.
   */
  timer.start();
    
  // compute Mandelbrot set on block and copy back result using OpenMP policy
  MandelbrotByBlockAndCopyFunctor functor(image, imageHost, constants, block_size);
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0,num_blocks), functor);

  /*
   * the following is not safe, because computation of Mandelbrot sets and
   * copy of result on host are done in different cuda streams, without
   * ordering waranty.
   */
  // MandelbrotByBlockFunctor functor1(image, constants, block_size);
  // Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0,num_blocks), functor1);

  // MandelbrotCopyFunctor    functor2(image, imageHost, block_size);
  // Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::OpenMP>(0,num_blocks), functor2);

  timer.stop();
  printf("Time: %lf seconds.\n", timer.elapsed());
  
  // print aesthetically, dont read this part
  int xmax=80;
  int ymax=60;
  for(int y=0;y<ymax;y++) {
    printf("\n");
    for(int x=0;x<xmax;x++) {
      int index = (y*constants.HEIGHT/ymax)*constants.WIDTH+(x*constants.WIDTH/xmax);
      int val = imageHost(index);
      
      if (val==200) printf("&");
      else if (val==42) printf("X");
      else if(val>64) printf("#");
      else if(val>32) printf(":");
      else if(val>8) printf(".");
      else printf(" ");
    }
  }
  
  printf("\n");

  // save color ppm file
  if (1) {
    FILE* myfile = fopen("mandelbrot.ppm","w");
    
    fprintf(myfile, "P6 %d %d 255\n", constants.WIDTH , constants.HEIGHT);
    for(unsigned int index=0;index<constants.WIDTH*constants.HEIGHT; ++index) {
      unsigned char data;
      // create an arbitrary RBG code mapping values taken by imageHost
      data = imageHost(index) % 4 * 64;
      fwrite(&data,1,1,myfile);
      data = imageHost(index) % 8 * 32;
      fwrite(&data,1,1,myfile);
      data = imageHost(index) % 16 * 16;
      fwrite(&data,1,1,myfile);
    }

    fclose(myfile);
  }
  
   
#ifdef CUDA
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
#else
  Kokkos::finalize();
#endif

  return 0;
}
