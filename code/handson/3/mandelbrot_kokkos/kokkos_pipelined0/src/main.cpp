#include <cstdio>
#include <iostream>
#include <fstream>

#include <omp.h>

#include "constants.h"
#include "kokkos_shared.h"
#include "mandelbrot.h"

#ifdef KOKKOS_ENABLE_CUDA
#include "CudaTimer.h"
#else // OpenMP
#include "OpenMPTimer.h"
#endif

#include <unistd.h>

using namespace std;

// ============================================================
// ============================================================
void compute_mandelbrot(int argc, char* argv[])
{
#ifdef KOKKOS_ENABLE_CUDA
  CudaTimer timer;
#else
  OpenMPTimer timer;
#endif
  
  int default_size = 1024;
  if (argc>1)
    default_size = std::atoi(argv[1]);
  Constants constants = Constants(default_size);

  // do computation by block
  int num_blocks = 8;
  int block_size = (constants.HEIGHT/num_blocks)*constants.WIDTH;

  // prepare data array for Mandelbrot set computation
  DataArray     image     = DataArray("image", constants.WIDTH*constants.HEIGHT);
  DataArrayHost imageHost = Kokkos::create_mirror_view(image);

  /*
   * Actual computation (GPU with CUDA or CPU with OpenMP)
   */
  timer.start();
  
  for (int block = 0; block < num_blocks; ++block) {
    int start = block * block_size;
    int end   = start + block_size;

    // create subviews
    typedef DataArray::size_type size_type;
    typedef Kokkos::pair<size_type, size_type> pair;
    DataArray     imageBlock     = Kokkos::subview(image,
						   pair(start,end));
    
    DataArrayHost imageBlockHost = Kokkos::subview(imageHost,
						   pair(start,end));

    // compute Mandelbrot set on block
    MandelbrotFunctor functor(image, constants);
    Kokkos::parallel_for(Kokkos::RangePolicy<Device>(start,end), functor);
        
    // copy back results from device to host
    Kokkos::deep_copy(imageBlockHost,imageBlock);
    
  }


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

} // compute_mandelbrot

int main(int argc, char* argv[]) {

  /*
   * Initialize kokkos (host + device)
   */
  Kokkos::initialize(argc, argv);

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
    Kokkos::print_configuration( msg );
    std::cout << msg.str();
    std::cout << "##########################\n";
  }

  compute_mandelbrot(argc,argv);
   
  Kokkos::finalize();

  return 0;
}
