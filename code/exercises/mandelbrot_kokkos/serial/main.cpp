#include <cstdio>
#include <iostream>
#include <fstream>

#include <vector>
#include <omp.h>

#include "constants.h"

#include "SimpleTimer.h"

#include <unistd.h>

using namespace std;

/* 
 * compute Mandelbrot pixelwise.
 * some global variables are used as parameters
 */
unsigned char mandelbrot(int Px, int Py, const Constants& c)
{
  
  float x0=c.xmin+Px*c.dx;
  float y0=c.ymin+Py*c.dy;
  float x=0.0;
  float y=0.0;
  int i;
  for(i=0;x*x+y*y<4.0 && i<c.MAX_ITERS;i++) {
    float xtemp=x*x-y*y+x0;
    y=2*x*y+y0;
    x=xtemp;
  }
  return (float) c.MAX_COLOR*i/c.MAX_ITERS;
}

/* convert a linear index to cartesian coordinate */
inline
void index2coord(int index, int &i, int &j, int Nx, int Ny) {
  i = index / Ny;
  j = index - i*Ny;
}

int main(int argc, char* argv[]) {
  
  int default_size = 1024;
  if (argc>1)
    default_size = std::atoi(argv[1]);
  
  Constants c = Constants(default_size);
  
  // Allocate data array for Mandelbrot set computation
  std::vector<unsigned char> image(c.WIDTH*c.HEIGHT);

  /*
   * Actual computation
   */
  SimpleTimer timer;
  timer.start();

  for(int index=0; index<c.WIDTH*c.HEIGHT; ++index) {
    
    int i,j;
    index2coord(index,i,j,c.WIDTH,c.HEIGHT);
    
    image[index]=mandelbrot(i,j,c);
    
  }
  
  printf("end of loop reached ...\n");

  timer.stop();
  printf("Time: %lf seconds.\n", timer.elapsed());
  
  // print aesthetically, dont read this part
  int xmax=80;
  int ymax=60;
  for(int y=0;y<ymax;y++) {
    printf("\n");
    for(int x=0;x<xmax;x++) {
      int index = (y*c.HEIGHT/ymax)*c.WIDTH+(x*c.WIDTH/xmax);
      int val = image[index];
      
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
  if (0) {
    FILE* myfile = fopen("mandelbrot.ppm","w");
    
    fprintf(myfile, "P6 %d %d 255\n", c.WIDTH , c.HEIGHT);
    for(unsigned int index=0;index<c.WIDTH*c.HEIGHT; ++index) {
      unsigned char data;
      // create an arbitrary RBG code mapping values taken by image
      data = image[index] % 4 * 64;
      fwrite(&data,1,1,myfile);
      data = image[index] % 8 * 32;
      fwrite(&data,1,1,myfile);
      data = image[index] % 16 * 16;
      fwrite(&data,1,1,myfile);
    }

    fclose(myfile);
  }
  
  return 0;
}
