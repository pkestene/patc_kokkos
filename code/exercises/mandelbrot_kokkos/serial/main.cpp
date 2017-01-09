#include <cstdio>
#include <iostream>
#include <fstream>

#include <vector>
#include <omp.h>

#include "constants.h"

#include "OpenMPTimer.h"

#include <unistd.h>

using namespace std;

/* some global constants */
unsigned int WIDTH;
unsigned int HEIGHT;
unsigned int MAX_ITERS;
unsigned int MAX_COLOR;
double       xmin;
double       xmax;
double       ymin;
double       ymax;
double         dx;
double         dy;

/* 
 * compute Mandelbrot pixelwise.
 * some global variables are used as parameters
 */
unsigned char mandelbrot(int Px, int Py)
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

/* convert a linear index to cartesian coordinate */
inline
void index2coord(int index, int &i, int &j, int Nx, int Ny) {
  i = index / Ny;
  j = index - i*Ny;
}

int main(int argc, char* argv[]) {

  OpenMPTimer timer;

  Constants constants = Constants();

  /* use constants to initialize the global constants */
  WIDTH     = constants.WIDTH;
  HEIGHT    = constants.HEIGHT;
  MAX_ITERS = constants.MAX_ITERS;
  MAX_COLOR = constants.MAX_COLOR;
  xmin      = constants.xmin;
  xmax      = constants.xmax;
  ymin      = constants.ymin;
  ymax      = constants.ymax;
  dx        = constants.dx;
  dy        = constants.dy;

  
  // Allocate data array for Mandelbrot set computation
  std::vector<unsigned char> image(constants.WIDTH*constants.HEIGHT);

  /*
   * Actual computation
   */
  timer.start();

  for(int index=0; index<WIDTH*HEIGHT; ++index) {
    
    int i,j;
    index2coord(index,i,j,WIDTH,HEIGHT);
    
    image[index]=mandelbrot(i,j);
    
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
      int index = (y*constants.HEIGHT/ymax)*constants.WIDTH+(x*constants.WIDTH/xmax);
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
    
    fprintf(myfile, "P6 %d %d 255\n", constants.WIDTH , constants.HEIGHT);
    for(unsigned int index=0;index<constants.WIDTH*constants.HEIGHT; ++index) {
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
