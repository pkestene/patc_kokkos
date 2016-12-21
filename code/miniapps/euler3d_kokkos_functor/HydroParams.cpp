#include "HydroParams.h"

#include <cstdlib> // for exit
#include <cstdio>  // for fprintf
#include <cstring> // for strcmp
#include <iostream>

#include "config/inih/ini.h" // our INI file reader

const char * varNames[5] = { "rho", "E", "mx", "my", "mz" };

/*
 * Hydro Parameters (read parameter file)
 */
void hydroParams_setup(HydroParams& params, ConfigMap &configMap)
{

  /* initialize RUN parameters */
  params.nStepmax = configMap.getInteger("run","nstepmax",1000);
  params.tEnd     = configMap.getFloat  ("run","tend",0.0);
  params.nOutput  = configMap.getInteger("run","noutput",100);
  if (params.nOutput == -1)
    params.enableOutput = false;
  
  /* initialize MESH parameters */
  params.nx = configMap.getInteger("mesh","nx", 2);
  params.ny = configMap.getInteger("mesh","ny", 2);
  params.nz = configMap.getInteger("mesh","nz", 2);

  params.xmin = configMap.getFloat("mesh", "xmin", 0.0);
  params.ymin = configMap.getFloat("mesh", "ymin", 0.0);
  params.zmin = configMap.getFloat("mesh", "zmin", 0.0);

  params.xmax = configMap.getFloat("mesh", "xmax", 1.0);
  params.ymax = configMap.getFloat("mesh", "ymax", 1.0);
  params.zmax = configMap.getFloat("mesh", "zmax", 1.0);

  params.boundary_type_xmin  = static_cast<int>(configMap.getInteger("mesh","boundary_type_xmin", BC_DIRICHLET));
  params.boundary_type_xmax  = static_cast<int>(configMap.getInteger("mesh","boundary_type_xmax", BC_DIRICHLET));
  params.boundary_type_ymin  = static_cast<int>(configMap.getInteger("mesh","boundary_type_ymin", BC_DIRICHLET));
  params.boundary_type_ymax  = static_cast<int>(configMap.getInteger("mesh","boundary_type_ymax", BC_DIRICHLET));
  params.boundary_type_zmin  = static_cast<int>(configMap.getInteger("mesh","boundary_type_zmin", BC_DIRICHLET));
  params.boundary_type_zmax  = static_cast<int>(configMap.getInteger("mesh","boundary_type_zmax", BC_DIRICHLET));

  params.settings.gamma0         = configMap.getFloat("hydro","gamma0", 1.4);
  params.settings.cfl            = configMap.getFloat("hydro", "cfl", 0.5);
  params.settings.iorder         = configMap.getInteger("hydro","iorder", 2);
  params.settings.slope_type     = configMap.getFloat("hydro","slope_type",1.0);
  params.settings.smallc         = configMap.getFloat("hydro","smallc", 1e-10);
  params.settings.smallr         = configMap.getFloat("hydro","smallr", 1e-10);

  params.niter_riemann  = configMap.getInteger("hydro","niter_riemann", 10);
  std::string riemannSolverStr = std::string(configMap.getString("hydro","riemann", "approx"));
  if ( !riemannSolverStr.compare("approx") ) {
    params.riemannSolverType = RIEMANN_APPROX;
  } else if ( !riemannSolverStr.compare("hll") ) {
    params.riemannSolverType = RIEMANN_HLL;
  } else if ( !riemannSolverStr.compare("hllc") ) {
    params.riemannSolverType = RIEMANN_HLLC;
  } else {
    std::cout << "Riemann Solver specified in parameter file is invalid\n";
    std::cout << "Use the default one : approx\n";
    params.riemannSolverType = RIEMANN_APPROX;
  }
    
  std::string problemStr = std::string(configMap.getString("hydro","problem", "unknown"));
  if ( !problemStr.compare("implode") ) {
    params.problemType = PROBLEM_IMPLODE;
  } else if ( !problemStr.compare("blast") ) {
    params.problemType = PROBLEM_BLAST;
  } else {
    std::cout << "Problem is invalid\n";
    std::cout << "Use the default one : implode\n";
    params.problemType = PROBLEM_IMPLODE;
  }

  params.blast_radius   = configMap.getFloat("blast","radius", (params.xmin+params.xmax)/2.0/10);
  params.blast_center_x = configMap.getFloat("blast","center_x", (params.xmin+params.xmax)/2);
  params.blast_center_y = configMap.getFloat("blast","center_y", (params.ymin+params.ymax)/2);
  params.blast_center_z = configMap.getFloat("blast","center_z", (params.zmin+params.zmax)/2);
  params.blast_density_in  = configMap.getFloat("blast","density_in", 1.0);
  params.blast_density_out = configMap.getFloat("blast","density_out", 1.2);
  params.blast_pressure_in  = configMap.getFloat("blast","pressure_in", 10.0);
  params.blast_pressure_out = configMap.getFloat("blast","pressure_out", 0.1);

  params.implementationVersion  = configMap.getFloat("OTHER","implementationVersion", 0);
  if (params.implementationVersion != 0 and
      params.implementationVersion != 1) {
    std::cout << "Implementation version is invalid (must be 0 or 1)\n";
    std::cout << "Use the default : 0\n";
    params.implementationVersion = 0;
  }

  hydroParams_init(params);

} // hydroParams_setup

// =======================================================
// =======================================================
void hydroParams_init(HydroParams& params)
{

  // set other parameters
  params.imax = params.nx - 1 + 2*params.ghostWidth;
  params.jmax = params.ny - 1 + 2*params.ghostWidth;
  params.kmax = params.nz - 1 + 2*params.ghostWidth;
  
  params.isize = params.imax - params.imin + 1;
  params.jsize = params.jmax - params.jmin + 1;
  params.ksize = params.kmax - params.kmin + 1;

  params.ijsize  = params.isize*params.jsize;
  params.ijksize = params.isize*params.jsize*params.ksize;
  
  params.dx = (params.xmax - params.xmin) / params.nx;
  params.dy = (params.ymax - params.ymin) / params.ny;
  params.dz = (params.zmax - params.zmin) / params.nz;
  
  params.settings.smallp  = params.settings.smallc*params.settings.smallc/
    params.settings.gamma0;
  params.settings.smallpp = params.settings.smallr*params.settings.smallp;
  params.settings.gamma6  = (params.settings.gamma0 + ONE_F)/(TWO_F * params.settings.gamma0);
  
  // check that given parameters are valid
  if ( (params.implementationVersion != 0) && 
       (params.implementationVersion != 1) && 
       (params.implementationVersion != 2) ) {
    fprintf(stderr, "The implementation version parameter should 0,1 or 2 !!!");
    fprintf(stderr, "Check your parameter file, section OTHER");
    exit(EXIT_FAILURE);
  } else {
    fprintf(stdout, "Using Euler implementation version %d\n", params.implementationVersion);
  }

} // hydroParams_init


// =======================================================
// =======================================================
void hydroParams_print(HydroParams& params)
{
  
  printf( "##########################\n");
  printf( "Simulation run parameters:\n");
  printf( "##########################\n");
  printf( "nx         : %d\n", params.nx);
  printf( "ny         : %d\n", params.ny);
  printf( "nz         : %d\n", params.nz);
  
  printf( "dx         : %f\n", params.dx);
  printf( "dy         : %f\n", params.dy);
  printf( "dz         : %f\n", params.dz);

  printf( "imin       : %d\n", params.imin);
  printf( "imax       : %d\n", params.imax);

  printf( "jmin       : %d\n", params.jmin);      
  printf( "jmax       : %d\n", params.jmax);      

  printf( "kmin       : %d\n", params.kmin);      
  printf( "kmax       : %d\n", params.kmax);      

  printf( "nStepmax   : %d\n", params.nStepmax);
  printf( "tEnd       : %f\n", params.tEnd);
  printf( "nOutput    : %d\n", params.nOutput);
  printf( "gamma0     : %f\n", params.settings.gamma0);
  printf( "cfl        : %f\n", params.settings.cfl);
  printf( "smallr     : %12.10f\n", params.settings.smallr);
  printf( "smallc     : %12.10f\n", params.settings.smallc);
  //printf( "niter_riemann : %d\n", params.niter_riemann);
  printf( "iorder     : %d\n", params.settings.iorder);
  printf( "slope_type : %f\n", params.settings.slope_type);
  printf( "riemann    : %d\n", params.riemannSolverType);
  printf( "problem    : %d\n", params.problemType);
  printf( "implementation version : %d\n",params.implementationVersion);
  printf( "##########################\n");
  
} // hydroParams_print
