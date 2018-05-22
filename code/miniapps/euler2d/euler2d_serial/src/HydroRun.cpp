#include <string> 
#include <cstdio>
#include <cstdbool>
#include <sstream>
#include <fstream>

#include "HydroRun.h"
#include "HydroParams.h"
#include "Timer.h"

static bool isBigEndian()
{
  const int i = 1;
  return ( (*(char*)&i) == 0 );
}

void swapValues(real_t *a, real_t *b) {

  real_t tmp = *a;

  *a = *b;
  *b = tmp;

} // swapValues


// =======================================================
// =======================================================
/**
 *
 */
HydroRun::HydroRun(HydroParams& params, ConfigMap& configMap) :
  params(params),
  configMap(configMap),
  U(), U2(), Q(),
  Fluxes_x(), Fluxes_y(),
  Slopes_x(), Slopes_y(),
  Qm_x(), Qm_y(), Qp_x(), Qp_y()
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  
  /*
   * memory allocation (use sizes with ghosts included)
   */
  U.allocate(make_uint3(isize,jsize,NBVAR));
  U2.allocate(make_uint3(isize,jsize,NBVAR));
  Q.allocate(make_uint3(isize,jsize,NBVAR));

  if (params.implementationVersion == 0) {

    Fluxes_x.allocate(make_uint3(isize,jsize,NBVAR));
    Fluxes_y.allocate(make_uint3(isize,jsize,NBVAR));    
    
  } else if (params.implementationVersion == 1) {

    Slopes_x.allocate(make_uint3(isize,jsize,NBVAR));
    Slopes_y.allocate(make_uint3(isize,jsize,NBVAR));

  } else if (params.implementationVersion == 2) {

    Qm_x.allocate(make_uint3(isize,jsize,NBVAR));
    Qm_y.allocate(make_uint3(isize,jsize,NBVAR));
    Qp_x.allocate(make_uint3(isize,jsize,NBVAR));
    Qp_y.allocate(make_uint3(isize,jsize,NBVAR));

  }

  // default riemann solver
  // riemann_solver_fn = &HydroRun::riemann_approx;
  // if (!riemannSolverStr.compare("hllc"))
  //   riemann_solver_fn = &HydroRun::riemann_hllc;
  
  /*
   * initialize hydro array at t=0
   */
  if ( params.problemType == PROBLEM_IMPLODE) {

    init_implode(U);
    init_implode(U2);

  } else if (params.problemType == PROBLEM_BLAST) {

    init_blast(U);
    init_blast(U2);

  } else {

    std::cout << "Problem : " << params.problemType
	      << " is not recognized / implemented in initHydroRun."
	      << std::endl;
    std::cout <<  "Use default - implode" << std::endl;
    init_implode(U);
    init_implode(U2);

  }

  // copy U into U2
  memcpy(U2.data(), U.data(),
	 isize * jsize * NBVAR * sizeof(real_t));
} // HydroRun::HydroRun


// =======================================================
// =======================================================
/**
 *
 */
HydroRun::~HydroRun()
{

} // HydroRun::~HydroRun

// =======================================================
// =======================================================
/**
 * Compute time step satisfying CFL condition.
 *
 * \param[in] useU integer, if 0 use data in U else use U2
 *
 * \return dt time step
 */
real_t HydroRun::compute_dt(int useU)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  const real_t dx = params.dx;
  const real_t dy = params.dy;

  real_t dt;
  real_t invDt = ZERO_F;
  int    i,j;
  real_t *Udata;
  
  // which array is the current one ?
  if (useU == 0)
    Udata = U.data();
  else
    Udata = U2.data();

  // for loop over inner region
  for (j=ghostWidth; j<=jsize-ghostWidth; j++) {
    for (i=ghostWidth; i<=isize-ghostWidth; i++) {

      HydroState uLoc; // conservative    variables in current cell
      HydroState qLoc; // primitive    variables in current cell
      real_t c;
      real_t vx, vy;
  
      int index = i + j*isize;

      // get local conservative variable
      uLoc[ID] = Udata[index+ID*ijsize];
      uLoc[IP] = Udata[index+IP*ijsize];
      uLoc[IU] = Udata[index+IU*ijsize];
      uLoc[IV] = Udata[index+IV*ijsize];

      // get primitive variables in current cell
      computePrimitives(uLoc, &c, qLoc);
      vx = c+FABS(qLoc[IU]);
      vy = c+FABS(qLoc[IV]);

      invDt = FMAX(invDt, vx/dx + vy/dy);

    } // end for i
  } // end for j

  dt = params.settings.cfl/invDt;

  return dt;

} // HydroRun::compute_dt

// =======================================================
// =======================================================
// ///////////////////////////////////////////
// Wrapper to the actual computation routine
// ///////////////////////////////////////////
void HydroRun::godunov_unsplit(int nStep, real_t dt)
{
  
  if ( nStep % 2 == 0 ) {
    godunov_unsplit_cpu(U , U2, dt, nStep);
  } else {
    godunov_unsplit_cpu(U2, U , dt, nStep);
  }
  
} // HydroRun::godunov_unsplit

// =======================================================
// =======================================================
// ///////////////////////////////////////////
// Actual CPU computation of Godunov scheme
// ///////////////////////////////////////////
void HydroRun::godunov_unsplit_cpu(DataArray& data_in, 
				   DataArray& data_out, 
				   real_t dt, 
				   int nStep)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  const real_t dx = params.dx;
  const real_t dy = params.dy;

  
  // local variables
  int i, j;
  
  real_t dtdx;
  real_t dtdy;
  
  real_t* Qdata = Q.data();
  real_t *Udata = data_out.data();

  dtdx = dt / dx;
  dtdy = dt / dy;

  // fill ghost cell in data_in
  boundaries_timer.start();
  make_boundaries(data_in);
  boundaries_timer.stop();
    
  // copy data_in into data_out (not necessary)
  // data_out = data_in;
  memcpy(data_out.data(), data_in.data(),
	 isize * jsize * NBVAR * sizeof(real_t));

  // start main computation
  godunov_timer.start();

  // convert conservative variable into primitives ones for the entire domain
  convertToPrimitives(data_in, Q);
  
  if (params.implementationVersion == 0) {
    
    for ( j=ghostWidth; j<=jsize-ghostWidth; j++ ) {
      for ( i=ghostWidth; i<=isize-ghostWidth; i++ ) {
	
	// local primitive variables
	HydroState qLoc; // local primitive variables

	// local primitive variables in neighbor cell
	HydroState qLocNeighbor;

	// local primitive variables in neighborbood
	HydroState qNeighbors_0;
	HydroState qNeighbors_1;
	HydroState qNeighbors_2;
	HydroState qNeighbors_3;

	// Local slopes and neighbor slopes
	HydroState dqX;
	HydroState dqY;
	HydroState dqX_neighbor;
	HydroState dqY_neighbor;

	// Local variables for Riemann problems solving
	HydroState qleft;
	HydroState qright;
	HydroState qgdnv;
	HydroState flux_x;
	HydroState flux_y;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// deal with left interface along X !
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int indexBase = i  + isize * j;
      
	// get primitive variables state vector
	{
	  int index = indexBase + ijsize*ID;
	  qLoc[ID]         = Qdata[index];
	  qNeighbors_0[ID] = Qdata[index+1];
	  qNeighbors_1[ID] = Qdata[index-1];
	  qNeighbors_2[ID] = Qdata[index+isize];
	  qNeighbors_3[ID] = Qdata[index-isize];
	  
	  index = indexBase + ijsize*IP;
	  qLoc[IP]         = Qdata[index];
	  qNeighbors_0[IP] = Qdata[index+1];
	  qNeighbors_1[IP] = Qdata[index-1];
	  qNeighbors_2[IP] = Qdata[index+isize];
	  qNeighbors_3[IP] = Qdata[index-isize];
	  
	  index = indexBase + ijsize*IU;
	  qLoc[IU]         = Qdata[index];
	  qNeighbors_0[IU] = Qdata[index+1];
	  qNeighbors_1[IU] = Qdata[index-1];
	  qNeighbors_2[IU] = Qdata[index+isize];
	  qNeighbors_3[IU] = Qdata[index-isize];
	  
	  index = indexBase + ijsize*IV;
	  qLoc[IV]         = Qdata[index];
	  qNeighbors_0[IV] = Qdata[index+1];
	  qNeighbors_1[IV] = Qdata[index-1];
	  qNeighbors_2[IV] = Qdata[index+isize];
	  qNeighbors_3[IV] = Qdata[index-isize];
	}
	slope_unsplit_hydro_2d(qLoc, 
			       qNeighbors_0, qNeighbors_1, 
			       qNeighbors_2, qNeighbors_3,
			       dqX, dqY);
	
	// slopes at left neighbor along X
	{
	  int index = indexBase + ijsize*ID;
	  qLocNeighbor[ID] = Qdata[index-1];
	  qNeighbors_0[ID] = Qdata[index];
	  qNeighbors_1[ID] = Qdata[index-2];
	  qNeighbors_2[ID] = Qdata[index-1+isize];
	  qNeighbors_3[ID] = Qdata[index-1-isize];

	  index = indexBase + ijsize*IP;
	  qLocNeighbor[IP] = Qdata[index-1];
	  qNeighbors_0[IP] = Qdata[index];
	  qNeighbors_1[IP] = Qdata[index-2];
	  qNeighbors_2[IP] = Qdata[index-1+isize];
	  qNeighbors_3[IP] = Qdata[index-1-isize];

	  index = indexBase + ijsize*IU;
	  qLocNeighbor[IU] = Qdata[index-1];
	  qNeighbors_0[IU] = Qdata[index];
	  qNeighbors_1[IU] = Qdata[index-2];
	  qNeighbors_2[IU] = Qdata[index-1+isize];
	  qNeighbors_3[IU] = Qdata[index-1-isize];

	  index = indexBase + ijsize*IV;
	  qLocNeighbor[IV] = Qdata[index-1];
	  qNeighbors_0[IV] = Qdata[index];
	  qNeighbors_1[IV] = Qdata[index-2];
	  qNeighbors_2[IV] = Qdata[index-1+isize];
	  qNeighbors_3[IV] = Qdata[index-1-isize];
	} 
	slope_unsplit_hydro_2d(qLocNeighbor, 
			       qNeighbors_0, qNeighbors_1, 
			       qNeighbors_2, qNeighbors_3,
			       dqX_neighbor, dqY_neighbor);

	//
	// compute reconstructed states at left interface along X
	//
	
	// left interface : right state
	trace_unsplit_2d_along_dir(qLoc,
				   dqX, dqY,
				   dtdx, dtdy, FACE_XMIN, qright);

	// left interface : left state
	trace_unsplit_2d_along_dir(qLocNeighbor,
				   dqX_neighbor,dqY_neighbor,
				   dtdx, dtdy, FACE_XMAX, qleft);

	// Solve Riemann problem at X-interfaces and compute X-fluxes
	//riemann_2d(qleft,qright,qgdnv,flux_x);
	riemann_hllc(qleft,qright,qgdnv,flux_x);
	
	//
	// store fluxes X
	//
	{

	  real_t *Flux = Fluxes_x.data();
	  Flux[i   + j*isize + ID*ijsize] = flux_x[ID] * dtdx;
	  Flux[i   + j*isize + IP*ijsize] = flux_x[IP] * dtdx;
	  Flux[i   + j*isize + IU*ijsize] = flux_x[IU] * dtdx;
	  Flux[i   + j*isize + IV*ijsize] = flux_x[IV] * dtdx;
	  
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// deal with left interface along Y !
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// slopes at left neighbor along Y
	{
	  int index = indexBase + ijsize*ID;
	  qLocNeighbor[ID] = Qdata[index-isize];
	  qNeighbors_0[ID] = Qdata[index+1-isize];
	  qNeighbors_1[ID] = Qdata[index-1-isize];
	  qNeighbors_2[ID] = Qdata[index];
	  qNeighbors_3[ID] = Qdata[index-2*isize];

	  index = indexBase + ijsize*IP;
	  qLocNeighbor[IP] = Qdata[index-isize];
	  qNeighbors_0[IP] = Qdata[index+1-isize];
	  qNeighbors_1[IP] = Qdata[index-1-isize];
	  qNeighbors_2[IP] = Qdata[index];
	  qNeighbors_3[IP] = Qdata[index-2*isize];

	  index = indexBase + ijsize*IU;
	  qLocNeighbor[IU] = Qdata[index-isize];
	  qNeighbors_0[IU] = Qdata[index+1-isize];
	  qNeighbors_1[IU] = Qdata[index-1-isize];
	  qNeighbors_2[IU] = Qdata[index];
	  qNeighbors_3[IU] = Qdata[index-2*isize];

	  index = indexBase + ijsize*IV;
	  qLocNeighbor[IV] = Qdata[index-isize];
	  qNeighbors_0[IV] = Qdata[index+1-isize];
	  qNeighbors_1[IV] = Qdata[index-1-isize];
	  qNeighbors_2[IV] = Qdata[index];
	  qNeighbors_3[IV] = Qdata[index-2*isize];

	}
	
	slope_unsplit_hydro_2d(qLocNeighbor, 
			       qNeighbors_0, qNeighbors_1, 
			       qNeighbors_2, qNeighbors_3,
			       dqX_neighbor, dqY_neighbor);

	//
	// compute reconstructed states at left interface along Y
	//
	
	// left interface : right state
	trace_unsplit_2d_along_dir(qLoc,
				   dqX, dqY,
				   dtdx, dtdy, FACE_YMIN, qright);

	// left interface : left state
	trace_unsplit_2d_along_dir(qLocNeighbor,
				   dqX_neighbor,dqY_neighbor,
				   dtdx, dtdy, FACE_YMAX, qleft);

	// Solve Riemann problem at Y-interfaces and compute Y-fluxes
	swapValues(&(qleft[IU]) ,&(qleft[IV]) );
	swapValues(&(qright[IU]),&(qright[IV]));
	//riemann_2d(qleft,qright,&qgdnv,&flux_y);
	riemann_hllc(qleft,qright,qgdnv,flux_y);

	//
	// store fluxes Y
	//
	{

	  real_t *Flux = Fluxes_y.data();
	  Flux[i   + j*isize + ID*ijsize] = flux_y[ID] * dtdy;
	  Flux[i   + j*isize + IP*ijsize] = flux_y[IP] * dtdy;
	  Flux[i   + j*isize + IU*ijsize] = flux_y[IU] * dtdy;
	  Flux[i   + j*isize + IV*ijsize] = flux_y[IV] * dtdy;

	}

      } // end for j
    } // end for i

    // actual hydro update
    for ( j=ghostWidth; j<jsize-ghostWidth; j++ ) {
      for ( i=ghostWidth; i<isize-ghostWidth; i++ ) {

	real_t *Flux_x = Fluxes_x.data();
	real_t *Flux_y = Fluxes_y.data();
		
	Udata[i + j*isize + ID*ijsize] +=  Flux_x[i   + j*isize + ID*ijsize];
	Udata[i + j*isize + IP*ijsize] +=  Flux_x[i   + j*isize + IP*ijsize];
	Udata[i + j*isize + IU*ijsize] +=  Flux_x[i   + j*isize + IU*ijsize];
	Udata[i + j*isize + IV*ijsize] +=  Flux_x[i   + j*isize + IV*ijsize];
	
	Udata[i + j*isize + ID*ijsize] -=  Flux_x[i+1 + j*isize + ID*ijsize];
	Udata[i + j*isize + IP*ijsize] -=  Flux_x[i+1 + j*isize + IP*ijsize];
	Udata[i + j*isize + IU*ijsize] -=  Flux_x[i+1 + j*isize + IU*ijsize];
	Udata[i + j*isize + IV*ijsize] -=  Flux_x[i+1 + j*isize + IV*ijsize];

	Udata[i + j*isize + ID*ijsize] +=  Flux_y[i + j*isize + ID*ijsize];
	Udata[i + j*isize + IP*ijsize] +=  Flux_y[i + j*isize + IP*ijsize];
	Udata[i + j*isize + IU*ijsize] +=  Flux_y[i + j*isize + IV*ijsize]; //
	Udata[i + j*isize + IV*ijsize] +=  Flux_y[i + j*isize + IU*ijsize]; //
	
	Udata[i + j*isize + ID*ijsize] -=  Flux_y[i + (j+1)*isize + ID*ijsize];
	Udata[i + j*isize + IP*ijsize] -=  Flux_y[i + (j+1)*isize + IP*ijsize];
	Udata[i + j*isize + IU*ijsize] -=  Flux_y[i + (j+1)*isize + IV*ijsize]; //
	Udata[i + j*isize + IV*ijsize] -=  Flux_y[i + (j+1)*isize + IU*ijsize]; //

      } // end for j
    } // end for i

    
  } else if (params.implementationVersion == 1) {

    real_t *slopes_x = Slopes_x.data();
    real_t *slopes_y = Slopes_y.data();

    // let's store the slopes
    for ( j=ghostWidth-1; j<=jsize-ghostWidth; j++ ) {
      for ( i=ghostWidth-1; i<=isize-ghostWidth; i++ ) {
	
	// local primitive variables
	HydroState qLoc; // local primitive variables

	// local primitive variables in neighborbood
	HydroState qNeighbors_0;
	HydroState qNeighbors_1;
	HydroState qNeighbors_2;
	HydroState qNeighbors_3;

	// Local slopes and neighbor slopes
	HydroState dqX; for (int index=0; index<NBVAR; ++index) dqX[index]=0.0;
	HydroState dqY; for (int index=0; index<NBVAR; ++index) dqY[index]=0.0;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// deal with left interface along X !
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int indexBase = i  + isize * j;
      
	// get primitive variables state vector
	{

	  int index = indexBase + ijsize*ID;
	  qLoc[ID]         = Qdata[index];
	  qNeighbors_0[ID] = Qdata[index+1];
	  qNeighbors_1[ID] = Qdata[index-1];
	  qNeighbors_2[ID] = Qdata[index+isize];
	  qNeighbors_3[ID] = Qdata[index-isize];

	  index = indexBase + ijsize*IP;
	  qLoc[IP]         = Qdata[index];
	  qNeighbors_0[IP] = Qdata[index+1];
	  qNeighbors_1[IP] = Qdata[index-1];
	  qNeighbors_2[IP] = Qdata[index+isize];
	  qNeighbors_3[IP] = Qdata[index-isize];

	  index = indexBase + ijsize*IU;
	  qLoc[IU]         = Qdata[index];
	  qNeighbors_0[IU] = Qdata[index+1];
	  qNeighbors_1[IU] = Qdata[index-1];
	  qNeighbors_2[IU] = Qdata[index+isize];
	  qNeighbors_3[IU] = Qdata[index-isize];

	  index = indexBase + ijsize*IV;
	  qLoc[IV]         = Qdata[index];
	  qNeighbors_0[IV] = Qdata[index+1];
	  qNeighbors_1[IV] = Qdata[index-1];
	  qNeighbors_2[IV] = Qdata[index+isize];
	  qNeighbors_3[IV] = Qdata[index-isize];

	}

	slope_unsplit_hydro_2d(qLoc, 
			       qNeighbors_0, qNeighbors_1, 
			       qNeighbors_2, qNeighbors_3,
			       dqX, dqY);

	// copy back slopes in global arrays
	{

	  int index = indexBase + ijsize*ID;
	  slopes_x[index] = dqX[ID];
	  slopes_y[index] = dqY[ID];

	  index = indexBase + ijsize*IP;
	  slopes_x[index] = dqX[IP];
	  slopes_y[index] = dqY[IP];

	  index = indexBase + ijsize*IU;
	  slopes_x[index] = dqX[IU];
	  slopes_y[index] = dqY[IU];
	  
	  index = indexBase + ijsize*IV;
	  slopes_x[index] = dqX[IV];
	  slopes_y[index] = dqY[IV];

	}

      } // end for i
    } // end for j


    // now trace and update along X axis
    for ( j=ghostWidth; j<=jsize-ghostWidth; j++ ) {
      for ( i=ghostWidth; i<=isize-ghostWidth; i++ ) {

	// local primitive variables
	HydroState qLoc; // local primitive variables

	// local primitive variables in neighbor cell
	HydroState qLocNeighbor;

	// Local slopes and neighbor slopes
	HydroState dqX;
	HydroState dqY;
	HydroState dqX_neighbor;
	HydroState dqY_neighbor;

	// Local variables for Riemann problems solving
	HydroState qleft;
	HydroState qright;
	HydroState qgdnv;
	HydroState flux_x;
	//HydroState flux_y;

	int indexBase = i  + isize * j;	

	//
	// compute reconstructed states at left interface along X
	//
	{

	  int index = indexBase + ijsize*ID;
	  qLoc[ID] = Qdata[index];
	  dqX[ID]  = slopes_x[index];
	  dqY[ID]  = slopes_y[index];

	  index = indexBase + ijsize*IP;
	  qLoc[IP] = Qdata[index];
	  dqX[IP]  = slopes_x[index];
	  dqY[IP]  = slopes_y[index];

	  index = indexBase + ijsize*IU;
	  qLoc[IU] = Qdata[index];
	  dqX[IU]  = slopes_x[index];
	  dqY[IU]  = slopes_y[index];

	  index = indexBase + ijsize*IV;
	  qLoc[IV] = Qdata[index];
	  dqX[IV]  = slopes_x[index];
	  dqY[IV]  = slopes_y[index];

	}
	
	// left interface : right state
	trace_unsplit_2d_along_dir(qLoc,
				   dqX, dqY,
				   dtdx, dtdy, FACE_XMIN, qright);

	{

	  int index = indexBase + ijsize*ID;
	  qLocNeighbor[ID] = Qdata[index-1];
	  dqX_neighbor[ID]  = slopes_x[index-1];
	  dqY_neighbor[ID]  = slopes_y[index-1];

	  index = indexBase + ijsize*IP;
	  qLocNeighbor[IP] = Qdata[index-1];
	  dqX_neighbor[IP]  = slopes_x[index-1];
	  dqY_neighbor[IP]  = slopes_y[index-1];

	  index = indexBase + ijsize*IU;
	  qLocNeighbor[IU] = Qdata[index-1];
	  dqX_neighbor[IU]  = slopes_x[index-1];
	  dqY_neighbor[IU]  = slopes_y[index-1];

	  index = indexBase + ijsize*IV;
	  qLocNeighbor[IV] = Qdata[index-1];
	  dqX_neighbor[IV]  = slopes_x[index-1];
	  dqY_neighbor[IV]  = slopes_y[index-1];

	}

	// left interface : left state
	trace_unsplit_2d_along_dir(qLocNeighbor,
				   dqX_neighbor,dqY_neighbor,
				   dtdx, dtdy, FACE_XMAX, qleft);

	// Solve Riemann problem at X-interfaces and compute X-fluxes
	riemann_hllc(qleft,qright,qgdnv,flux_x);
	
	//
	// update hydro array
	//
	{
	  Udata[i-1 + j*isize + ID*ijsize] += -flux_x[ID]*dtdx;
	  Udata[i-1 + j*isize + IP*ijsize] += -flux_x[IP]*dtdx;
	  Udata[i-1 + j*isize + IU*ijsize] += -flux_x[IU]*dtdx;
	  Udata[i-1 + j*isize + IV*ijsize] += -flux_x[IV]*dtdx;
	  
	  Udata[i   + j*isize + ID*ijsize] +=  flux_x[ID]*dtdx;
	  Udata[i   + j*isize + IP*ijsize] +=  flux_x[IP]*dtdx;
	  Udata[i   + j*isize + IU*ijsize] +=  flux_x[IU]*dtdx;
	  Udata[i   + j*isize + IV*ijsize] +=  flux_x[IV]*dtdx;
	}

      } // end for i
    } // end for j 

    // now trace and update along Y axis
    for ( j=ghostWidth; j<=jsize-ghostWidth; j++ ) {
      for ( i=ghostWidth; i<=isize-ghostWidth; i++ ) {

	// local primitive variables
	HydroState qLoc; // local primitive variables

	// local primitive variables in neighbor cell
	HydroState qLocNeighbor;

	// Local slopes and neighbor slopes
	HydroState dqX;
	HydroState dqY;
	HydroState dqX_neighbor;
	HydroState dqY_neighbor;

	// Local variables for Riemann problems solving
	HydroState qleft;
	HydroState qright;
	HydroState qgdnv;
	//HydroState flux_x;
	HydroState flux_y;

	int indexBase = i  + isize * j;	

	//
	// compute reconstructed states at left interface along Y
	//
	{

	  int index = indexBase + ijsize*ID;
	  qLoc[ID] = Qdata[index];
	  dqX[ID]  = slopes_x[index];
	  dqY[ID]  = slopes_y[index];

	  index = indexBase + ijsize*IP;
	  qLoc[IP] = Qdata[index];
	  dqX[IP]  = slopes_x[index];
	  dqY[IP]  = slopes_y[index];

	  index = indexBase + ijsize*IU;
	  qLoc[IU] = Qdata[index];
	  dqX[IU]  = slopes_x[index];
	  dqY[IU]  = slopes_y[index];

	  index = indexBase + ijsize*IV;
	  qLoc[IV] = Qdata[index];
	  dqX[IV]  = slopes_x[index];
	  dqY[IV]  = slopes_y[index];

	}
	
	// left interface : right state
	trace_unsplit_2d_along_dir(qLoc,
				   dqX, dqY,
				   dtdx, dtdy, FACE_YMIN, qright);

	{

	  int index = indexBase + ijsize*ID;
	  qLocNeighbor[ID] = Qdata[index-isize];
	  dqX_neighbor[ID] = slopes_x[index-isize];
	  dqY_neighbor[ID] = slopes_y[index-isize];

	  index = indexBase + ijsize*IP;
	  qLocNeighbor[IP] = Qdata[index-isize];
	  dqX_neighbor[IP] = slopes_x[index-isize];
	  dqY_neighbor[IP] = slopes_y[index-isize];

	  index = indexBase + ijsize*IU;
	  qLocNeighbor[IU] = Qdata[index-isize];
	  dqX_neighbor[IU] = slopes_x[index-isize];
	  dqY_neighbor[IU] = slopes_y[index-isize];

	  index = indexBase + ijsize*IV;
	  qLocNeighbor[IV] = Qdata[index-isize];
	  dqX_neighbor[IV] = slopes_x[index-isize];
	  dqY_neighbor[IV] = slopes_y[index-isize];
	}

	// left interface : left state
	trace_unsplit_2d_along_dir(qLocNeighbor,
				   dqX_neighbor,dqY_neighbor,
				   dtdx, dtdy, FACE_YMAX, qleft);

	// Solve Riemann problem at Y-interfaces and compute Y-fluxes
	swapValues(&(qleft[IU]) ,&(qleft[IV]) );
	swapValues(&(qright[IU]),&(qright[IV]));
	riemann_hllc(qleft,qright,qgdnv,flux_y);

	//
	// update hydro array
	//
	{
	  Udata[i + (j-1)*isize + ID*ijsize] += -flux_y[ID]*dtdy;
	  Udata[i + (j-1)*isize + IP*ijsize] += -flux_y[IP]*dtdy;
	  Udata[i + (j-1)*isize + IU*ijsize] += -flux_y[IV]*dtdy; // IU/IV swapped
	  Udata[i + (j-1)*isize + IV*ijsize] += -flux_y[IU]*dtdy; // IU/IV swapped
	  
	  Udata[i + j*isize     + ID*ijsize] +=  flux_y[ID]*dtdy;
	  Udata[i + j*isize     + IP*ijsize] +=  flux_y[IP]*dtdy;
	  Udata[i + j*isize     + IU*ijsize] +=  flux_y[IV]*dtdy; // IU/IV swapped
	  Udata[i + j*isize     + IV*ijsize] +=  flux_y[IU]*dtdy; // IU/IV swapped
	}

      } // end for i
    } // end for j 
    
  } else if (params.implementationVersion == 2) {
    
    // trace computation: fill arrays qm_x, qm_y, qp_x, qp_y
    computeTrace(data_in, dt);

    // Compute flux via Riemann solver and update (time integration)
    computeFluxesAndUpdate(data_out, dt);

  } // end params.implementationVersion == 2
  
  godunov_timer.stop();
  
} // HydroRun::godunov_unsplit_cpu

// =======================================================
// =======================================================
// ///////////////////////////////////////////////////////////////////
// Convert conservative variables array U into primitive var array Q
// ///////////////////////////////////////////////////////////////////
void HydroRun::convertToPrimitives(DataArray &u,
				   DataArray &q)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;

  int    i,j;
    
  real_t* Udata = u.data();
  real_t* Qdata = q.data();

  for ( j=0; j<jsize; j++) {
    for ( i=0; i<isize; i++) {
      
      HydroState uLoc; // conservative    variables in current cell
      HydroState qLoc; // primitive    variables in current cell
      real_t c;

      int index = i + j*isize;

      // get local conservative variable
      uLoc[ID] = Udata[index+ID*ijsize];
      uLoc[IP] = Udata[index+IP*ijsize];
      uLoc[IU] = Udata[index+IU*ijsize];
      uLoc[IV] = Udata[index+IV*ijsize];

      computePrimitives(uLoc, &c, qLoc);
      
      // copy q state in q global
      Qdata[i+isize*j+ijsize*ID] = qLoc[ID];
      Qdata[i+isize*j+ijsize*IP] = qLoc[IP];
      Qdata[i+isize*j+ijsize*IU] = qLoc[IU];
      Qdata[i+isize*j+ijsize*IV] = qLoc[IV];

    } // end for i
  } // end for j
  
} // HydroRun::convertToPrimitives

// =======================================================
// =======================================================
// ///////////////////////////////////////////////////////////////////
// Compute trace (only used in implementation version 2), i.e.
// fill global array qm_x, qmy, qp_x, qp_y
// ///////////////////////////////////////////////////////////////////
void HydroRun::computeTrace(DataArray &U, real_t dt)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  const real_t dx = params.dx;
  const real_t dy = params.dy;

  // local variables
  int i,j;
  real_t dtdx;
  real_t dtdy;
  real_t* Qdata= Q.data();
  real_t* qm_x = Qm_x.data();
  real_t* qm_y = Qm_y.data();
  real_t* qp_x = Qp_x.data();
  real_t* qp_y = Qp_y.data();

  dtdx = dt / dx;
  dtdy = dt / dy;

  for ( j=1; j<=jsize-2; j++) {
    for ( i=1; i<=isize-2; i++) {

      HydroState qLoc   ; // local primitive variables
      HydroState qPlusX ;
      HydroState qMinusX;
      HydroState qPlusY ;
      HydroState qMinusY;

      HydroState dqX; for (int index=0; index<NBVAR; ++index) dqX[index]=0.0;
      HydroState dqY; for (int index=0; index<NBVAR; ++index) dqY[index]=0.0;

      HydroState qmX;
      HydroState qmY;
      HydroState qpX;
      HydroState qpY;
      
      int indexBase = i  + isize * j;
      
      // get primitive variables state vector
      {
	int index = indexBase + ijsize*ID;
	qLoc   [ID] = Qdata[index];
	qPlusX [ID] = Qdata[index+1];
	qMinusX[ID] = Qdata[index-1];
	qPlusY [ID] = Qdata[index+isize];
	qMinusY[ID] = Qdata[index-isize];

	index = indexBase + ijsize*IP;
	qLoc   [IP] = Qdata[index];
	qPlusX [IP] = Qdata[index+1];
	qMinusX[IP] = Qdata[index-1];
	qPlusY [IP] = Qdata[index+isize];
	qMinusY[IP] = Qdata[index-isize];

	index = indexBase + ijsize*IU;
	qLoc   [IU] = Qdata[index];
	qPlusX [IU] = Qdata[index+1];
	qMinusX[IU] = Qdata[index-1];
	qPlusY [IU] = Qdata[index+isize];
	qMinusY[IU] = Qdata[index-isize];

	index = indexBase + ijsize*IV;
	qLoc   [IV] = Qdata[index];
	qPlusX [IV] = Qdata[index+1];
	qMinusX[IV] = Qdata[index-1];
	qPlusY [IV] = Qdata[index+isize];
	qMinusY[IV] = Qdata[index-isize];

      } // 
      
      // get hydro slopes dq
      slope_unsplit_hydro_2d(qLoc, 
			     qPlusX, qMinusX, 
			     qPlusY, qMinusY, 
			     dqX, dqY);
      
      // compute qm, qp
      trace_unsplit_hydro_2d(qLoc, 
			     dqX, dqY,
			     dtdx,dtdy, 
			     qmX, qmY,
			     qpX, qpY);

      // store qm, qp : only what is really needed
      {
	int index = indexBase + ijsize*ID;
	qm_x[index] = qmX[ID];
	qp_x[index] = qpX[ID];
	qm_y[index] = qmY[ID];
	qp_y[index] = qpY[ID];

	index = indexBase + ijsize*IP;
	qm_x[index] = qmX[IP];
	qp_x[index] = qpX[IP];
	qm_y[index] = qmY[IP];
	qp_y[index] = qpY[IP];

	index = indexBase + ijsize*IU;
	qm_x[index] = qmX[IU];
	qp_x[index] = qpX[IU];
	qm_y[index] = qmY[IU];
	qp_y[index] = qpY[IU];

	index = indexBase + ijsize*IV;
	qm_x[index] = qmX[IV];
	qp_x[index] = qpX[IV];
	qm_y[index] = qmY[IV];
	qp_y[index] = qpY[IV];
      } // end 

    } // end for i
  } // end for j

} // HydroRun::computeTrace

// =======================================================
// =======================================================
// //////////////////////////////////////////////////////////////////
// Compute flux via Riemann solver and update (time integration)
// //////////////////////////////////////////////////////////////////
void HydroRun::computeFluxesAndUpdate(DataArray &u, 
				      real_t dt)
{
 
  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  const real_t dx = params.dx;
  const real_t dy = params.dy;

  // local variables
  int i,j;
  
  real_t dtdx = dt / dx;
  real_t dtdy = dt / dy;

  real_t* udata = u.data();
  real_t* qm_x = Qm_x.data();
  real_t* qm_y = Qm_y.data();
  real_t* qp_x = Qp_x.data();
  real_t* qp_y = Qp_y.data();

  for ( j=ghostWidth; j<=jsize-ghostWidth; j++ ) {
    for ( i=ghostWidth; i<=isize-ghostWidth; i++ ) {

      HydroState qleft, qright;
      HydroState flux_x, flux_y;
      HydroState qgdnv;
      int index;

      //
      // Solve Riemann problem at X-interfaces and compute
      // X-fluxes
      //
      index = i-1 + isize * j;
      qleft[ID]   = qm_x[index + ijsize * ID];
      qleft[IP]   = qm_x[index + ijsize * IP];
      qleft[IU]   = qm_x[index + ijsize * IU];
      qleft[IV]   = qm_x[index + ijsize * IV];
      
      index = i + isize * j;
      qright[ID]  = qp_x[index + ijsize * ID];
      qright[IP]  = qp_x[index + ijsize * IP];
      qright[IU]  = qp_x[index + ijsize * IU];
      qright[IV]  = qp_x[index + ijsize * IV];
      
      // compute hydro flux_x
      riemann_hllc(qleft,qright,qgdnv,flux_x);

      //
      // Solve Riemann problem at Y-interfaces and compute Y-fluxes
      //
      index = i + isize * (j-1);
      qleft[ID]   = qm_y[index + ijsize * ID];
      qleft[IP]   = qm_y[index + ijsize * IP];
      qleft[IU]   = qm_y[index + ijsize * IV]; // watchout IU, IV permutation
      qleft[IV]   = qm_y[index + ijsize * IU]; // watchout IU, IV permutation

      index = i + isize * j;
      qright[ID]  = qp_y[index + ijsize * ID];
      qright[IP]  = qp_y[index + ijsize * IP];
      qright[IU]  = qp_y[index + ijsize * IV]; // watchout IU, IV permutation
      qright[IV]  = qp_y[index + ijsize * IU]; // watchout IU, IV permutation
      
      // compute hydro flux_y
      riemann_hllc(qleft,qright,qgdnv,flux_y);
            
      //
      // update hydro array
      //
      index = i-1 + isize*j;
      udata[index + ijsize * ID] += - flux_x[ID]*dtdx;
      udata[index + ijsize * IP] += - flux_x[IP]*dtdx;
      udata[index + ijsize * IU] += - flux_x[IU]*dtdx;
      udata[index + ijsize * IV] += - flux_x[IV]*dtdx;

      index = i + isize*j;
      udata[index + ijsize * ID] +=   flux_x[ID]*dtdx;
      udata[index + ijsize * IP] +=   flux_x[IP]*dtdx;
      udata[index + ijsize * IU] +=   flux_x[IU]*dtdx;
      udata[index + ijsize * IV] +=   flux_x[IV]*dtdx;

      index = i + isize*(j-1);
      udata[index + ijsize * ID] += - flux_y[ID]*dtdy;
      udata[index + ijsize * IP] += - flux_y[IP]*dtdy;
      udata[index + ijsize * IU] += - flux_y[IV]*dtdy; // watchout IU and IV swapped
      udata[index + ijsize * IV] += - flux_y[IU]*dtdy; // watchout IU and IV swapped

      index = i + isize*j;
      udata[index + ijsize * ID] +=   flux_y[ID]*dtdy;
      udata[index + ijsize * IP] +=   flux_y[IP]*dtdy;
      udata[index + ijsize * IU] +=   flux_y[IV]*dtdy; // watchout IU and IV swapped
      udata[index + ijsize * IV] +=   flux_y[IU]*dtdy; // watchout IU and IV swapped

    } // end for i
  } // end for j
 
} // computeFluxesAndUpdate

// =======================================================
// =======================================================
/**
 * Hydrodynamical Implosion Test.
 * http://www.astro[IP]rinceton.edu/~jstone/Athena/tests/implode/Implode.html
 */
void HydroRun::init_implode(DataArray &u)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  
  const real_t xmin = params.xmin;
  const real_t ymin = params.ymin;
  const real_t dx = params.dx;
  const real_t dy = params.dy;

  const real_t gamma0 = params.settings.gamma0;
  
  real_t* udata = u.data();

  for (int j=0; j<jsize; j++) {
    real_t y = ymin + dy/2 + (j-ghostWidth)*dy;

    for (int i=0; i<isize; i++) {
      real_t x = xmin + + dx/2 + (i-ghostWidth)*dx;

      int index = i + isize * j;

      real_t tmp = x + y;
      if (tmp > 0.5) {
	udata[index + ID*ijsize] = 1.0;
	udata[index + IP*ijsize] = 1.0/(gamma0-1.0);
	udata[index + IU*ijsize] = 0.0;
	udata[index + IV*ijsize] = 0.0;
      } else {
	udata[index + ID*ijsize] = 0.125;
	udata[index + IP*ijsize] = 0.14/(gamma0-1.0);
	udata[index + IU*ijsize] = 0.0;
	udata[index + IV*ijsize] = 0.0;
      }
    }
  }

} // init_implode

// =======================================================
// =======================================================
/**
 * Hydrodynamical blast Test.
 * http://www.astro.princeton.edu/~jstone/Athena/tests/blast/blast.html
 */
void HydroRun::init_blast(DataArray &u)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  const real_t xmin = params.xmin;
  const real_t ymin = params.ymin;
  const real_t dx = params.dx;
  const real_t dy = params.dy;
    
  // blast problem parameters
  const real_t blast_radius = params.blast_radius;
  const real_t radius2      = blast_radius*blast_radius;
  const real_t blast_center_x    = params.blast_center_x;
  const real_t blast_center_y    = params.blast_center_y;
  const real_t blast_density_in  = params.blast_density_in;
  const real_t blast_density_out = params.blast_density_out;
  const real_t blast_pressure_in = params.blast_pressure_in;
  const real_t blast_pressure_out= params.blast_pressure_out;

  const real_t gamma0 = params.settings.gamma0;

  real_t* udata = u.data();

  for (int j=0; j<jsize; j++) {
    real_t y = ymin + dy/2 + (j-ghostWidth)*dy;

    for (int i=0; i<isize; i++) {
      real_t x = xmin + dx/2 + (i-ghostWidth)*dx;

      int index = i + isize * j;

      real_t d2 = 
	(x-blast_center_x)*(x-blast_center_x)+
	(y-blast_center_y)*(y-blast_center_y);

      if (d2 < radius2) {
	udata[index + ID*ijsize] = blast_density_in;
	udata[index + IP*ijsize] = blast_pressure_in/(gamma0-1.0);
	udata[index + IU*ijsize] = 0.0;
	udata[index + IV*ijsize] = 0.0;
      } else {
	udata[index + ID*ijsize] = blast_density_out;
	udata[index + IP*ijsize] = blast_pressure_out/(gamma0-1.0);
	udata[index + IU*ijsize] = 0.0;
	udata[index + IV*ijsize] = 0.0;
      }
    }
  }

} // HydroRun::init_blast

// =======================================================
// =======================================================
// ///////////////////////////////////////////////////////
// output routine (VTK file format, ASCII, VtkImageData)
// ///////////////////////////////////////////////////////
void HydroRun::saveVTK(DataArray &u,
		       int iStep,
		       std::string name)
{

  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int nx = params.nx;
  const int ny = params.ny;
  const int imin = params.imin;
  const int imax = params.imax;
  const int jmin = params.jmin;
  const int jmax = params.jmax;
  const int ghostWidth = params.ghostWidth;
  
  // local variables
  int i,j,iVar;
  std::string outputDir    = configMap.getString("output", "outputDir", "./");
  std::string outputPrefix = configMap.getString("output", "outputPrefix", "output");
    
  // check scalar data type
  bool useDouble = false;

  real_t* udata = u.data();

  if (sizeof(real_t) == sizeof(double)) {
    useDouble = true;
  }
  
  // write iStep in string stepNum
  std::ostringstream stepNum;
  stepNum.width(7);
  stepNum.fill('0');
  stepNum << iStep;
  
  // concatenate file prefix + file number + suffix
  std::string filename     = outputDir + "/" + outputPrefix+"_"+stepNum.str() + ".vti";
  
  // open file 
  std::fstream outFile;
  outFile.open(filename.c_str(), std::ios_base::out);
  
  // write header
  outFile << "<?xml version=\"1.0\"?>\n";
  if (isBigEndian()) {
    outFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n";
  } else {
    outFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  }

  // write mesh extent
  outFile << "  <ImageData WholeExtent=\""
	  << 0 << " " << nx << " "
	  << 0 << " " << ny << " "
	  << 0 << " " << 0  << "\" "
	  <<  "Origin=\""
	  << params.xmin << " " << params.ymin << " " << 0.0 << "\" "
	  << "Spacing=\""
    	  << params.dx << " " << params.dy << " " << 0.0 << "\">\n";
  outFile << "  <Piece Extent=\""
	  << 0 << " " << nx << " "
	  << 0 << " " << ny << " "
	  << 0 << " " << 0  << " "    
	  << "\">\n";

  outFile << "    <PointData>\n";
  outFile << "    </PointData>\n";
  outFile << "    <CellData>\n";

  // write data array (ascii), remove ghost cells
  for ( iVar=0; iVar<NBVAR; iVar++) {
    outFile << "    <DataArray type=\"";
    if (useDouble)
      outFile << "Float64";
    else
      outFile << "Float32";
    outFile << "\" Name=\"" << varNames[iVar] << "\" format=\"ascii\" >\n";
    
    for ( j=jmin+ghostWidth; j<=jmax-ghostWidth; j++) {
      for ( i=imin+ghostWidth; i<=imax-ghostWidth; i++) {
	outFile << udata[i + isize * j + ijsize * iVar] << " ";
      }
    }
    outFile << "\n    </DataArray>\n";
  } // end for iVar

  outFile << "    </CellData>\n";

  // write footer
  outFile << "  </Piece>\n";
  outFile << "  </ImageData>\n";
  outFile << "</VTKFile>\n";
  
  outFile.close();

} // HydroRun::saveVTK

// =======================================================
// =======================================================
// //////////////////////////////////////////////////
// Fill ghost cells according to border condition :
// absorbant, reflexive or periodic
// //////////////////////////////////////////////////
void HydroRun::make_boundaries(DataArray &u)
{

  const int nx = params.nx;
  const int ny = params.ny;
  
  const int isize = params.isize;
  const int jsize = params.jsize;
  const int ijsize = params.isize * params.jsize;
  const int ghostWidth = params.ghostWidth;
  
  const int imin = params.imin;
  const int imax = params.imax;
  
  const int jmin = params.jmin;
  const int jmax = params.jmax;
  
  // local variables
  int i,j,i0,j0,iVar;
  real_t sign;

  real_t* udata = u.data();

  {
    // boundary xmin
    int boundary_type = params.boundary_type_xmin;
    
    for ( iVar=0; iVar<NBVAR; iVar++ ) {
      for ( i=0; i<ghostWidth; i++ ) {
	sign=ONE_F;
	if ( boundary_type == BC_DIRICHLET ) {
	  i0=2*ghostWidth-1-i;
	  if (iVar==IU) sign=-ONE_F;
	} else if( boundary_type == BC_NEUMANN ) {
	  i0=ghostWidth;
	} else { // periodic
	  i0=nx+i;
	}
	
	for ( j=jmin+ghostWidth; j<=jmax-ghostWidth; j++ ) {
	  udata[i + isize*j + ijsize*iVar] = udata[i0 + isize*j + ijsize*iVar]*sign;
	} // end for j
      } // end for i
    } // end for iVar
  }

  {
    // boundary xmax
    int boundary_type = params.boundary_type_xmax;

    for ( iVar=0; iVar<NBVAR; iVar++ ) {
      for ( i=nx+ghostWidth; i<=nx+2*ghostWidth-1; i++ ) {
	sign=ONE_F;
	if ( boundary_type == BC_DIRICHLET ) {
	  i0=2*nx+2*ghostWidth-1-i;
	  if (iVar==IU) sign=-ONE_F;
	} else if ( boundary_type == BC_NEUMANN ) {
	  i0=nx+ghostWidth-1;
	} else { // periodic
	  i0=i-nx;
	}
	for ( j=jmin+ghostWidth; j<=jmax-ghostWidth; j++) {
	  udata[i + isize*j + ijsize*iVar] = udata[i0 + isize*j + ijsize*iVar]*sign;
	} // end for j
      } // end for i
    } // end for iVar
  }

  {
    // boundary ymin
    int boundary_type = params.boundary_type_ymin;
    
    for ( iVar=0; iVar<NBVAR; iVar++) {
      for ( j=0; j<ghostWidth; j++ ) {
	sign=ONE_F;
	if ( boundary_type == BC_DIRICHLET ) {
	  j0=2*ghostWidth-1-j;
	  if (iVar==IV) sign=-ONE_F;
	} else if ( boundary_type == BC_NEUMANN ) {
	  j0=ghostWidth;
	} else { // periodic
	  j0=ny+j;
	}
	
	for ( i=imin+ghostWidth; i<=imax-ghostWidth; i++ ) {
	  udata[i + isize*j + ijsize*iVar] = udata[i + isize*j0 + ijsize*iVar]*sign;
	} // end for i
      } // end for j
    } // end for iVar
  }

  {
    // boundary ymax
    int boundary_type = params.boundary_type_ymax;

    for ( iVar=0; iVar<NBVAR; iVar++ ) {
      for ( j=ny+ghostWidth; j<=ny+2*ghostWidth-1; j++ ) {
	sign=ONE_F;
	if ( boundary_type == BC_DIRICHLET ) {
	  j0=2*ny+2*ghostWidth-1-j;
	  if (iVar==IV) sign=-ONE_F;
	} else if ( boundary_type == BC_NEUMANN ) {
	  j0=ny+ghostWidth-1;
	} else { // periodic
	  j0=j-ny;
	}
	
	for ( i=imin+ghostWidth; i<=imax-ghostWidth; i++ ) {
	  udata[i + isize*j + ijsize*iVar] = udata[i + isize*j0 + ijsize*iVar]*sign;
	} // end for i
      } //end for j
    } // end for iVar
  }
  
} // HydroRun::make_boundaries


/**
 * Equation of state:
 * compute pressure p and speed of sound c, from density rho and
 * internal energy eint using the "calorically perfect gas" equation
 * of state : \f$ eint=\frac{p}{\rho (\gamma-1)} \f$
 * Recall that \f$ \gamma \f$ is equal to the ratio of specific heats \f$ \left[
 * c_p/c_v \right] \f$.
 * 
 * @param[in]  rho  density
 * @param[in]  eint internal energy
 * @param[out] p    pressure
 * @param[out] c    speed of sound
 */
void HydroRun::eos(real_t rho, real_t eint, real_t* p, real_t* c)
{
  const real_t gamma0 = params.settings.gamma0;
  const real_t smallp = params.settings.smallp;
    
  *p = FMAX((gamma0 - ONE_F) * rho * eint, rho * smallp);
  *c = SQRT(gamma0 * (*p) / rho);
} // HydroRun::eos

/**
 * Convert conservative variables (rho, rho*u, rho*v, e) to 
 * primitive variables (rho,u,v,p)
 * @param[in]  u  conservative variables array
 * @param[out] q  primitive    variables array (allocated in calling routine, size is constant NBVAR)
 * @param[out] c  local speed of sound
 */
void HydroRun::computePrimitives(const HydroState& u, 
				 real_t* c, 
				 HydroState& q)
{
  const real_t gamma0 = params.settings.gamma0;
  const real_t smallr = params.settings.smallr;
  const real_t smallp = params.settings.smallp;
  
  real_t d, p, ux, uy;

  d = fmax(u[ID], smallr);
  ux = u[IU] / d;
  uy = u[IV] / d;
  
  real_t eken = HALF_F * (ux*ux + uy*uy);
  real_t e = u[IP] / d - eken;
  // if (e < 0) {
  //   printf("FATAL ERROR : hydro eint < 0  : e %f eken %f d %f u %f v %f\n",u.p,eken,u.d,u[IU],u[IV]);
  //   exit(0);
  // }
  
  // compute pressure and speed of sound
  //eos(q.d, e, &(q.p), c);
  p = fmax((gamma0 - 1.0) * d * e, d * smallp);
  *c = sqrt(gamma0 * (p) / d);

  q[ID] = d;
  q[IP] = p;
  q[IU] = ux;
  q[IV] = uy;

} // HydroRun::computePrimitive

/**
 * Trace computations for unsplit Godunov scheme.
 *
 * \param[in] q          : Primitive variables state.
 * \param[in] qNeighbors : state in the neighbor cells (2 neighbors
 * per dimension, in the following order x+, x-, y+, y-, z+, z-)
 * \param[in] c          : local sound speed.
 * \param[in] dtdx       : dt over dx
 * \param[out] qm        : qm state (one per dimension)
 * \param[out] qp        : qp state (one per dimension)
 */
void HydroRun::trace_unsplit_2d(const HydroState& q, 
				const HydroState& qNeighbors_0,
				const HydroState& qNeighbors_1,
				const HydroState& qNeighbors_2,
				const HydroState& qNeighbors_3,
				real_t c, 
				real_t dtdx, 
				real_t dtdy, 
				HydroState& qm_x,
				HydroState& qm_y,
				HydroState& qp_x,
				HydroState& qp_y)
{

  const real_t gamma0 = params.settings.gamma0;
  const real_t smallr = params.settings.smallr;
  
  // first compute slopes
  HydroState dqX, dqY;
  dqX[ID] = 0.0;
  dqX[IP] = 0.0;
  dqX[IU] = 0.0;
  dqX[IV] = 0.0;
  dqY[ID] = 0.0;
  dqY[IP] = 0.0;
  dqY[IU] = 0.0;
  dqY[IV] = 0.0;

  slope_unsplit_hydro_2d(q, 
			 qNeighbors_0, qNeighbors_1, 
			 qNeighbors_2, qNeighbors_3,
			 dqX, dqY);
  
  // Cell centered values
  real_t r =  q[ID];
  real_t p =  q[IP];
  real_t u =  q[IU];
  real_t v =  q[IV];
  
  // TVD slopes in all directions
  real_t drx = dqX[ID];
  real_t dpx = dqX[IP];
  real_t dux = dqX[IU];
  real_t dvx = dqX[IV];
  
  real_t dry = dqY[ID];
  real_t dpy = dqY[IP];
  real_t duy = dqY[IU];
  real_t dvy = dqY[IV];
  
  // source terms (with transverse derivatives)
  real_t sr0 = -u*drx-v*dry - (dux+dvy)*r;
  real_t sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma0*p;
  real_t su0 = -u*dux-v*duy - (dpx    )/r;
  real_t sv0 = -u*dvx-v*dvy - (dpy    )/r;
  
  // Right state at left interface
  qp_x[ID] = r - HALF_F*drx + sr0*dtdx*HALF_F;
  qp_x[IP] = p - HALF_F*dpx + sp0*dtdx*HALF_F;
  qp_x[IU] = u - HALF_F*dux + su0*dtdx*HALF_F;
  qp_x[IV] = v - HALF_F*dvx + sv0*dtdx*HALF_F;
  qp_x[ID] = fmax(smallr, qp_x[ID]);
  
  // Left state at right interface
  qm_x[ID] = r + HALF_F*drx + sr0*dtdx*HALF_F;
  qm_x[IP] = p + HALF_F*dpx + sp0*dtdx*HALF_F;
  qm_x[IU] = u + HALF_F*dux + su0*dtdx*HALF_F;
  qm_x[IV] = v + HALF_F*dvx + sv0*dtdx*HALF_F;
  qm_x[ID] = fmax(smallr, qm_x[ID]);
  
  // Top state at bottom interface
  qp_y[ID] = r - HALF_F*dry + sr0*dtdy*HALF_F;
  qp_y[IP] = p - HALF_F*dpy + sp0*dtdy*HALF_F;
  qp_y[IU] = u - HALF_F*duy + su0*dtdy*HALF_F;
  qp_y[IV] = v - HALF_F*dvy + sv0*dtdy*HALF_F;
  qp_y[ID] = fmax(smallr, qp_y[ID]);
  
  // Bottom state at top interface
  qm_y[ID] = r + HALF_F*dry + sr0*dtdy*HALF_F;
  qm_y[IP] = p + HALF_F*dpy + sp0*dtdy*HALF_F;
  qm_y[IU] = u + HALF_F*duy + su0*dtdy*HALF_F;
  qm_y[IV] = v + HALF_F*dvy + sv0*dtdy*HALF_F;
  qm_y[ID] = fmax(smallr, qm_y[ID]);
  
} // HydroRun::trace_unsplit_2d

/**
 * Trace computations for unsplit Godunov scheme.
 *
 * \param[in] q          : Primitive variables state.
 * \param[in] dqX        : slope along X
 * \param[in] dqY        : slope along Y
 * \param[in] c          : local sound speed.
 * \param[in] dtdx       : dt over dx
 * \param[in] dtdy       : dt over dy
 * \param[in] faceId     : which face will be reconstructed
 * \param[out] qface     : q reconstructed state at cell interface
 */
void HydroRun::trace_unsplit_2d_along_dir(const HydroState& q, 
					  const HydroState& dqX,
					  const HydroState& dqY,
					  real_t dtdx, 
					  real_t dtdy, 
					  int    faceId,
					  HydroState& qface)
{
  
  const real_t gamma0 = params.settings.gamma0;
  const real_t smallr = params.settings.smallr;
  
  // Cell centered values
  real_t r =  q[ID];
  real_t p =  q[IP];
  real_t u =  q[IU];
  real_t v =  q[IV];
  
  // TVD slopes in all directions
  real_t drx = dqX[ID];
  real_t dpx = dqX[IP];
  real_t dux = dqX[IU];
  real_t dvx = dqX[IV];
  
  real_t dry = dqY[ID];
  real_t dpy = dqY[IP];
  real_t duy = dqY[IU];
  real_t dvy = dqY[IV];
  
  // source terms (with transverse derivatives)
  real_t sr0 = -u*drx-v*dry - (dux+dvy)*r;
  real_t sp0 = -u*dpx-v*dpy - (dux+dvy)*gamma0*p;
  real_t su0 = -u*dux-v*duy - (dpx    )/r;
  real_t sv0 = -u*dvx-v*dvy - (dpy    )/r;
  
  if (faceId == FACE_XMIN) {
    // Right state at left interface
    qface[ID] = r - HALF_F*drx + sr0*dtdx*HALF_F;
    qface[IP] = p - HALF_F*dpx + sp0*dtdx*HALF_F;
    qface[IU] = u - HALF_F*dux + su0*dtdx*HALF_F;
    qface[IV] = v - HALF_F*dvx + sv0*dtdx*HALF_F;
    qface[ID] = fmax(smallr, qface[ID]);
  }

  if (faceId == FACE_XMAX) {
    // Left state at right interface
    qface[ID] = r + HALF_F*drx + sr0*dtdx*HALF_F;
    qface[IP] = p + HALF_F*dpx + sp0*dtdx*HALF_F;
    qface[IU] = u + HALF_F*dux + su0*dtdx*HALF_F;
    qface[IV] = v + HALF_F*dvx + sv0*dtdx*HALF_F;
    qface[ID] = fmax(smallr, qface[ID]);
  }
  
  if (faceId == FACE_YMIN) {
    // Top state at bottom interface
    qface[ID] = r - HALF_F*dry + sr0*dtdy*HALF_F;
    qface[IP] = p - HALF_F*dpy + sp0*dtdy*HALF_F;
    qface[IU] = u - HALF_F*duy + su0*dtdy*HALF_F;
    qface[IV] = v - HALF_F*dvy + sv0*dtdy*HALF_F;
    qface[ID] = fmax(smallr, qface[ID]);
  }

  if (faceId == FACE_YMAX) {
    // Bottom state at top interface
    qface[ID] = r + HALF_F*dry + sr0*dtdy*HALF_F;
    qface[IP] = p + HALF_F*dpy + sp0*dtdy*HALF_F;
    qface[IU] = u + HALF_F*duy + su0*dtdy*HALF_F;
    qface[IV] = v + HALF_F*dvy + sv0*dtdy*HALF_F;
    qface[ID] = fmax(smallr, qface[ID]);
  }

} // HydroRun::trace_unsplit_2d_along_dir

/**
 * This another implementation of trace computations for 2D data; it
 * is used when unsplitVersion = 1
 *
 * Note that :
 * - hydro slopes computations are done outside this routine
 *
 * \param[in]  q  primitive variable state vector
 * \param[in]  dq primitive variable slopes
 * \param[in]  dtdx dt divided by dx
 * \param[in]  dtdy dt divided by dy
 * \param[out] qm
 * \param[out] qp
 *
 */
void HydroRun::trace_unsplit_hydro_2d(const HydroState& q,
				      const HydroState& dqX,
				      const HydroState& dqY,
				      real_t dtdx,
				      real_t dtdy,
				      HydroState& qm_x,
				      HydroState& qm_y,
				      HydroState& qp_x,
				      HydroState& qp_y)
{
  
  const real_t gamma0 = params.settings.gamma0;
  const real_t smallr = params.settings.smallr;
  const real_t smallp = params.settings.smallp;

  // Cell centered values
  real_t r = q[ID];
  real_t p = q[IP];
  real_t u = q[IU];
  real_t v = q[IV];

  // Cell centered TVD slopes in X direction
  real_t drx = dqX[ID];  drx *= HALF_F;
  real_t dpx = dqX[IP];  dpx *= HALF_F;
  real_t dux = dqX[IU];  dux *= HALF_F;
  real_t dvx = dqX[IV];  dvx *= HALF_F;
  
  // Cell centered TVD slopes in Y direction
  real_t dry = dqY[ID];  dry *= HALF_F;
  real_t dpy = dqY[IP];  dpy *= HALF_F;
  real_t duy = dqY[IU];  duy *= HALF_F;
  real_t dvy = dqY[IV];  dvy *= HALF_F;

  // Source terms (including transverse derivatives)
  real_t sr0, su0, sv0, sp0;

  /*only true for cartesian grid */
  {
    sr0 = (-u*drx-dux*r)       *dtdx + (-v*dry-dvy*r)       *dtdy;
    su0 = (-u*dux-dpx/r)       *dtdx + (-v*duy      )       *dtdy;
    sv0 = (-u*dvx      )       *dtdx + (-v*dvy-dpy/r)       *dtdy;
    sp0 = (-u*dpx-dux*gamma0*p)*dtdx + (-v*dpy-dvy*gamma0*p)*dtdy;    
  } // end cartesian

  // Update in time the  primitive variables
  r = r + sr0;
  u = u + su0;
  v = v + sv0;
  p = p + sp0;

  // Face averaged right state at left interface
  qp_x[ID] = r - drx;
  qp_x[IU] = u - dux;
  qp_x[IV] = v - dvx;
  qp_x[IP] = p - dpx;
  qp_x[ID] = fmax(smallr,  qp_x[ID]);
  qp_x[IP] = fmax(smallp * qp_x[ID], qp_x[IP]);
  
  // Face averaged left state at right interface
  qm_x[ID] = r + drx;
  qm_x[IU] = u + dux;
  qm_x[IV] = v + dvx;
  qm_x[IP] = p + dpx;
  qm_x[ID] = fmax(smallr,  qm_x[ID]);
  qm_x[IP] = fmax(smallp * qm_x[ID], qm_x[IP]);

  // Face averaged top state at bottom interface
  qp_y[ID] = r - dry;
  qp_y[IU] = u - duy;
  qp_y[IV] = v - dvy;
  qp_y[IP] = p - dpy;
  qp_y[ID] = fmax(smallr,  qp_y[ID]);
  qp_y[IP] = fmax(smallp * qp_y[ID], qp_y[IP]);
  
  // Face averaged bottom state at top interface
  qm_y[ID] = r + dry;
  qm_y[IU] = u + duy;
  qm_y[IV] = v + dvy;
  qm_y[IP] = p + dpy;
  qm_y[ID] = fmax(smallr,  qm_y[ID]);
  qm_y[IP] = fmax(smallp * qm_y[ID], qm_y[IP]);
  
} // HydroRun::trace_unsplit_hydro_2d

/**
 * Compute primitive variables slopes (dqX,dqY) for one component from q and its neighbors.
 * This routine is only used in the 2D UNSPLIT integration and slope_type = 0,1 and 2.
 * 
 * Only slope_type 1 and 2 are supported.
 *
 * \param[in]  q       : current primitive variable
 * \param[in]  qPlusX  : value in the next neighbor cell along XDIR
 * \param[in]  qMinusX : value in the previous neighbor cell along XDIR
 * \param[in]  qPlusY  : value in the next neighbor cell along YDIR
 * \param[in]  qMinusY : value in the previous neighbor cell along YDIR
 * \param[out] dqX     : reference to an array returning the X slopes
 * \param[out] dqY     : reference to an array returning the Y slopes
 *
 */
void HydroRun::slope_unsplit_hydro_2d_scalar(real_t q, 
					     real_t qPlusX,
					     real_t qMinusX,
					     real_t qPlusY,
					     real_t qMinusY,
					     real_t *dqX,
					     real_t *dqY)
{
  const real_t slope_type = params.settings.slope_type;
  
  real_t dlft, drgt, dcen, dsgn, slop, dlim;

  // slopes in first coordinate direction
  dlft = slope_type*(q      - qMinusX);
  drgt = slope_type*(qPlusX - q      );
  dcen = HALF_F * (qPlusX - qMinusX);
  dsgn = (dcen >= ZERO_F) ? ONE_F : -ONE_F;
  slop = fmin( FABS(dlft), FABS(drgt) );
  dlim = slop;
  if ( (dlft*drgt) <= ZERO_F )
    dlim = ZERO_F;
  *dqX = dsgn * fmin( dlim, FABS(dcen) );
  
  // slopes in second coordinate direction
  dlft = slope_type*(q      - qMinusY);
  drgt = slope_type*(qPlusY - q      );
  dcen = HALF_F * (qPlusY - qMinusY);
  dsgn = (dcen >= ZERO_F) ? ONE_F : -ONE_F;
  slop = fmin( FABS(dlft), FABS(drgt) );
  dlim = slop;
  if ( (dlft*drgt) <= ZERO_F )
    dlim = ZERO_F;
  *dqY = dsgn * fmin( dlim, FABS(dcen) );

} // HydroRun::slope_unsplit_hydro_2d_scalar

/**
 * Compute primitive variables slope (vector dq) from q and its neighbors.
 * This routine is only used in the 2D UNSPLIT integration and slope_type = 0,1 and 2.
 * 
 * Only slope_type 1 and 2 are supported.
 *
 * \param[in]  q       : current primitive variable state
 * \param[in]  qPlusX  : state in the next neighbor cell along XDIR
 * \param[in]  qMinusX : state in the previous neighbor cell along XDIR
 * \param[in]  qPlusY  : state in the next neighbor cell along YDIR
 * \param[in]  qMinusY : state in the previous neighbor cell along YDIR
 * \param[out] dqX     : reference to an array returning the X slopes
 * \param[out] dqY     : reference to an array returning the Y slopes
 *
 */
void HydroRun::slope_unsplit_hydro_2d(const HydroState& q, 
				      const HydroState& qPlusX, 
				      const HydroState& qMinusX,
				      const HydroState& qPlusY,
				      const HydroState& qMinusY,
				      HydroState& dqX,
				      HydroState& dqY)
{
  
  const real_t slope_type = params.settings.slope_type;
  
  if (slope_type==0) {
    
    dqX[ID] = ZERO_F;
    dqX[IP] = ZERO_F;
    dqX[IU] = ZERO_F;
    dqX[IV] = ZERO_F;

    dqY[ID] = ZERO_F;
    dqY[IP] = ZERO_F;
    dqY[IU] = ZERO_F;
    dqY[IV] = ZERO_F;

    return;
  }

  if (slope_type==1 || slope_type==2) {  // minmod or average

    slope_unsplit_hydro_2d_scalar( q[ID], qPlusX[ID], qMinusX[ID], qPlusY[ID], qMinusY[ID], &(dqX[ID]), &(dqY[ID]));
    slope_unsplit_hydro_2d_scalar( q[IP], qPlusX[IP], qMinusX[IP], qPlusY[IP], qMinusY[IP], &(dqX[IP]), &(dqY[IP]));
    slope_unsplit_hydro_2d_scalar( q[IU], qPlusX[IU], qMinusX[IU], qPlusY[IU], qMinusY[IU], &(dqX[IU]), &(dqY[IU]));
    slope_unsplit_hydro_2d_scalar( q[IV], qPlusX[IV], qMinusX[IV], qPlusY[IV], qMinusY[IV], &(dqX[IV]), &(dqY[IV]));

  } // end slope_type == 1 or 2
  
} // HydroRun::slope_unsplit_hydro_2d

/**
 * Compute cell fluxes from the Godunov state
 * \param[in]  qgdnv input Godunov state
 * \param[out] flux  output flux vector
 */
void HydroRun::cmpflx(const HydroState& qgdnv, 
		      HydroState& flux)
{
  const real_t gamma0 = params.settings.gamma0;
  
  // Compute fluxes
  // Mass density
  flux[ID] = qgdnv[ID] * qgdnv[IU];
  
  // Normal momentum
  flux[IU] = flux[ID] * qgdnv[IU] + qgdnv[IP];
  
  // Transverse momentum
  flux[IV] = flux[ID] * qgdnv[IV];

  // Total energy
  real_t entho = ONE_F / (gamma0 - ONE_F);
  real_t ekin;
  ekin = HALF_F * qgdnv[ID] * (qgdnv[IU]*qgdnv[IU] + qgdnv[IV]*qgdnv[IV]);
  
  real_t etot = qgdnv[IP] * entho + ekin;
  flux[IP] = qgdnv[IU] * (etot + qgdnv[IP]);

} // HydroRun::cmpflx

/** 
 * Riemann solver, equivalent to riemann_approx in RAMSES (see file
 * godunov_utils.f90 in RAMSES).
 * 
 * @param[in] qleft  : input left state
 * @param[in] qright : input right state
 * @param[out] qgdnv : output Godunov state
 * @param[out] flux  : output flux
 */
void HydroRun::riemann_approx(const HydroState& qleft, 
			      const HydroState& qright,
			      HydroState& qgdnv, 
			      HydroState& flux)		    
{

  real_t gamma0  = params.settings.gamma0;
  real_t gamma6  = params.settings.gamma6;
  real_t smallr  = params.settings.smallr;
  real_t smallc  = params.settings.smallc;
  real_t smallp  = params.settings.smallp;
  real_t smallpp = params.settings.smallpp;

  // Pressure, density and velocity
  real_t rl = fmax(qleft [ID], smallr);
  real_t ul =      qleft [IU];
  real_t pl = fmax(qleft [IP], rl*smallp);
  real_t rr = fmax(qright[ID], smallr);
  real_t ur =      qright[IU];
  real_t pr = fmax(qright[IP], rr*smallp);
  
  // Lagrangian sound speed
  real_t cl = gamma0*pl*rl;
  real_t cr = gamma0*pr*rr;
  
  // First guess
  real_t wl = SQRT(cl);
  real_t wr = SQRT(cr);
  real_t pstar = fmax(((wr*pl+wl*pr)+wl*wr*(ul-ur))/(wl+wr), (real_t) ZERO_F);
  real_t pold = pstar;
  real_t conv = ONE_F;
  
  // Newton-Raphson iterations to find pstar at the required accuracy
  for(int iter = 0; (iter < 10 /*niter_riemann*/) && (conv > 1e-6); ++iter)
    {
      real_t wwl = SQRT(cl*(ONE_F+gamma6*(pold-pl)/pl));
      real_t wwr = SQRT(cr*(ONE_F+gamma6*(pold-pr)/pr));
      real_t ql = 2.0f*wwl*wwl*wwl/(wwl*wwl+cl);
      real_t qr = 2.0f*wwr*wwr*wwr/(wwr*wwr+cr);
      real_t usl = ul-(pold-pl)/wwl;
      real_t usr = ur+(pold-pr)/wwr;
      real_t delp = fmax(qr*ql/(qr+ql)*(usl-usr),-pold);
      
      pold = pold+delp;
      conv = FABS(delp/(pold+smallpp));	 // Convergence indicator
    }
  
  // Star region pressure
  // for a two-shock Riemann problem
  pstar = pold;
  wl = SQRT(cl*(ONE_F+gamma6*(pstar-pl)/pl));
  wr = SQRT(cr*(ONE_F+gamma6*(pstar-pr)/pr));
  
  // Star region velocity
  // for a two shock Riemann problem
  real_t ustar = HALF_F * (ul + (pl-pstar)/wl + ur - (pr-pstar)/wr);
  
  // Left going or right going contact wave
  real_t sgnm = COPYSIGN(ONE_F, ustar);
  
  // Left or right unperturbed state
  real_t ro, uo, po, wo;
  if(sgnm > ZERO_F)
    {
      ro = rl;
      uo = ul;
      po = pl;
      wo = wl;
    }
  else
    {
      ro = rr;
      uo = ur;
      po = pr;
      wo = wr;
    }
  real_t co = fmax(smallc, SQRT(FABS(gamma0*po/ro)));
  
  // Star region density (Shock, fmax prevents vacuum formation in star region)
  real_t rstar = fmax((real_t) (ro/(ONE_F+ro*(po-pstar)/(wo*wo))), (real_t) (smallr));
  // Star region sound speed
  real_t cstar = fmax(smallc, SQRT(FABS(gamma0*pstar/rstar)));
  
  // Compute rarefaction head and tail speed
  real_t spout  = co    - sgnm*uo;
  real_t spin   = cstar - sgnm*ustar;
  // Compute shock speed
  real_t ushock = wo/ro - sgnm*uo;
  
  if(pstar >= po)
    {
      spin  = ushock;
      spout = ushock;
    }
  
  // Sample the solution at x/t=0
  real_t scr = fmax(spout-spin, smallc+FABS(spout+spin));
  real_t frac = HALF_F * (ONE_F + (spout + spin)/scr);

  if (frac != frac) /* Not a Number */
    frac = 0.0;
  else
    frac = frac >= 1.0 ? 1.0 : frac <= 0.0 ? 0.0 : frac;
  
  qgdnv[ID] = frac*rstar + (ONE_F-frac)*ro;
  qgdnv[IU] = frac*ustar + (ONE_F-frac)*uo;
  qgdnv[IP] = frac*pstar + (ONE_F-frac)*po;
  
  if(spout < ZERO_F)
    {
      qgdnv[ID] = ro;
      qgdnv[IU] = uo;
      qgdnv[IP] = po;
    }
  
  if(spin > ZERO_F)
    {
      qgdnv[ID] = rstar;
      qgdnv[IU] = ustar;
      qgdnv[IP] = pstar;
    }
  
  // transverse velocity
  if(sgnm > ZERO_F)
    {
      qgdnv[IV] = qleft[IV];
    }
  else
    {
      qgdnv[IV] = qright[IV];
    }
  
  cmpflx(qgdnv, flux);
  
} // HydroRun::riemann_approx

/** 
 * Riemann solver HLLC
 *
 * @param[in] qleft : input left state
 * @param[in] qright : input right state
 * @param[out] qgdnv : output Godunov state
 * @param[out] flux  : output flux
 */
void HydroRun::riemann_hllc(const HydroState& qleft, 
			    const HydroState& qright,
			    HydroState& qgdnv, 
			    HydroState& flux)
{

  real_t gamma0  = params.settings.gamma0;
  real_t gamma6  = params.settings.gamma6;
  real_t smallr  = params.settings.smallr;
  real_t smallc  = params.settings.smallc;
  real_t smallp  = params.settings.smallp;
  real_t smallpp = params.settings.smallpp;

  const real_t entho = ONE_F / (gamma0 - ONE_F);
  
  // Left variables
  real_t rl = fmax(qleft[ID], smallr);
  real_t pl = fmax(qleft[IP], rl*smallp);
  real_t ul =      qleft[IU];
    
  real_t ecinl = HALF_F*rl*ul*ul;
  ecinl += HALF_F*rl*qleft[IV]*qleft[IV];

  real_t etotl = pl*entho+ecinl;
  real_t ptotl = pl;

  // Right variables
  real_t rr = fmax(qright[ID], smallr);
  real_t pr = fmax(qright[IP], rr*smallp);
  real_t ur =      qright[IU];

  real_t ecinr = HALF_F*rr*ur*ur;
  ecinr += HALF_F*rr*qright[IV]*qright[IV];
  
  real_t etotr = pr*entho+ecinr;
  real_t ptotr = pr;
    
  // Find the largest eigenvalues in the normal direction to the interface
  real_t cfastl = SQRT(fmax(gamma0*pl/rl,smallc*smallc));
  real_t cfastr = SQRT(fmax(gamma0*pr/rr,smallc*smallc));

  // Compute HLL wave speed
  real_t SL = fmin(ul,ur) - fmax(cfastl,cfastr);
  real_t SR = fmax(ul,ur) + fmax(cfastl,cfastr);

  // Compute lagrangian sound speed
  real_t rcl = rl*(ul-SL);
  real_t rcr = rr*(SR-ur);
    
  // Compute acoustic star state
  real_t ustar    = (rcr*ur   +rcl*ul   +  (ptotl-ptotr))/(rcr+rcl);
  real_t ptotstar = (rcr*ptotl+rcl*ptotr+rcl*rcr*(ul-ur))/(rcr+rcl);

  // Left star region variables
  real_t rstarl    = rl*(SL-ul)/(SL-ustar);
  real_t etotstarl = ((SL-ul)*etotl-ptotl*ul+ptotstar*ustar)/(SL-ustar);
    
  // Right star region variables
  real_t rstarr    = rr*(SR-ur)/(SR-ustar);
  real_t etotstarr = ((SR-ur)*etotr-ptotr*ur+ptotstar*ustar)/(SR-ustar);
    
  // Sample the solution at x/t=0
  real_t ro, uo, ptoto, etoto;
  if (SL > ZERO_F) {
    ro=rl;
    uo=ul;
    ptoto=ptotl;
    etoto=etotl;
  } else if (ustar > ZERO_F) {
    ro=rstarl;
    uo=ustar;
    ptoto=ptotstar;
    etoto=etotstarl;
  } else if (SR > ZERO_F) {
    ro=rstarr;
    uo=ustar;
    ptoto=ptotstar;
    etoto=etotstarr;
  } else {
    ro=rr;
    uo=ur;
    ptoto=ptotr;
    etoto=etotr;
  }
      
  // Compute the Godunov flux
  flux[ID] = ro*uo;
  flux[IU] = ro*uo*uo+ptoto;
  flux[IP] = (etoto+ptoto)*uo;
  if (flux[ID] > ZERO_F) {
    flux[IV] = flux[ID]*qleft[IV];
  } else {
    flux[IV] = flux[ID]*qright[IV];
  }
  
} // HydroRun::riemann_hllc
