/**
 *
 */
#ifndef HYDRO_RUN_H_
#define HYDRO_RUN_H_

#include <array>

#include "real_type.h" // for type real_t
#include "Arrays.h"
#include "Timer.h"
#include "HydroParams.h"
#include "hydro_common.h" // for NBVAR


/**
 * Main hydrodynamics data structure.
 */
class HydroRun
{

public:

  HydroRun(HydroParams& params, ConfigMap& configMap);
  virtual ~HydroRun();
  
  // hydroParams
  HydroParams& params;
  ConfigMap&   configMap;
  

  DataArray U;    /*!< hydrodynamics conservative variables arrays */
  DataArray U2;   /*!< hydrodynamics conservative variables arrays */
  DataArray Q;    /*!< hydrodynamics primitive    variables array  */

  /* implementation 0 */
  DataArray Fluxes_x;
  DataArray Fluxes_y;
  
  /* implementation 1 only */
  DataArray Slopes_x; /*!< implementation 1 only */
  DataArray Slopes_y; /*!< implementation 1 only */

  /* implementation 2 only */
  DataArray Qm_x; /*!< hydrodynamics Riemann states array implementation 2 */
  DataArray Qm_y; /*!< hydrodynamics Riemann states array */
  DataArray Qp_x; /*!< hydrodynamics Riemann states array */
  DataArray Qp_y; /*!< hydrodynamics Riemann states array */

  //riemann_solver_t riemann_solver_fn; /*!< riemann solver function pointer */

  Timer boundaries_timer, godunov_timer;
  
  // methods
  real_t compute_dt(int useU);
  
  void godunov_unsplit(int nStep, real_t dt);
  
  void godunov_unsplit_cpu(DataArray& data_in, 
			   DataArray& data_out, 
			   real_t dt, 
			   int nStep);
  
  void convertToPrimitives(DataArray &u, DataArray &q);
  
  void computeTrace(DataArray &u, real_t dt);
  
  void computeFluxesAndUpdate(DataArray &u, 
			      real_t dt);
  
  void init_implode(DataArray &u);
  void init_blast(DataArray &u);
  
  void saveVTK(DataArray &u, int iStep, std::string name);
  
  void make_boundaries(DataArray &u);

  // utility routines
  void eos(real_t rho, real_t eint, real_t* p, real_t* c);
  
  void computePrimitives(const HydroState& u, real_t* c, HydroState& q);
  
  void trace_unsplit_2d(const HydroState& q, 
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
			HydroState& qp_y);
  
  void trace_unsplit_2d_along_dir(const HydroState& q, 
				  const HydroState& dqX,
				  const HydroState& dqY,
				  real_t dtdx, 
				  real_t dtdy, 
				  int    faceId,
				  HydroState& qface);
  
  void trace_unsplit_hydro_2d(const HydroState& q,
			      const HydroState& dqX,
			      const HydroState& dqY,
			      real_t dtdx,
			      real_t dtdy,
			      HydroState& qm_x,
			      HydroState& qm_y,
			      HydroState& qp_x,
			      HydroState& qp_y);
  
  void slope_unsplit_hydro_2d_scalar(real_t q, 
				     real_t qPlusX,
				     real_t qMinusX,
				     real_t qPlusY,
				     real_t qMinusY,
				     real_t *dqX,
				     real_t *dqY);
  
  void slope_unsplit_hydro_2d(const HydroState& q, 
			      const HydroState& qPlusX, 
			      const HydroState& qMinusX,
			      const HydroState& qPlusY,
			      const HydroState& qMinusY,
			      HydroState& dqX,
			      HydroState& dqY);
  
  void cmpflx(const HydroState& qgdnv, 
	      HydroState& flux);
  
  void riemann_approx(const HydroState& qleft,
		      const HydroState& qright,
		      HydroState& qgdnv, 
		      HydroState& flux);
  
  void riemann_hllc(const HydroState& qleft,
		    const HydroState& qright,
		    HydroState& qgdnv,
		    HydroState& flux);
  
}; // class HydroRun

#endif // HYDRO_RUN_H_
