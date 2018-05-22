#ifndef HYDRO_COMMON_H_
#define HYDRO_COMMON_H_

#include "real_type.h"
#include "Arrays.h"

#define NDIM 2
#define NBVAR 4

/// a POD data structure to store local conservative / primitive variables
using HydroState = std::array<real_t, NBVAR>;

/// the DataArray alias
using DataArray = HostArray<real_t>;

#endif // HYDRO_COMMON_H_
