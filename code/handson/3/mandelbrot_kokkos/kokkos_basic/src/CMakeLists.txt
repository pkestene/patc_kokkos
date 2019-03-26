#
# mandelbrot_kokkos executable
#
if (KOKKOS_ENABLE_CUDA)
  set(EXE_NAME "mandelbrot_kokkos_basic.cuda")
else()
  set(EXE_NAME "mandelbrot_kokkos_basic.openmp")
endif()

set (EXE_NAME_SRC
  kokkos_shared.h
  constants.h
  mandelbrot.h
  main.cpp)

if (KOKKOS_ENABLE_CUDA)
  set (EXE_NAME_SRC
    ${EXE_NAME_SRC}
    CudaTimer.h)
else()
  set (EXE_NAME_SRC
    ${EXE_NAME_SRC}
    OpenMPTimer.cpp
    OpenMPTimer.h)
endif()

add_executable(${EXE_NAME}
  ${EXE_NAME_SRC})

target_link_libraries(${EXE_NAME} kokkos)
if (KOKKOS_ENABLE_HWLOC)
  target_link_libraries(${EXE_NAME} hwloc)
endif()
