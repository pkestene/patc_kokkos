cmake_minimum_required(VERSION 3.1)

project(mandelbrot_kokkos_cmake CXX)

# C++11 is for Kokkos
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

#
# sources
#

# first buid kokkos
add_subdirectory(kokkos)

# pass Kokkos include directories to our target application
include_directories(${Kokkos_INCLUDE_DIRS_RET})

# include header for main application
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# mandelbrot_kokkos executable
#
if (KOKKOS_ENABLE_CUDA)
  set(EXE_NAME "mandelbrot_kokkos_pipeline0.cuda")
else()
  set(EXE_NAME "mandelbrot_kokkos_pipeline0.openmp")
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
