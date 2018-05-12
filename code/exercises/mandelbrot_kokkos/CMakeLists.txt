cmake_minimum_required(VERSION 3.1)

# project(mandelbrot_kokkos_cmake CXX)

# # C++11 is for Kokkos
# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_EXTENSIONS OFF)

# # first buid kokkos
# add_subdirectory(kokkos)

# # pass Kokkos include directories to our target application
# include_directories(${Kokkos_INCLUDE_DIRS_RET})

add_subdirectory(serial)
#add_subdirectory(kokkos_basic)
#add_subdirectory(kokkos_mdrange)
#add_subdirectory(kokkos_pipelined0)
#add_subdirectory(kokkos_pipelined1)
