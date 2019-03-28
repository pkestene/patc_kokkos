cmake_minimum_required(VERSION 3.3)

project(kokkos-handson-4 CXX)

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

# build the miniapp
option (BUILD_SOLUTION "Enable / disable building solution" OFF)

if (BUILD_SOLUTION)
   add_subdirectory(src_solution)
else()
   add_subdirectory(src)
endif()
