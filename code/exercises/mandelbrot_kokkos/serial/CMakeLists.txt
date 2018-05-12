project(mandelbrot_serial CXX)

# C++11 is for Kokkos
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(mandelbrot.serial
  main.cpp
  constants.h
  SimpleTimer.h
  SimpleTimer.cpp)
