# set minimal version
cmake_minimum_required(VERSION 3.18)

project(12_FiniteDifference2DWave
  LANGUAGES C CXX)

# make sure the user doesn't play dirty with symlinks
get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

# disallow in-source builds
if("${srcdir}" STREQUAL "${bindir}")
  message("######################################################")
  message("# ${PROJECT_NAME} should not be configured and built in the source directory")
  message("# You must run cmake from a build directory.")
  message("# For example:")
  message("# mkdir _build ; cd _build")
  message("# run cmake from the build directory.")
  message("######################################################")
  message(FATAL_ERROR "Quitting configuration")
endif()

# set default build type to "Release"
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#
# build exe
#
add_executable(finite_difference_2dwave
  finite_difference_2dwave.cpp)

configure_file(plot_imshow_csv.py plot_imshow_csv.py COPYONLY)

file(GLOB CSV_LIST "${CMAKE_CURRENT_BINARY_DIR}/*.csv")
set_property(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  APPEND
  PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${CSV_LIST})

##################### PRINT CONFIGURE STATUS ######################
message("//===================================================")
message("// ${PROJECT_NAME} build configuration:")
message("//===================================================")
message("")
message("  CMake version          : ${CMAKE_VERSION}")
if (NOT CMAKE_BUILD_TYPE)
  message("  CMake build type       : NOT SET !")
else()
  message("  CMake build type       : ${CMAKE_BUILD_TYPE}")
endif()
message("  CMake install prefix   : ${CMAKE_INSTALL_PREFIX}")
message("  CMake system processor : ${CMAKE_SYSTEM_PROCESSOR}")
message("  CMake system name (OS) : ${CMAKE_SYSTEM_NAME}")
message("")
message("  C++ Compiler : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
if(DEFINED CMAKE_CXX_COMPILER_WRAPPER)
  message("                 ${CMAKE_CXX_COMPILER_WRAPPER}")
endif()
message("                 ${CMAKE_CXX_COMPILER}")
message("  CMAKE_CXX_STANDARD : ${CMAKE_CXX_STANDARD}")
message("")
