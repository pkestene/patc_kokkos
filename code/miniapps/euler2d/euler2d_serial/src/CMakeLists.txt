add_executable(euler2d
  ../config/inih/ini.h
  ../config/inih/ini.cpp
  ../config/inih/INIReader.h
  ../config/inih/INIReader.cpp
  ../config/ConfigMap.h
  ../config/ConfigMap.cpp
  real_type.h
  common_types.h
  HydroParams.h
  HydroParams.cpp
  HydroRun.h
  HydroRun.cpp
  hydro_common.h
  Timer.h
  Timer.cpp
  main.cpp)

target_include_directories(euler2d
  PUBLIC
  ${CMAKE_SOURCE_DIR})
target_link_libraries(euler2d)

configure_file(test_blast.ini test_blast.ini COPYONLY)
configure_file(test_implode.ini test_implode.ini COPYONLY)
configure_file(test_implode_big.ini test_implode_big.ini COPYONLY)
