#
# define a macro to build an executable
#
macro(build_exec name)

  # define executable name
  if (KOKKOS_ENABLE_CUDA)
    set(EXE_NAME "${name}.cuda")
  else()
    set(EXE_NAME "${name}.openmp")
  endif()
  
  # define list of source files
  set (EXE_NAME_SRC
    ${name}.cpp
    )

  # add timer
  if (KOKKOS_ENABLE_CUDA)
    set (EXE_NAME_SRC
      ${EXE_NAME_SRC}
      CudaTimer.h
      )
  else()
    set (EXE_NAME_SRC
      ${EXE_NAME_SRC}
      OpenMPTimer.h
      OpenMPTimer.cpp
      )  
  endif()
  
  # define how to build and link
  add_executable(${EXE_NAME}
    ${EXE_NAME_SRC})
  
  target_link_libraries(${EXE_NAME} kokkos)
  if (KOKKOS_ENABLE_HWLOC)
    target_link_libraries(${EXE_NAME} hwloc)
  endif()
  
endmacro(build_exec)

####################################################
build_exec(stencil)


  

