Euler2d is a fluid dynamics mini application solving the Euler system
(compressible fluid dynamics without viscosity) on a cartesian grid.

You have 3 directories:
- euler2d_serial provides a serial implementation. See Readme.md
  to get build and run instruction

- euler2d_kokkos provides an incomplete kokkos parallel port of euler2d_serial
  Your job is fill the missing part (just grep "TODO" in source file)

- euler2d_kokkos_solution provide a fully functionnal solution
