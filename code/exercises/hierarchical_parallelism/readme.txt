The exercise is to parallelize a dot product using a league of teams of threads.

1. Have a look at how this is done in regular OpenMP, e.g. :
https://github.com/OpenMP/Examples/blob/v4.5.0/sources/Example_teams.1.c

2. Use the kokkos project template, and example 01_thread_teams_lambda (from kokkos source), to adapt this simple algorithm to kokkos.


additionnal documentation:
https://github.com/kokkos/kokkos/wiki/HierarchicalParallelism#822-basic-kernels
