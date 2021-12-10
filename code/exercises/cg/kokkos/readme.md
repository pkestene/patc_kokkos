# How to build ?

## serial

```shell
make KOKKOS_DEVICE=Serial
```

## OpenMP

```shell
make KOKKOS_DEVICE=OpenMP
```

## Cuda

```shell
# be sure to nvcc in your PATH
# e.g.
# module load cuda/11.5
make KOKKOS_DEVICE=Cuda
```
