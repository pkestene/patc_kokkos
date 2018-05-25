# KOKKOS setup (OpenMP or CUDA backend)
# here we assume KOKKOS_PATH is setup by an environment variable
KOKKOS_PATH ?= .

default: build
	echo "Start Build"

EXE_PREFIX=laplace2d_kokkos_mpi

# assuming OpenMPI
MPI_CXXFLAGS = $(shell mpicxx -showme:compile)
MPI_LDFLAGS  = $(shell mpicxx -showme:link)

NVCC_WRAPPER = $(KOKKOS_PATH)/bin/nvcc_wrapper

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = $(NVCC_WRAPPER)
CXXFLAGS = -ccbin g++ -I. -Wall -O3 $(MPI_CXXFLAGS)
LDFLAGS = -lgomp $(MPI_LDFLAGS)
EXE = $(EXE_PREFIX).cuda
KOKKOS_DEVICES = "Cuda,OpenMP"
KOKKOS_CUDA_OPTIONS += "enable_lambda"
else
# Use either GNU or XLC compiler (for Power8)
CXX = mpicxx
CXXFLAGS = -O3 -I. -Wall --std=c++11
#CXX = xlc++
#CXXFLAGS = -O5 -I. -qarch=pwr8 -qsmp=omp -std=c++11
EXE = $(EXE_PREFIX).omp
KOKKOS_DEVICES = "OpenMP"
endif

SRCDIR = $(shell pwd)
SRC = \
	OpenMPTimer.cpp \
	main.cpp

HEADERS = common.h kokkos_shared.h \
	OpenMPTimer.h CudaTimer.h \
	DataContext.h DataContextKokkos.h \
	laplace2d_kokkos_kernel.h \
	laplace2d_serial_kernel.h \
	params.h check_results.h

TMP_OBJ = $(SRC:.c=.o)
OBJ     = $(TMP_OBJ:.cpp=.o)

$(OBJ): $(HEADERS)

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	echo $@
	$(CXX) $(KOKKOS_LDFLAGS) $(LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $@

clean: kokkos-clean
	rm -f $(OBJ) $(EXE)

cleandata:
	rm -f *.ppm

cleanall: clean
	rm -f *.ppm

# Compilation rules

%.o:    %.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $@
