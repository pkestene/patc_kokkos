# KOKKOS setup (OpenMP or CUDA backend)
# here we assume KOKKOS_PATH is setup by "module load kokkos"
KOKKOS_PATH = $(HOME)/Kokkos/kokkos
KOKKOS_DEVICES ?= "OpenMP"
EXE_NAME=cg

SRC = \
	OpenMPTimer.cpp \
	main.cpp

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = $(EXE_NAME).cuda
KOKKOS_ARCH = "Turing75"
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else ifneq (,$(findstring OpenMP,$(KOKKOS_DEVICES)))
# Use either GNU or XLC compiler (for Power8)
CXX = g++
EXE = $(EXE_NAME).omp
KOKKOS_ARCH = "BDW"
else
# Use either GNU or XLC compiler (for Power8)
CXX = g++
EXE = $(EXE_NAME).serial
KOKKOS_ARCH = "BDW"
endif

CXXFLAGS = -O3
LINK = ${CXX}
LINKFLAGS = -lgomp

DEPFLAGS = -M

OBJ = $(SRC:.cpp=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE) : $(OBJ)  $(KOKKOS_LINK_DEPENDS)
	echo $@
	$(LINK) $(KOKKOS_LDFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE) $(LINKFLAGS)

clean: kokkos-clean
	rm -f *.o

cleanexe :
	rm -f *.cuda *.serial *.omp

cleandata:
	rm -f *.ppm

cleanall: clean cleandata cleanexe

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<

test: $(EXE)
	./$(EXE)
