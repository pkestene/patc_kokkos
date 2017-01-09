EXE_PREFIX=mandelbrot

# Use either GNU or XLC compiler (for Power8)
CXX = g++
CXXFLAGS = -g -O3 -I. -Wall --std=c++11
LDFLAGS=-lgomp
#CXX = xlc++
#CXXFLAGS = -O5 -I. -qarch=pwr8 -qsmp=omp -std=c++11
EXE = $(EXE_PREFIX).serial

SRCDIR = $(shell pwd)
SRC = \
	OpenMPTimer.cpp \
	main.cpp

TMP_OBJ = $(SRC:.c=.o)
OBJ     = $(TMP_OBJ:.cpp=.o)

.DEFAULT_GOAL=$(EXE)

$(EXE) : $(OBJ)
	echo $@
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJ) $(EXE)

cleandata:
	rm -f *.ppm

cleanall: clean
	rm -f *.ppm

%.o:    %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
