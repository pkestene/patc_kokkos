CPPFLAGS = -DUSE_DOUBLE -fopenmp

CC=g++

all: run

laplace2d: main.cpp params.h common.h DataContext.h laplace2d_serial_kernel.h OpenMPTimer.h OpenMPTimer.cpp
	$(CC) $(CPPFLAGS) main.cpp OpenMPTimer.cpp -o $@

run: laplace2d
	./laplace2d
