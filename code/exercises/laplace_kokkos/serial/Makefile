CPPFLAGS = -DUSE_DOUBLE

CC=g++

all: run

laplace2d: main.cpp params.h
	$(CC) $(CPPFLAGS) main.cpp -o $@

run: laplace2d
	./laplace2d
