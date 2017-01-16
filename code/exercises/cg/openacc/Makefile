CC=pgc++
CFLAGS_ACC=-fast -acc -ta=tesla:managed -Minfo=accel
CFLAGS_OMP=-fast -mp -ta=host -Minfo

all: 
	${CC} ${CFLAGS_ACC} main.cpp -o cg.gpu
	${CC} ${CFLAGS_OMP} main.cpp -o cg.cpu
clean:
	rm -Rf cg.cpu cg.gpu pgprof* gprof.out *.o core
