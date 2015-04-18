#======================================================================
# test
#======================================================================

CC = clang
CXX = clang++

BIN = test.out
OBJS = Misc.o Boundary.o Mesh.o PoissonSolver.o NSSolverEuler.o io.o main.o

.PHONY: clean debug release

debug: CFLAGS = -std=c++11 -g
debug:
	make ${BIN} CFLAGS="${CFLAGS}"

release: CFLAGS = -std=c++11 -O3
release:
	make ${BIN} CFLAGS="${CFLAGS}"

${BIN}: ${OBJS}
	${CXX} ${CFLAGS} -o ${BIN} ${OBJS}

%.o: %.cpp
	${CXX} ${CFLAGS} -o $@ -c $<

clean:
	-@rm ${BIN} ${OBJS}
