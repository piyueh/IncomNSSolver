#======================================================================
# test
#======================================================================

CC = clang
CXX = clang++


BIN = test.out
OBJS = Boundary.o PoissonSolver.o main.o

.PHONY: clean debug release

debug: CFLAGS = -std=c++11 -g
debug:
	make clean
	make ${BIN} CFLAGS="${CFLAGS}"

release: CFLAGS = -std=c++11 -O3 -DNDEBUG -DEIGEN_NO_DEBUG
release:
	make clean
	make ${BIN} CFLAGS="${CFLAGS}"

${BIN}: ${OBJS}
	${CXX} ${CFLAGS} -o ${BIN} ${OBJS}

%.o: %.cpp
	${CXX} ${CFLAGS} -o $@ -c $<

clean:
	-@rm ${BIN} ${OBJS}
