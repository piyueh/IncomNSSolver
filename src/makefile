#======================================================================
# test
#======================================================================

CC = clang
CXX = clang++

CFLAGS = -std=c++11 -g

BIN = test.out
OBJS = Boundary.o PoissonSolver.o main.o

.PHONY: clean

${BIN}: ${OBJS}
	${CXX} ${CFLAGS} -o ${BIN} ${OBJS}

%.o: %.cpp
	${CXX} ${CFLAGS} -o $@ -c $<

clean:
	-@rm ${BIN} ${OBJS}
