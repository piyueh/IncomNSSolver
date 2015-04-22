#======================================================================
# test
#======================================================================

CC = clang
CXX = clang++

SRC = src
BPATH = bin
OPATH = objs

BIN = IncomNSSolver
OBJS = Misc.o Boundary.o Mesh.o Data.o PoissonSolver.o NSSolver.o io.o main.o

.PHONY: clean debug release

debug: CFLAGS = -std=c++11 -g
debug:
	@if [ ! -e ${OPATH} ]; then mkdir ${OPATH}; fi
	@if [ ! -e ${BPATH} ]; then mkdir ${BPATH}; fi
	make ${BPATH}/${BIN} CFLAGS="${CFLAGS}"

release: CFLAGS = -std=c++11 -O3 -march=native -DNDEBUG -DEIGEN_NO_DEBUG
release:
	@if [ ! -e ${OPATH} ]; then mkdir ${OPATH}; fi
	@if [ ! -e ${BPATH} ]; then mkdir ${BPATH}; fi
	make ${BPATH}/${BIN} CFLAGS="${CFLAGS}"

${BPATH}/${BIN}: $(foreach i, ${OBJS}, ${OPATH}/${i})
	${CXX} ${CFLAGS} -o ${BPATH}/${BIN} $(foreach i, ${OBJS}, ${OPATH}/${i})

${OPATH}/%.o: ${SRC}/%.cpp
	${CXX} ${CFLAGS} -o $@ -c $<

clean:
	-@rm -r ${OPATH}
	-@rm -r ${BPATH}
