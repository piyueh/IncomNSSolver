#======================================================================
# test
#======================================================================

CC = clang
CXX = clang++

SRC = src
BPATH = bin
OPATH = objs

IPATH = -I${MKLROOT}/include
LPATH = -L${MKLROOT}/lib/intel64
LIB = -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -lcholmod

BIN = IncomNSSolver
OBJS = Misc.o Boundary.o Mesh.o Data.o Solid.o PoissonSolver.o NSSolver.o io.o main.o

.PHONY: clean debug release

debug: CFLAGS = -std=c++11 -g -DCYLINDER
debug:
	@if [ ! -e ${OPATH} ]; then mkdir ${OPATH}; fi
	@if [ ! -e ${BPATH} ]; then mkdir ${BPATH}; fi
	make ${BPATH}/${BIN} CFLAGS="${CFLAGS}"

release: CFLAGS = -std=c++11 -O3 -fopenmp -march=native -mtune=native -m64 \
	-DNDEBUG -DEIGEN_NO_DEBUG
release:
	echo ${MKLROOT}
	@if [ ! -e ${OPATH} ]; then mkdir ${OPATH}; fi
	@if [ ! -e ${BPATH} ]; then mkdir ${BPATH}; fi
	make ${BPATH}/${BIN} CFLAGS="${CFLAGS}"

${BPATH}/${BIN}: $(foreach i, ${OBJS}, ${OPATH}/${i})
	${CXX} ${CFLAGS} ${IPATH} ${LPATH} ${LIB} \
		-o ${BPATH}/${BIN} $(foreach i, ${OBJS}, ${OPATH}/${i})

${OPATH}/%.o: ${SRC}/%.cpp
	${CXX} ${CFLAGS} ${IPATH} -o $@ -c $<

clean:
	-@rm -r ${OPATH}
	-@rm -r ${BPATH}
