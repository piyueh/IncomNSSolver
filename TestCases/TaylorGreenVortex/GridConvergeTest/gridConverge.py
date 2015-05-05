import numpy
from matplotlib import pyplot


def u_ext(x, y, t):
    '''
    '''
    return - numpy.exp(- 2. * t) * numpy.cos(x) * numpy.sin(y)


def v_ext(x, y, t):
    '''
    '''
    return numpy.exp(- 2. * t) * numpy.sin(x) * numpy.cos(y)


def p_ext(x, y, t):
    '''
    '''
    return - numpy.exp(-4 * t) * \
            (numpy.cos(2 * x) + numpy.cos(2 * y)) * 0.25



Lx = Ly = 2 * numpy.pi

N = ["25", "51", "101", "151", "201", "251", "301", "351", "401"]

File = {n: "N"+n+"/Data.txt" for n in N}

ErrU = numpy.zeros(0, dtype=numpy.float)
ErrV = numpy.zeros(0, dtype=numpy.float)
ErrP = numpy.zeros(0, dtype=numpy.float)

for i, n in enumerate(N):

    f = open(File[n], "r")

    t = float(f.readline())
    print(n, t)

    uN = numpy.array([int(x) for x in f.readline().split()])
    u = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(uN))

    vN = numpy.array([int(x) for x in f.readline().split()])
    v = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(vN))

    wN = numpy.array([int(x) for x in f.readline().split()])
    w = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(wN))

    pN = numpy.array([int(x) for x in f.readline().split()])
    p = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(pN))

    f.close()

    usave = u.copy()
    vsave = v.copy()

    u = u[:, :, 1].T
    v = v[:, :, 1].T
    p = p[:, :, 1].T

    p -= numpy.average(p[1:-1, 1:-1])

    Nx = vN[0] - 2
    Ny = uN[1] - 2
    dx = Lx / Nx
    dy = Ly / Ny

    xp = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
    yp = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
    Xp, Yp = numpy.meshgrid(xp, yp)

    xu = numpy.linspace(-dx, Lx+dx, Nx+3)
    yu = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
    Xu, Yu = numpy.meshgrid(xu, yu)

    xv = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
    yv = numpy.linspace(-dy, Ly+dy, Ny+3)
    Xv, Yv = numpy.meshgrid(xv, yv)

    u_e = u_ext(Xu, Yu, t)
    v_e = v_ext(Xv, Yv, t)
    p_e = p_ext(Xp, Yp, t)

    ErrU = numpy.append(ErrU, numpy.abs(u-u_e)[1:-1, 1:-1].max())
    ErrV = numpy.append(ErrV, numpy.abs(v-v_e)[1:-1, 1:-1].max())
    ErrP = numpy.append(ErrP, numpy.abs(p-p_e)[1:-1, 1:-1].max())


for i, n in enumerate(N):
    print(n, ErrU[i], ErrV[i], ErrP[i])

DL = numpy.array([Lx / float(n) for n in N])


pyplot.figure()
pyplot.title(r"$L_\infty$ of absolute error " + 
              "\n relative to the exact solution",
             fontsize=16)
pyplot.xlabel(r"$\Delta x$(=$\Delta y$)", fontsize=18)
pyplot.ylabel(r"$L_\infty$", fontsize=18)
pyplot.loglog(DL, ErrU, "kx", markersize=10, label="u-velocity")
pyplot.loglog(DL, ErrV, "k^", markersize=10, label="v-velocity")
pyplot.loglog(DL, ErrP, "kd", markersize=10, label="Pressure")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrU[0]/(2**n) for n in range(5)]), 
        "r--", label="1st order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrU[0]/(4**n) for n in range(5)]), 
        "g--", label="2nd order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrU[0]/(8**n) for n in range(5)]), 
        "b--", label="3rd order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrP[0]/(2**n) for n in range(5)]), 
        "r--")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrP[0]/(4**n) for n in range(5)]), 
        "g--")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrP[0]/(8**n) for n in range(5)]), 
        "b--")
pyplot.xlim(DL.min()/10, DL.max()*10)
pyplot.ylim(ErrU.min()/10, 10**(numpy.log10(ErrP.max())+2))
pyplot.legend(loc=9, ncol=3, mode="expand", numpoints=1)
pyplot.grid(True)
pyplot.savefig("GridErrExact.png", format="png")


pyplot.show()
