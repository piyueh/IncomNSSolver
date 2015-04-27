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

NBase = [101, 201, 401]

DT_Base = {
        101: ["0.0064", "0.0032", "0.0016", "0.0008", "0.0004", "0.0002", 
              "0.0001", "0.00005", "0.000025", "0.0000125", "0.00000625"],

        201: ["0.0064", "0.0032", "0.0016", "0.0008", "0.0004", "0.0002", 
              "0.0001", "0.00005", "0.000025", "0.0000125", "0.00000625"],

        401: ["0.0008", "0.0004", "0.0002", "0.0001", "0.00005", "0.000025", 
              "0.0000125", "0.00000625", "0.000003125", "0.0000015625", "0.00000078125"],
           }

for N in NBase:

    DT = DT_Base[N]

    File = {dt: "N"+str(N)+"/DT"+dt+"/Data.txt" for dt in DT}

    uAll = []
    vAll = []
    pAll = []

    ErrU = numpy.zeros(0, dtype=numpy.float)
    ErrV = numpy.zeros(0, dtype=numpy.float)
    ErrP = numpy.zeros(0, dtype=numpy.float)

    ErrReU = numpy.zeros(0, dtype=numpy.float)
    ErrReV = numpy.zeros(0, dtype=numpy.float)
    ErrReP = numpy.zeros(0, dtype=numpy.float)

    for dt in DT:

        f = open(File[dt], "r")

        t = float(f.readline())
        print(N, dt, t)

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

        uAll.append(u)
        vAll.append(v)
        pAll.append(p)

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


    for i, dt in enumerate(DT):

        ErrReU = numpy.append(ErrReU, numpy.abs(uAll[i] - uAll[-1])[Ny/4+1, Nx/2+1])
        ErrReV = numpy.append(ErrReV, numpy.abs(vAll[i] - vAll[-1])[Ny/2+1, Nx/4+1])
        ErrReP = numpy.append(ErrReP, numpy.abs(pAll[i] - pAll[-1])[Ny/2+1, Nx/2+1])

    for i, dt in enumerate(DT):
        print(dt, ErrU[i], ErrV[i], ErrP[i])

    for i, dt in enumerate(DT):
        print(dt, ErrReU[i], ErrReV[i], ErrReP[i])


    DT = numpy.array([float(dt) for dt in DT])


    pyplot.figure()
    pyplot.title(r"$L_\infty$ of absolute error " + 
                  "\n relative to the exact solution, " +
                 r"$N_x=N_y=$" + str(N),
                 fontsize=16)
    pyplot.xlabel(r"$\Delta t$", fontsize=18)
    pyplot.ylabel(r"$L_\infty$", fontsize=18)
    pyplot.loglog(DT, ErrU, "kx", markersize=10, label="u-velocity")
    pyplot.loglog(DT, ErrV, "k^", markersize=10, label="v-velocity")
    pyplot.loglog(DT, ErrP, "kd", markersize=10, label="pressure")
    pyplot.loglog(DT, 
            numpy.array([ErrU[0]/(2**n) for n in range(DT.size)]), 
            "r--", label="1st order")
    pyplot.loglog(DT, 
            numpy.array([ErrU[0]/(4**n) for n in range(DT.size)]), 
            "g--", label="2nd order")
    pyplot.loglog(DT, 
            numpy.array([ErrU[0]/(8**n) for n in range(DT.size)]), 
            "b--", label="3rd order")
    pyplot.loglog(DT, 
            numpy.array([ErrP[0]/(2**n) for n in range(DT.size)]), 
            "r--")
    pyplot.loglog(DT, 
            numpy.array([ErrP[0]/(4**n) for n in range(DT.size)]), 
            "g--")
    pyplot.loglog(DT, 
            numpy.array([ErrP[0]/(8**n) for n in range(DT.size)]), 
            "b--")
    pyplot.xlim(DT.min()/10, DT.max()*10)
    pyplot.ylim(ErrU.min()/10, 10**(numpy.log10(ErrP.max())+2))
    pyplot.legend(loc=9, ncol=3, mode="expand", numpoints=1)
    pyplot.grid(True)
    pyplot.savefig("N"+str(N)+"_TimeErrExact.png", format="png")


    pyplot.figure()
    pyplot.title(r"$L_\infty$ of absolute error " + 
                  "\n relative to the finest time step, " +
                 r"$N_x=N_y=$" + str(N),
                 fontsize=16)
    pyplot.xlabel(r"$\Delta t$", fontsize=18)
    pyplot.ylabel(r"$L_\infty$", fontsize=18)
    pyplot.loglog(DT[:-1], ErrReU[:-1], "kx", markersize=10, label="u-velocity")
    pyplot.loglog(DT[:-1], ErrReV[:-1], "k^", markersize=10, label="v-velocity")
    pyplot.loglog(DT[:-1], ErrReP[:-1], "kd", markersize=10, label="pressure")
    pyplot.loglog(DT, 
            numpy.array([ErrReU[0]/(2**n) for n in range(DT.size)]), 
            "r--", label="1st order")
    pyplot.loglog(DT, 
            numpy.array([ErrReU[0]/(4**n) for n in range(DT.size)]), 
            "g--", label="2nd order")
    pyplot.loglog(DT, 
            numpy.array([ErrReU[0]/(8**n) for n in range(DT.size)]), 
            "b--", label="3rd order")
    pyplot.loglog(DT, 
            numpy.array([ErrReP[0]/(2**n) for n in range(DT.size)]), 
            "r--")
    pyplot.loglog(DT, 
            numpy.array([ErrReP[0]/(4**n) for n in range(DT.size)]), 
            "g--")
    pyplot.loglog(DT, 
            numpy.array([ErrReP[0]/(8**n) for n in range(DT.size)]), 
            "b--")
    pyplot.xlim(DT.min()/10, DT.max()*10)
    pyplot.ylim(ymax=10**(numpy.log10(ErrReP.max())+6))
    pyplot.legend(loc=9, ncol=3, mode="expand", numpoints=1)
    pyplot.grid(True)
    pyplot.savefig("N"+str(N)+"_TimeErrFinestTime.png", format="png")

pyplot.show()
