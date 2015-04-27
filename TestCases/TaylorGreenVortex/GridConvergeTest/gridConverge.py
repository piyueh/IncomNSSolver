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

ErrReU = numpy.zeros(0, dtype=numpy.float)
ErrReV = numpy.zeros(0, dtype=numpy.float)
ErrReP = numpy.zeros(0, dtype=numpy.float)

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
    p = p[:, :, 0].T


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

    p += p_ext(Xp[0, 0], Yp[0, 0], t)

    ErrU = numpy.append(ErrU, numpy.sqrt(numpy.sum((u-u_e)[1:-1, 1:-1]**2)/(Nx*Ny)))
    ErrV = numpy.append(ErrV, numpy.sqrt(numpy.sum((v-v_e)[1:-1, 1:-1]**2)/(Nx*Ny)))
    ErrP = numpy.append(ErrP, numpy.sqrt(numpy.sum((p-p_e)[1:-1, 1:-1]**2)/(Nx*Ny)))

    ErrReU = numpy.append(ErrReU, numpy.max(u[1:-1, 1:-1]))
    ErrReV = numpy.append(ErrReV, numpy.max(v[1:-1, 1:-1]))
    ErrReP = numpy.append(ErrReP, numpy.max(p[1:-1, 1:-1]))


for i, n in enumerate(N):
    print(n, ErrReU[i], ErrReV[i], ErrReP[i])

for i, n in enumerate(N):
    ErrReU[i] = numpy.abs(ErrReU[i] - ErrReU[-1])
    ErrReV[i] = numpy.abs(ErrReV[i] - ErrReV[-1])
    ErrReP[i] = numpy.abs(ErrReP[i] - ErrReP[-1])


for i, n in enumerate(N):
    print(n, ErrU[i], ErrV[i], ErrP[i])

for i, n in enumerate(N):
    print(n, ErrReU[i], ErrReV[i], ErrReP[i])


DL = numpy.array([Lx / float(n) for n in N])


pyplot.figure()
pyplot.title(r"$L_2$ norm of absolute error to the exact solution")
pyplot.loglog(DL, ErrU, "kx", markersize=10, label="Err(u)")
pyplot.loglog(DL, ErrV, "k^", markersize=10, label="Err(v)")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrU[0]/(2**n) for n in range(5)]), 
        "r--", label="1st order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrU[0]/(4**n) for n in range(5)]), 
        "r--", label="2nd order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrU[0]/(8**n) for n in range(5)]), 
        "r--", label="3rd order")
#pyplot.axis("equal")
pyplot.legend(loc=0)
pyplot.grid(True)


pyplot.figure()
pyplot.title(r"$L_2$ norm of absolute error to the finest grid")
pyplot.loglog(DL[:-1], ErrReU[:-1], "kx", markersize=10, label="ErrRe(u)")
pyplot.loglog(DL[:-1], ErrReV[:-1], "k^", markersize=10, label="ErrRe(v)")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrReU[0]/(2**n) for n in range(5)]), 
        "r--", label="1st order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrReU[0]/(4**n) for n in range(5)]), 
        "r--", label="2nd order")
pyplot.loglog(numpy.array([DL[0]/(2**n) for n in range(5)]), 
        numpy.array([ErrReU[0]/(8**n) for n in range(5)]), 
        "r--", label="3rd order")
#pyplot.axis("equal")
pyplot.legend(loc=0)
pyplot.grid(True)

pyplot.show()
