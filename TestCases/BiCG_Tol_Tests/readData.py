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

N = [101, 201]
DT = [0.0001, 0.0002]
TOL = ["1E-15", "1E-13", "1E-11", "1E-9", "1E-7", "1E-5", "1E-3"] 

File = {n: {dt: {tol: "N"+str(n)+"DT"+str(dt)+"/Tol"+tol+"/Data.txt" for tol in TOL} for dt in DT} for n in N}

ErrU = {n: {dt: numpy.zeros(0, dtype=numpy.float) for dt in DT} for n in N}
ErrV = {n: {dt: numpy.zeros(0, dtype=numpy.float) for dt in DT} for n in N}
ErrP = {n: {dt: numpy.zeros(0, dtype=numpy.float) for dt in DT} for n in N}

for i, n in enumerate(N):
    for j, dt in enumerate(DT): 
        for k, tol in enumerate(TOL):

            f = open(File[n][dt][tol], "r")

            t = float(f.readline())
            print(n, dt, tol, t)

            uN = numpy.array([int(x) for x in f.readline().split()])
            u = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(uN))

            vN = numpy.array([int(x) for x in f.readline().split()])
            v = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(vN))

            wN = numpy.array([int(x) for x in f.readline().split()])
            w = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(wN))

            pN = numpy.array([int(x) for x in f.readline().split()])
            p = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(pN))

            f.close()

            u = u[:, :, 1].T
            v = v[:, :, 1].T
            p = p[:, :, 1].T

            p -= numpy.average(p[1:-1, 1:-1])

            Nx = vN[0] - 2
            Ny = uN[1] - 2

            assert(Nx == n and Ny == n)

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

            ErrU[n][dt] = numpy.append(ErrU[n][dt], numpy.abs(u-u_e)[1:-1, 1:-1].max())
            ErrV[n][dt] = numpy.append(ErrV[n][dt], numpy.abs(v-v_e)[1:-1, 1:-1].max())
            ErrP[n][dt] = numpy.append(ErrP[n][dt], numpy.abs(p-p_e)[1:-1, 1:-1].max())



for i, n in enumerate(N):
    for j, dt in enumerate(DT): 
        for k, tol in enumerate(TOL):
            print(n, dt, tol, ErrU[n][dt][k], ErrV[n][dt][k], ErrP[n][dt][k])


TOL = numpy.array([float(tol) for tol in TOL])

pyplot.figure()
pyplot.title(r"$L_\infty$ of absolute error " + 
        "\n relative to the exact solution: u velocity",
             fontsize=16)
pyplot.xlabel(r"Tolerance", fontsize=18)
pyplot.ylabel(r"$L_\infty$", fontsize=18)

pyplot.loglog(TOL, ErrU[101][0.0001], "b^-", markersize=10, label="N=101, dt=0.0001")
pyplot.loglog(TOL, ErrU[101][0.0002], "r.-", markersize=10, label="N=101, dt=0.0002")
pyplot.loglog(TOL, ErrU[201][0.0001], "b^-.", markersize=10, label="N=201, dt=0.0001")
pyplot.loglog(TOL, ErrU[201][0.0002], "r.-.", markersize=10, label="N=201, dt=0.0002")

pyplot.xlim(1e-16, 1e-2)
pyplot.ylim(1e-6, 1e-2)
pyplot.legend(loc=9, ncol=2, mode="expand", numpoints=1)
pyplot.grid(True)
pyplot.savefig("PoissonTolTest_U.png", format="png")


pyplot.figure()
pyplot.title(r"$L_\infty$ of absolute error " + 
        "\n relative to the exact solution: v velocity",
             fontsize=16)
pyplot.xlabel(r"Tolerance", fontsize=18)
pyplot.ylabel(r"$L_\infty$", fontsize=18)

pyplot.loglog(TOL, ErrV[101][0.0001], "b^-", markersize=10, label="N=101, dt=0.0001")
pyplot.loglog(TOL, ErrV[101][0.0002], "r.-", markersize=10, label="N=101, dt=0.0002")
pyplot.loglog(TOL, ErrV[201][0.0001], "b^-.", markersize=10, label="N=201, dt=0.0001")
pyplot.loglog(TOL, ErrV[201][0.0002], "r.-.", markersize=10, label="N=201, dt=0.0002")

pyplot.xlim(1e-16, 1e-2)
pyplot.ylim(1e-6, 1e-2)
pyplot.legend(loc=9, ncol=2, mode="expand", numpoints=1)
pyplot.grid(True)
pyplot.savefig("PoissonTolTest_V.png", format="png")

pyplot.figure()
pyplot.title(r"$L_\infty$ of absolute error " + 
        "\n relative to the exact solution: Pressure",
             fontsize=16)
pyplot.xlabel(r"Tolerance", fontsize=18)
pyplot.ylabel(r"$L_\infty$", fontsize=18)

pyplot.loglog(TOL, ErrP[101][0.0001], "b^-", markersize=10, label="N=101, dt=0.0001")
pyplot.loglog(TOL, ErrP[101][0.0002], "r.-", markersize=10, label="N=101, dt=0.0002")
pyplot.loglog(TOL, ErrP[201][0.0001], "b^-.", markersize=10, label="N=201, dt=0.0001")
pyplot.loglog(TOL, ErrP[201][0.0002], "r.-.", markersize=10, label="N=201, dt=0.0002")

pyplot.xlim(1e-16, 1e-2)
pyplot.ylim(1e-6, 1e-2)
pyplot.legend(loc=9, ncol=2, mode="expand", numpoints=1)
pyplot.grid(True)
pyplot.savefig("PoissonTolTest_P.png", format="png")

pyplot.show()
