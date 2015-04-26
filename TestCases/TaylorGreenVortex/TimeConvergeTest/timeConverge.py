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

Case = [0.00025, 0.000125, 0.0000625, 0.00003125]

File = {0.00025: "DT25E-5/2000.txt", 0.000125: "DT12.5E-5/4000.txt", 
        0.0000625: "DT6.25E-5/8000.txt", 0.00003125: "DT3.125E-5/16000.txt"}

L2norm = {i: 0 for i in Case}

for dt in Case:

    f = open(File[dt], "r")

    t = float(f.readline())
    print(t)

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

    xp = numpy.linspace(dx/2, Lx-dx/2, Nx)
    yp = numpy.linspace(dy/2, Ly-dy/2, Ny)
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


    L2norm[dt] = numpy.sqrt(numpy.sum((u - u_e)**2)) / (Nx * Ny)


print(L2norm)
