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


def TaylorGreen(t, Nx, Ny, Nz=1):

    Lx = Ly = 2 * numpy.pi
    dx = Lx / Nx
    dy = Ly / Ny

    xp = numpy.linspace(dx/2, Lx-dx/2, Nx)
    yp = numpy.linspace(dy/2, Ly-dy/2, Ny)
    Xp, Yp = numpy.meshgrid(xp, yp)
    Xp = Xp.T
    Yp = Yp.T

    xu = numpy.linspace(-dx, Lx+dx, Nx+3)
    yu = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
    Xu, Yu = numpy.meshgrid(xu, yu)
    Xu = Xu.T
    Yu = Yu.T

    xv = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
    yv = numpy.linspace(-dy, Ly+dy, Ny+3)
    Xv, Yv = numpy.meshgrid(xv, yv)
    Xv = Xv.T
    Yv = Yv.T

    u_e = numpy.zeros((Nx+3, Ny+2, Nz+2))
    v_e = numpy.zeros((Nx+2, Ny+3, Nz+2))
    w_e = numpy.zeros((Nx+2, Ny+2, Nz+3))
    p_e = numpy.zeros((Nx, Ny, Nz))

    u_e[:, :, 1] = u_ext(Xu, Yu, t)
    v_e[:, :, 1] = v_ext(Xv, Yv, t)
    p_e[:, :, 0] = p_ext(Xp, Yp, t)

    return p_e, u_e, v_e, w_e


def outputInitData(t, u, v, w, p, filename):
    '''
    '''
    f = open(filename, 'wb')

    f.write(bytes("TIME " + str(t), "UTF-8"))

    f.write(bytes("\nNu ", "UTF-8"))
    numpy.savetxt(f, u.shape, fmt='%i', delimiter=' ', newline=' ')
    f.write(bytes("\nu ", "UTF-8"))
    numpy.savetxt(f, u.flatten(), delimiter=' ', newline=' ')

    f.write(bytes("\nNv ", "UTF-8"))
    numpy.savetxt(f, v.shape, fmt='%i', delimiter=' ', newline=' ')
    f.write(bytes("\nv ", "UTF-8"))
    numpy.savetxt(f, v.flatten(), delimiter=' ', newline=' ')

    f.write(bytes("\nNw ", "UTF-8"))
    numpy.savetxt(f, w.shape, fmt='%i', delimiter=' ', newline=' ')
    f.write(bytes("\nw ", "UTF-8"))
    numpy.savetxt(f, w.flatten(), delimiter=' ', newline=' ')

    f.write(bytes("\nNp ", "UTF-8"))
    numpy.savetxt(f, p.shape, fmt='%i', delimiter=' ', newline=' ')

    f.close()


t = 0
p, u, v, w = TaylorGreen(t, 160, 160)
outputInitData(t, u, v, w, p, "initData.dat")
