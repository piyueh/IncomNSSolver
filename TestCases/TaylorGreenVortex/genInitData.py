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


def TaylorGreen(t, Nx, Ny):

    Lx = Ly = 2 * numpy.pi
    dx = Lx / Nx
    dy = Ly / Ny

    xp = numpy.linspace(dx/2, Lx-dx/2, Nx)
    yp = numpy.linspace(dy/2, Ly-dy/2, Ny)
    Xp, Yp = numpy.meshgrid(xp, yp)

    xu = numpy.linspace(0, Lx, Nx+1)
    yu = numpy.linspace(dy/2, Ly-dy/2, Ny)
    Xu, Yu = numpy.meshgrid(xu, yu)

    xv = numpy.linspace(dx/2, Lx-dx/2, Nx)
    yv = numpy.linspace(0, Ly, Ny+1)
    Xv, Yv = numpy.meshgrid(xv, yv)

    u_e = u_ext(Xu, Yu, t)
    v_e = v_ext(Xv, Yv, t)
    p_e = p_ext(Xp, Yp, t)

    return Xp, Yp, p_e, Xu, Yu, u_e, Xv, Yv, v_e


def outputInitData(
