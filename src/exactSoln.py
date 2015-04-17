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


Lx = 2 * numpy.pi
Ly = 2 * numpy.pi

Nx = 50
Ny = 50

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


t = 2

p = p_ext(Xp, Yp, t)
u = u_ext(Xu, Yu, t)
v = v_ext(Xv, Yv, t)

uc = (u[:, 1:] + u[:, :-1]) / 2
vc = (v[1:, :] + v[:-1, :]) / 2

print(uc.shape)
print(vc.shape)

pyplot.figure()
pyplot.contour(Xu, Yu, u)

pyplot.figure()
pyplot.contour(Xv, Yv, v)

pyplot.figure()
pyplot.streamplot(Xp, Yp, uc, vc)


pyplot.show()
