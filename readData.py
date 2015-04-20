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


f = open("2700.txt", "r")

t = float(f.readline())

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


uc = (u[1:-1, 2:-1] + u[1:-1, 1:-2])*0.5
vc = (v[2:-1, 1:-1] + v[1:-2, 1:-1])*0.5

p = p[:, :, 0].T

Nx = vN[0] - 2
Ny = uN[1] - 2
Lx = 2 * numpy.pi
Ly = 2 * numpy.pi
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

t = 0.4

u_e = u_ext(Xu, Yu, t)
v_e = v_ext(Xv, Yv, t)
p_e = p_ext(Xp, Yp, t)

uc_e = (u_e[:, 1:] + u_e[:, :-1]) / 2
vc_e = (v_e[1:, :] + v_e[:-1, :]) / 2


pyplot.figure()
fig = pyplot.contourf(Xu, Yu, u[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("u simulation")

pyplot.figure()
fig = pyplot.contourf(Xv, Yv, v[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("v simulation")

pyplot.figure()
fig = pyplot.contourf(Xp, Yp, p, 100)
pyplot.colorbar(fig)
pyplot.title("p simulation")


pyplot.figure()
fig = pyplot.streamplot(Xp, Yp, uc, vc,
                        density=3)
pyplot.axis("equal")
pyplot.xlim(0, 2*numpy.pi)
pyplot.ylim(0, 2*numpy.pi)

'''
pyplot.figure()
fig = pyplot.contourf(Xu, Yu, u_e, 100)
pyplot.colorbar(fig)
pyplot.title("u exact")

pyplot.figure()
fig = pyplot.contourf(Xv, Yv, v_e, 100)
pyplot.colorbar(fig)
pyplot.title("v exact")

pyplot.figure()
fig = pyplot.contourf(Xp, Yp, p_e, 100)
pyplot.colorbar(fig)
pyplot.title("p exact")



pyplot.figure()
fig = pyplot.contourf(Xu, Yu, u[1:-1, 1:-1] - u_e)
pyplot.colorbar(fig)

pyplot.figure()
fig = pyplot.contourf(Xv, Yv, v[1:-1, 1:-1] - v_e)
pyplot.colorbar(fig)
'''

pyplot.show()
