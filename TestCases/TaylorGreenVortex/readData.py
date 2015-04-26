import numpy
from genTGinitData import u_ext as u_ext
from genTGinitData import v_ext as v_ext
from genTGinitData import p_ext as p_ext
from matplotlib import pyplot

f = open("TimeConvergeTest/Data.txt", "r")

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


uc_e = (u_e[:, 1:] + u_e[:, :-1]) / 2
vc_e = (v_e[1:, :] + v_e[:-1, :]) / 2








'''
pyplot.figure()
fig = pyplot.streamplot(Xp, Yp, uc, vc,
                        density=3)
pyplot.axis("equal")
pyplot.xlim(0, 2*numpy.pi)
pyplot.ylim(0, 2*numpy.pi)
'''

pyplot.figure()
fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], u[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("u simulation")

pyplot.figure()
fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], u_e[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("u exact")


pyplot.figure()
fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], v[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("v simulation")

pyplot.figure()
fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], v_e[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("v exact")


pyplot.figure()
fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], p[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("p simulation")

pyplot.figure()
fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], p_e[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("p exact")



pyplot.figure()
fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1],
                      (u - u_e)[1:-1, 1:-1], 100)
pyplot.colorbar(fig)

pyplot.figure()
fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1],
                      (v - v_e)[1:-1, 1:-1], 100)
pyplot.colorbar(fig)

pyplot.figure()
fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], (p - p_e)[1:-1, 1:-1], 100)
pyplot.colorbar(fig)

pyplot.show()
