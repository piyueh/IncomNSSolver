import numpy
from matplotlib import pyplot

f = open("N512/Data.txt", "r")

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
Lx = Ly = 1
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

pyplot.figure()
fig = pyplot.contourf(Xu, Yu, u[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("u simulation at time = " + str(t) + " sec")

pyplot.figure()
fig = pyplot.contourf(Xv, Yv, v[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.title("v simulation at time = " + str(t) + " sec")

pyplot.figure()
fig = pyplot.contourf(Xp, Yp, p, 100)
pyplot.colorbar(fig)
pyplot.title("p simulation at time = " + str(t) + " sec")

pyplot.figure()
fig = pyplot.streamplot(Xp, Yp, uc, vc, density=3)
pyplot.title("stream lines at time = " + str(t) + " sec")
pyplot.axis("equal")
pyplot.xlim(0, Lx)
pyplot.ylim(0, Ly)
pyplot.show()
