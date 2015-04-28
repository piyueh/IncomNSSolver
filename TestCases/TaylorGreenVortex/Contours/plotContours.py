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

File = "Data.txt"

f = open(File, "r")

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

u = u[:, :, 1].T
v = v[:, :, 1].T
p = p[:, :, 1].T

p -= numpy.average(p[1:-1, 1:-1])

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


pyplot.figure()
pyplot.title("u velocity @ T=2 sec", fontsize=16)
pyplot.xlabel("x")
pyplot.ylabel("y")
fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], u[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.savefig("TGvortexContour_u.png", format="png")

pyplot.figure()
pyplot.title("v velocity @ T=2 sec", fontsize=16)
pyplot.xlabel("x")
pyplot.ylabel("y")
fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], v[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.savefig("TGvortexContour_v.png", format="png")

pyplot.figure()
pyplot.title("p velocity @ T=2 sec", fontsize=16)
pyplot.xlabel("x")
pyplot.ylabel("y")
fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], p[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.savefig("TGvortexContour_p.png", format="png")


uc = (u[:, 1:] + u[:, :-1]) / 2
vc = (v[1:, :] + v[:-1, :]) / 2

uVor = (u[1:, 1:-1] - u[:-1, 1:-1]) / dy
vVor = (v[1:-1, 1:] - v[1:-1, :-1]) / dx
Vor = vVor - uVor

XVor, YVor = numpy.meshgrid(xu[1:-1], yv[1:-1]) 

pyplot.figure()
pyplot.title("Vorticity and Streamlines @ T=2 sec", fontsize=16)
pyplot.xlabel("x")
pyplot.ylabel("y")
fig = pyplot.contourf(XVor, YVor, Vor, 100)
pyplot.colorbar(fig)
pyplot.streamplot(Xp, Yp, uc, vc, density=2)
pyplot.xlim(0, 2*numpy.pi)
pyplot.ylim(0, 2*numpy.pi)
pyplot.savefig("TGvortexContour_VorStream.png", format="png")




pyplot.show()
