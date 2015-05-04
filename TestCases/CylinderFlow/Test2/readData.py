import numpy
from matplotlib import pyplot

f = open("Data.txt", "r")

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


u = u[:, :, 1].T
v = v[:, :, 1].T
p = p[:, :, 1].T

uc = (u[1:-1, 2:-1] + u[1:-1, 1:-2]) * 0.5
vc = (v[2:-1, 1:-1] + v[1:-2, 1:-1]) * 0.5

Nx = vN[0] - 2
Ny = uN[1] - 2
Lx = 40
Ly = 20
dx = Lx / Nx
dy = Ly / Ny

uD = (u[1:-1, 2:-1] - u[1:-1, 1:-2]) / dx
vD = (v[2:-1, 1:-1] - v[1:-2, 1:-1]) / dy
Div = uD + vD

xp = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
yp = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
Xp, Yp = numpy.meshgrid(xp, yp)


xu = numpy.linspace(-dx, Lx+dx, Nx+3)
yu = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
Xu, Yu = numpy.meshgrid(xu, yu)


xv = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
yv = numpy.linspace(-dy, Ly+dy, Ny+3)
Xv, Yv = numpy.meshgrid(xv, yv)

th = numpy.linspace(0, 2 * numpy.pi, 360)
xCirc = 0.5 * numpy.cos(th) + 10
yCirc = 0.5 * numpy.sin(th) + 10

pyplot.figure()
pyplot.title("Streamlines @ T=" + str(t) + "2 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.streamplot(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], uc, vc, density=4)
pyplot.fill(xCirc, yCirc, fc='w', ec='k')
pyplot.axis("equal")
pyplot.xlim(0, Lx)
pyplot.ylim(0, Ly)
#pyplot.savefig("Lid_Streamlines.png", format="png")

pyplot.figure()
pyplot.title("u velocity @ T=" + str(t) + "2 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], u[1:-1, 1:-1],
                      extend="both", levels=numpy.linspace(-0, 1.5, 100))
pyplot.colorbar(fig)
pyplot.fill(xCirc, yCirc, fc='w', ec='k')
#pyplot.savefig("Lid_uVelocity.png", format="png")

pyplot.figure()
pyplot.title("v velocity @ T=" + str(t) + "2 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], v[1:-1, 1:-1],
                      extend="both", levels=numpy.linspace(-1, 1, 100))
pyplot.colorbar(fig)
pyplot.fill(xCirc, yCirc, fc='w', ec='k')
#pyplot.savefig("Lid_vVelocity.png", format="png")

pyplot.figure()
pyplot.title("Pressure @ T=" + str(t) + "2 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], p[1:-1, 1:-1],
                      extend="both", levels=numpy.linspace(-0.5, 0.5, 100))
pyplot.colorbar(fig)
pyplot.fill(xCirc, yCirc, fc='w', ec='k')
#pyplot.savefig("Lid_Pressure.png", format="png")


uVor = (u[1:, 1:-1] - u[:-1, 1:-1]) / dy
vVor = (v[1:-1, 1:] - v[1:-1, :-1]) / dx
Vor = vVor - uVor

XVor, YVor = numpy.meshgrid(xu[1:-1], yv[1:-1])

pyplot.figure()
pyplot.title("Vorticity @ T=" + str(t) + "2 sec", fontsize=16)
pyplot.xlabel("x")
pyplot.ylabel("y")
fig = pyplot.contourf(XVor, YVor, Vor,
                      extend="both", levels=numpy.linspace(-2, 2, 100))
#pyplot.scatter(XVor, YVor)
pyplot.fill(xCirc, yCirc, fc='w', ec='k')
pyplot.colorbar(fig)
#pyplot.savefig("TGvortexContour_VorStream.png", format="png")


pyplot.show()

