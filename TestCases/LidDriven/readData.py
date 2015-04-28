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
Lx = Ly = 1
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



yGhia = numpy.array([0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 
                     0.2813, 0.4531, 0.5, 0.6172, 0.7344, 0.8516, 
                     0.9531, 0.9609, 0.9688, 0.9766, 1])
uGhia = numpy.array([0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, 
                     -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 
                     0.68717, 0.73722, 0.78871, 0.84123, 1])

xGhia = numpy.array([0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 
                     0.2266, 0.2344, 0.5, 0.8047, 0.8594, 0.9063, 
                     0.9453, 0.9531, 0.9609, 0.9688, 1])
vGhia = numpy.array([0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 
                     0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914, 
                     -0.10313, -0.08864, -0.07391, -0.05906, 0])


pyplot.figure()
pyplot.title("Streamlines @ t = 50 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.streamplot(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], uc, vc, density=2)
pyplot.axis("equal")
pyplot.xlim(0, Lx)
pyplot.ylim(0, Ly)
pyplot.savefig("Lid_Streamlines.png", format="png")

pyplot.figure()
pyplot.title("u velocity @ t = 50 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], u[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.savefig("Lid_uVelocity.png", format="png")

pyplot.figure()
pyplot.title("v velocity @ t = 50 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], v[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.savefig("Lid_vVelocity.png", format="png")

pyplot.figure()
pyplot.title("Pressure @ t = 50 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$y / L$", fontsize=18)
fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], p[1:-1, 1:-1], 100)
pyplot.colorbar(fig)
pyplot.savefig("Lid_Pressure.png", format="png")

pyplot.figure()
pyplot.title("u velocity along the vertical center line" +
             "\n @ t = 50 sec", fontsize=18)
pyplot.xlabel(r"$y / L$", fontsize=18)
pyplot.ylabel(r"$u / U_{lid}$", fontsize=18)
pyplot.plot(u[1:-1, Nx/2+1], Yu[1:-1, Nx/2+1], 'k-', lw=2, label="present project")
pyplot.plot(uGhia, yGhia, 'r^', markersize=10, label="Ghia et. al., 1982")
pyplot.grid(True)
pyplot.legend(loc=0)
pyplot.savefig("Lid_CenterUvelocity.png", format="png")

pyplot.figure()
pyplot.title("v velocity along the horizontal center line" +
             "\n @ t = 50 sec", fontsize=18)
pyplot.xlabel(r"$x / L$", fontsize=18)
pyplot.ylabel(r"$v / U_{lid}$", fontsize=18)
pyplot.plot(Xv[Ny/2+1, 1:-1], v[Ny/2+1, 1:-1], 'k-', lw=2, label="present project")
pyplot.plot(xGhia, vGhia, 'r^', lw=2, markersize=10, label="Ghia et. al., 1982")
pyplot.grid(True)
pyplot.legend(loc=0)
pyplot.savefig("Lid_CenterVvelocity.png", format="png")


pyplot.show()

