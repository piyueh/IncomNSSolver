# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 23:40:27 2015

@author: pychuang
"""

import numpy
from matplotlib import pyplot
from scipy import sparse


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


def TrivialFunc(*trivial):
    pass

def residual(u, dx, dy, f):
    '''
    Calculate residual.
    '''
    r = numpy.zeros_like(u)
    r[1:-1, 1:-1] = f[1:-1, 1:-1] - \
            (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / dy / dy - \
            (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2]) / dx / dx
    return r

def set_NeumannBC(u, refP):
    '''
    Adiabatic BC. refP is always at u[0, 1]
    '''
    u[0, 1:-1] = u[1, 1:-1]
    u[-1, 1:-1] = u[-2, 1:-1]
    u[1:-1, 0] = u[1:-1, 1]
    u[1:-1, -1] = u[1:-1, -2]
    u[0, 1] = refP

def L2norm(x):
    '''
    Return L2-norm of a vector x.
    '''
    return numpy.sqrt(numpy.sum(x**2))


def CG(Nx, Ny, u, f, dx, dy, tol, BCtype, refP=None, **opt):
    '''
    Conjugate gradient method. A general linear solver for problem Ax=b.
    '''
    assert dx == dy

    BCHandle = {'N': set_NeumannBC, 'D': TrivialFunc}

    A = sparse.diags([-4]*Nx*Ny, 0)
    A += sparse.diags([1]*(Nx-1)+([0]+[1]*(Nx-1))*(Ny-1), 1)
    A += sparse.diags([1]*(Nx-1)+([0]+[1]*(Nx-1))*(Ny-1), -1)
    A += sparse.diags([1]*(Nx*Ny-Nx), Nx)
    A += sparse.diags([1]*(Nx*Ny-Nx), -Nx)

    if BCtype == 'D':
        b = dx * dx * f[1:-1, 1:-1].reshape(Nx * Ny)
        print("b = \n")
        print(b)
        b[:Nx] += -u[0, 1:-1]
        b[::Nx] += -u[1:-1, 0]
        b[Nx-1::Nx] += -u[1:-1, -1]
        b[-Nx:] += -u[-1, 1:-1]
    elif BCtype == 'N':
        b = dx * dx * f[1:-1, 1:-1].reshape(Nx * Ny)
        print("b = \n")
        print(b)
        b[0] += - refP
        A[range(1, Nx), range(1, Nx)] += 1
        A[range(-Nx, 0), range(-Nx, 0)] += 1
        A[range(0, Nx*Ny, Nx), range(0, Nx*Ny, Nx)] += 1
        A[range(Nx-1, Nx*Ny, Nx), range(Nx-1, Nx*Ny, Nx)] += 1

    x = u[1:-1, 1:-1].reshape(Nx * Ny)

    err = 1000.
    itr_count = 0
    r = b - A.dot(x)
    p = r
    while (err > tol and itr_count < int(1e8)):

        itr_count += 1

        aph = r.dot(r) / p.dot(A.dot(p))
        aph *= p
        x += aph

        # aph currently is $aph * p = x^{k+1} - x^k$
        err = L2norm(aph)

        r_new = b - A.dot(x)
        beta = r_new.dot(r_new) / r.dot(r)
        p = r_new + beta * p
        r = r_new

        if itr_count % 50 == 0:
            print(itr_count, err, aph.max(), aph.min())

    u[1:-1, 1:-1] = x.reshape((Ny, Nx))

    BCHandle[BCtype](u, refP)

    return u, itr_count



f = open("G.txt", "r")

GuN = numpy.array([int(x) for x in f.readline().split()])
Gu = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(GuN))

GvN = numpy.array([int(x) for x in f.readline().split()])
Gv = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(GvN))

GwN = numpy.array([int(x) for x in f.readline().split()])
Gw = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(GwN))

f.close()



b = numpy.genfromtxt("b.txt")

f = open("Data.txt", "r")

t = float(f.readline())
t = 0

uN = numpy.array([int(x) for x in f.readline().split()])
u = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(uN))

vN = numpy.array([int(x) for x in f.readline().split()])
v = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(vN))

wN = numpy.array([int(x) for x in f.readline().split()])
w = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(wN))

pN = numpy.array([int(x) for x in f.readline().split()])
p = numpy.array([float(x) for x in f.readline().split()]).reshape(tuple(pN))

f.close()

Nx = vN[0] - 2
Ny = uN[1] - 2
Lx = Ly = 2 * numpy.pi
dx = Lx / Nx
dy = Ly / Ny

b = b.reshape((Nx, Ny))

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

p_e, u_e, v_e, w_e  = TaylorGreen(t, Nx, Ny, Nz=1)
u_e = u_e[:, :, 1]
v_e = v_e[:, :, 1]

Gu_e = numpy.zeros_like(u[:, :, 0])
Gv_e = numpy.zeros_like(v[:, :, 0])

dt = 0.0005

Gu_e[1:-1, 1:-1] = \
    (u_e[2:, 1:-1] - 2 * u_e[1:-1, 1:-1] + u_e[:-2, 1:-1]) / dx / dx \
  + (u_e[1:-1, 2:] - 2 * u_e[1:-1, 1:-1] + u_e[1:-1, :-2]) / dy / dy \
  - 0.25 * ((u_e[2:, 1:-1] + u_e[1:-1, 1:-1])**2 -
            (u_e[1:-1, 1:-1] + u_e[:-2, 1:-1])**2) / dx \
  - 0.25 * ((u_e[1:-1, 2:] + u_e[1:-1, 1:-1])*(v_e[:-1, 2:-1]+v_e[1:, 2:-1]) -
            (u_e[1:-1, :-2] + u_e[1:-1, 1:-1])*(v_e[:-1, 1:-2]+v_e[1:, 1:-2])) / dy


Gv_e[1:-1, 1:-1] = \
    (v_e[2:, 1:-1] - 2 * v_e[1:-1, 1:-1] + v_e[:-2, 1:-1]) / dx / dx \
  + (v_e[1:-1, 2:] - 2 * v_e[1:-1, 1:-1] + v_e[1:-1, :-2]) / dy / dy \
  - 0.25 * ((u_e[2:-1, 1:] + u_e[2:-1, :-1])*(v_e[2:, 1:-1]+v_e[1:-1,1:-1]) -
            (u_e[1:-2, 1:] + u_e[1:-2, :-1])*(v_e[1:-1, 1:-1]+v_e[:-2,1:-1])) / dx \
  - 0.25 * ((v_e[1:-1, 2:] + v_e[1:-1, 1:-1])*(v_e[1:-1, 2:] + v_e[1:-1, 1:-1]) -
            (v_e[1:-1, :-2]+ v_e[1:-1, 1:-1])*(v_e[1:-1, :-2]+v_e[1:-1, 1:-1])) / dy

u_e += dt * Gu_e
v_e += dt * Gv_e

b_e = numpy.zeros((Nx+2, Ny+2))
b_e[1:-1, 1:-1] = (u_e[2:-1, 1:-1] - u_e[1:-2, 1:-1]) / dx + \
      (v_e[1:-1, 2:-1] - v_e[1:-1, 1:-2]) / dy
b_e /= dt

p_e, Nitr = CG(Nx, Ny, numpy.zeros((Nx+2, Ny+2)),
               b_e, dx, dy, 1e-10, 'N', refP=0)

u_e[2:-2, 1:-1] -= dt * (p_e[2:-1, 1:-1] - p_e[1:-2, 1:-1]) / dx
v_e[1:-1, 2:-2] -= dt * (p_e[1:-1, 2:-1] - p_e[1:-1, 1:-2]) / dy





pyplot.figure()
fig = pyplot.contourf(Xu, Yu, u_e - u[:, :, 1], 100)
pyplot.colorbar(fig)

pyplot.figure()
fig = pyplot.contourf(Xv, Yv, v_e - v[:, :, 1], 100)
pyplot.colorbar(fig)

pyplot.figure()
fig = pyplot.contourf(Xp, Yp, p_e[1:-1, 1:-1], 100)
pyplot.colorbar(fig)

