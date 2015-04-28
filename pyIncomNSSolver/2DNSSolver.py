# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 23:40:27 2015

@author: pychuang
"""
import numpy
from matplotlib import pyplot
from scipy import sparse
from CG import *
import sys
import getopt


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


def genCoord(Lx, Ly, dx, dy, Nx, Ny):
    '''
    '''
    xu = numpy.linspace(-dx, Lx+dx, Nx+3)
    yu = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
    Xu, Yu = numpy.meshgrid(xu, yu)

    xv = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
    yv = numpy.linspace(-dy, Ly+dy, Ny+3)
    Xv, Yv = numpy.meshgrid(xv, yv)

    xp = numpy.linspace(-dx/2, Lx+dx/2, Nx+2)
    yp = numpy.linspace(-dy/2, Ly+dy/2, Ny+2)
    Xp, Yp = numpy.meshgrid(xp, yp)

    return Xu.T, Yu.T, Xv.T, Yv.T, Xp.T, Yp.T, 


def updGhosts(u, v):
    u[0, :] = u[2, :]
    u[-1, :] = u[-3, :]
    u[:, 0] = - u[:, 1]
    u[:, -1] = - u[:, -2]

    v[0, :] = - v[1, :]
    v[-1, :] = - v[-2, :]
    v[:, 0] = v[:, 2]
    v[:, -1] = v[:, -3]


def CalG(nu, u, v, dx, dy, Gu, Gv, coef):
    '''
    '''
    Gu[1:-1, 1:-1] = Gu[1:-1, 1:-1] * coef + \
        nu * ((u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx / dx \
            + (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2]) / dy / dy) \
      - 0.25 * ((u[2:, 1:-1] + u[1:-1, 1:-1])**2 -
                (u[1:-1, 1:-1] + u[:-2, 1:-1])**2) / dx \
      - 0.25 * ((u[1:-1, 2:] + u[1:-1, 1:-1])*(v[:-1, 2:-1]+v[1:, 2:-1]) -
                (u[1:-1, :-2] + u[1:-1, 1:-1])*(v[:-1, 1:-2]+v[1:, 1:-2])) / dy

    Gv[1:-1, 1:-1] = Gv[1:-1, 1:-1] * coef + \
        nu * ((v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[:-2, 1:-1]) / dx / dx \
            + (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, :-2]) / dy / dy) \
      - 0.25 * ((u[2:-1, 1:] + u[2:-1, :-1])*(v[2:, 1:-1]+v[1:-1,1:-1]) -
                (u[1:-2, 1:] + u[1:-2, :-1])*(v[1:-1, 1:-1]+v[:-2,1:-1])) / dx \
      - 0.25 * ((v[1:-1, 2:] + v[1:-1, 1:-1])*(v[1:-1, 2:] + v[1:-1, 1:-1]) -
                (v[1:-1, :-2]+ v[1:-1, 1:-1])*(v[1:-1, :-2]+v[1:-1, 1:-1])) / dy


def singleEuler(nu, u, v, Gu, Gv, dx, dy, dt, Nx, Ny, c1, c2):
    '''
    '''

    updGhosts(u, v)

    CalG(nu, u, v, dx, dy, Gu, Gv, c1)

    u += c2 * dt * Gu
    v += c2 * dt * Gv

    b = numpy.zeros((Nx+2, Ny+2))
    b[1:-1, 1:-1] = (u[2:-1, 1:-1] - u[1:-2, 1:-1]) / dx + \
                    (v[1:-1, 2:-1] - v[1:-1, 1:-2]) / dy

    p, Nitr = CG(Nx, Ny, numpy.zeros((Nx+2, Ny+2)),
                               b, dx, dy, 1e-15, 'N', refP=0)

    u[2:-2, 1:-1] -= (p[2:-1, 1:-1] - p[1:-2, 1:-1]) / dx
    v[1:-1, 2:-2] -= (p[1:-1, 2:-1] - p[1:-1, 1:-2]) / dy

    return u, v, p, Gu, Gv


def RK3_NSSolver(nu, u, v, dx, dy, dt, Nx, Ny, Nt):
    '''
    '''
    t = 0

    Gu = numpy.zeros((Nx+3, Ny+2))
    Gv = numpy.zeros((Nx+2, Ny+3))

    for n in range(Nt):

        u, v, p, Gu, Gv = \
            singleEuler(nu, u, v, Gu, Gv, dx, dy, dt, Nx, Ny, 0, 1.0/3.0)
        u, v, p, Gu, Gv = \
            singleEuler(nu, u, v, Gu, Gv, dx, dy, dt, Nx, Ny, -5./9., 15./16.)
        u, v, p, Gu, Gv = \
            singleEuler(nu, u, v, Gu, Gv, dx, dy, dt, Nx, Ny, -153./128., 8./15.)

        t = (n + 1) * dt
        print("N = ", n, " Time = ", t)


    updGhosts(u, v)

    return t, u, v, p



Nx = Ny = 50
Lx = Ly = 2 * numpy.pi
dx = Lx / Nx
dy = Ly / Ny


Xu, Yu, Xv, Yv, Xp, Yp = genCoord(Lx, Ly, dx, dy, Nx, Ny)


nu = 1
DT = numpy.array([0.005 / (2**n) for n in range(8)])
ErrU = numpy.zeros(0)
ErrV = numpy.zeros(0)


for dt in DT:

    Nt = int(0.1 / dt + 0.5)

    u = u_ext(Xu, Yu, 0)
    v = v_ext(Xv, Yv, 0)

    t, u, v, p = RK3_NSSolver(nu, u, v, dx, dy, dt, Nx, Ny, Nt)

    uc = (u[2:-1, 1:-1] + u[1:-2, 1:-1]) / 2
    vc = (v[1:-1, 2:-1] + v[1:-1, 1:-2]) / 2

    ErrU = numpy.append(ErrU, u[1:-1, 1:-1].max())
    ErrV = numpy.append(ErrV, v[1:-1, 1:-1].max())


ErrU -= ErrU[-1]
ErrV -= ErrV[-1]
ErrU = numpy.abs(ErrU)
ErrV = numpy.abs(ErrV)


pyplot.loglog(DT[:-1], ErrU[:-1], 'k.-', lw=2, label="Err(u) velocity")
pyplot.loglog(DT[:-1], ErrU[:-1], 'kx-', lw=2, label="Err(v) velocity")
pyplot.loglog(DT[:-1], numpy.array([ErrU[0]/(2**n) for n in range(DT.size-1)]), 'r--', lw=2, label="1st order")
pyplot.loglog(DT[:-1], numpy.array([ErrU[0]/(4**n) for n in range(DT.size-1)]), 'g--', lw=2, label="2nd order")
pyplot.loglog(DT[:-1], numpy.array([ErrU[0]/(8**n) for n in range(DT.size-1)]), 'b--', lw=2, label="3rd order")
pyplot.axis("equal")
pyplot.legend(loc=0)
pyplot.show()
