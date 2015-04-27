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


def updGhosts(u, v, p):
    u[0, :] = u[2, :]
    u[-1, :] = u[-3, :]
    u[:, 0] = - u[:, 1]
    u[:, -1] = - u[:, -2]

    v[0, :] = - v[1, :]
    v[-1, :] = - v[-2, :]
    v[:, 0] = v[:, 2]
    v[:, -1] = v[:, -3]

    p[0, :] = p[1, :]
    p[-1, :] = p[-2, :]
    p[:, 0] = p[:, 1]
    p[:, -1] = p[:, -2]


def CalG(nu, u, v, p, dp, dx, dy, Gu, Gv, coef):
    '''
    '''
    Gu[1:-1, 1:-1] = (Gu[1:-1, 1:-1] - (dp[1:, 1:-1] - dp[:-1, 1:-1]) / dx) * coef + \
        nu * ((u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx / dx \
            + (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2]) / dy / dy) \
      - 0.25 * ((u[2:, 1:-1] + u[1:-1, 1:-1])**2 -
                (u[1:-1, 1:-1] + u[:-2, 1:-1])**2) / dx \
      - 0.25 * ((u[1:-1, 2:] + u[1:-1, 1:-1])*(v[:-1, 2:-1]+v[1:, 2:-1]) -
                (u[1:-1, :-2] + u[1:-1, 1:-1])*(v[:-1, 1:-2]+v[1:, 1:-2])) / dy \
      - (p[1:, 1:-1] - p[:-1, 1:-1]) / dx

    Gv[1:-1, 1:-1] = (Gv[1:-1, 1:-1] - (dp[1:-1, 1:] - dp[1:-1, :-1]) / dy) * coef + \
        nu * ((v[2:, 1:-1] - 2 * v[1:-1, 1:-1] + v[:-2, 1:-1]) / dx / dx \
            + (v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, :-2]) / dy / dy) \
      - 0.25 * ((u[2:-1, 1:] + u[2:-1, :-1])*(v[2:, 1:-1]+v[1:-1,1:-1]) -
                (u[1:-2, 1:] + u[1:-2, :-1])*(v[1:-1, 1:-1]+v[:-2,1:-1])) / dx \
      - 0.25 * ((v[1:-1, 2:] + v[1:-1, 1:-1])*(v[1:-1, 2:] + v[1:-1, 1:-1]) -
                (v[1:-1, :-2]+ v[1:-1, 1:-1])*(v[1:-1, :-2]+v[1:-1, 1:-1])) / dy \
      - (p[1:-1, 1:] - p[1:-1, :-1]) / dy


def singleEuler(nu, u, v, p, dp, Gu, Gv, dx, dy, dt, Nx, Ny, c1, c2):
    '''
    '''

    updGhosts(u, v, p)

    CalG(nu, u, v, p, dp, dx, dy, Gu, Gv, c1)

    u += c2 * dt * Gu
    v += c2 * dt * Gv

    b = numpy.zeros((Nx+2, Ny+2))
    b[1:-1, 1:-1] = (u[2:-1, 1:-1] - u[1:-2, 1:-1]) / dx + \
                    (v[1:-1, 2:-1] - v[1:-1, 1:-2]) / dy

    dp, Nitr = CG(Nx, Ny, numpy.zeros((Nx+2, Ny+2)),
                               b, dx, dy, 1e-15, 'N', refP=0)

    u[2:-2, 1:-1] -= (dp[2:-1, 1:-1] - dp[1:-2, 1:-1]) / dx
    v[1:-1, 2:-2] -= (dp[1:-1, 2:-1] - dp[1:-1, 1:-2]) / dy

    dp /= (c2 * dt)
    p += dp

    return u, v, p, dp, Gu, Gv


def RK3_NSSolver(nu, u, v, p, dx, dy, dt, Nx, Ny, Nt):
    '''
    '''
    t = 0

    Gu = numpy.zeros((Nx+3, Ny+2))
    Gv = numpy.zeros((Nx+2, Ny+3))
    dp = numpy.zeros_like(p)

    for n in range(Nt):

        u, v, p, dp, Gu, Gv = \
            singleEuler(nu, u, v, p, dp, Gu, Gv, dx, dy, dt, Nx, Ny, 0, 1./3.)
        u, v, p, dp, Gu, Gv = \
            singleEuler(nu, u, v, p, dp, Gu, Gv, dx, dy, dt, Nx, Ny, -5./9., 15./16.)
        u, v, p, dp, Gu, Gv = \
            singleEuler(nu, u, v, p, dp, Gu, Gv, dx, dy, dt, Nx, Ny, -153./128., 8./15.)

        t = (n + 1) * dt
        print("N = ", n+1, " Time = ", t)


    updGhosts(u, v, p)

    return t, u, v, p



Nx = Ny = 50
Lx = Ly = 2 * numpy.pi
dx = Lx / Nx
dy = Ly / Ny


Xu, Yu, Xv, Yv, Xp, Yp = genCoord(Lx, Ly, dx, dy, Nx, Ny)


nu = 1
DT = numpy.array([0.005 / (2**n) for n in range(7)])

ErrU = numpy.zeros(0)
ErrV = numpy.zeros(0)
ErrP = numpy.zeros(0)


for dt in DT:

    Nt = int(0.005 / dt + 0.5)

    u = u_ext(Xu, Yu, 0)
    v = v_ext(Xv, Yv, 0)
    p = v_ext(Xp, Yp, 0)

    t, u, v, p = RK3_NSSolver(nu, u, v, p, dx, dy, dt, Nx, Ny, Nt)

    uc = (u[2:-1, 1:-1] + u[1:-2, 1:-1]) / 2
    vc = (v[1:-1, 2:-1] + v[1:-1, 1:-2]) / 2

    ErrU = numpy.append(ErrU, numpy.abs(u[1:-1, 1:-1] - u_ext(Xu, Yu, t)[1:-1, 1:-1]).max())
    ErrV = numpy.append(ErrV, numpy.abs(v[1:-1, 1:-1] - v_ext(Xv, Yv, t)[1:-1, 1:-1]).max())
    ErrP = numpy.append(ErrP, numpy.abs(p[1:-1, 1:-1] + p_ext(Xp[1, 1], Yp[1, 1], t) - p_ext(Xp, Yp, t)[1:-1, 1:-1]).max())
    '''
    ErrU = numpy.append(ErrU, u[1:-1, 1:-1].max())
    ErrV = numpy.append(ErrV, v[1:-1, 1:-1].max())
    ErrP = numpy.append(ErrP, p[1:-1, 1:-1].max())
    '''

    '''
    pyplot.figure()
    #fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], (u - u_ext(Xu, Yu, t))[1:-1, 1:-1], 100)
    fig = pyplot.contourf(Xu[1:-1, 1:-1], Yu[1:-1, 1:-1], u[1:-1, 1:-1], 100)
    pyplot.colorbar(fig)

    pyplot.figure()
    #fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], (v - v_ext(Xv, Yv, t))[1:-1, 1:-1], 100)
    fig = pyplot.contourf(Xv[1:-1, 1:-1], Yv[1:-1, 1:-1], v[1:-1, 1:-1], 100)
    pyplot.colorbar(fig)

    pyplot.figure()
    #fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], (p - p_ext(Xp, Yp, t))[1:-1, 1:-1], 100)
    fig = pyplot.contourf(Xp[1:-1, 1:-1], Yp[1:-1, 1:-1], p[1:-1, 1:-1] + p_ext(Xp[1, 1], Yp[1, 1], t), 100)
    pyplot.colorbar(fig)

    pyplot.show()
    '''


'''
ErrU = numpy.abs(ErrU-ErrU[-1])
ErrV = numpy.abs(ErrV-ErrV[-1])
ErrP = numpy.abs(ErrP-ErrP[-1])
'''


for i, dt in enumerate(DT):
    print(dt, ErrU[i], ErrV[i], ErrP[i])

pyplot.loglog(DT[:-1], ErrU[:-1], 'k.-', lw=2, label="Err(u) velocity")
pyplot.loglog(DT[:-1], ErrV[:-1], 'kx-', lw=2, label="Err(v) velocity")
pyplot.loglog(DT[:-1], ErrP[:-1], 'kd-', lw=2, label="Err(p) velocity")
pyplot.loglog(DT, numpy.array([ErrU[0]/(2**n) for n in range(DT.size)]), 'r--', lw=2, label="1st order")
pyplot.loglog(DT, numpy.array([ErrU[0]/(4**n) for n in range(DT.size)]), 'g--', lw=2, label="2nd order")
pyplot.loglog(DT, numpy.array([ErrU[0]/(8**n) for n in range(DT.size)]), 'b--', lw=2, label="3rd order")
pyplot.axis("equal")
pyplot.legend(loc=0)
pyplot.show()
