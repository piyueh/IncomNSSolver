# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 23:40:27 2015

@author: pychuang
"""
import numpy
from scipy import sparse
import sys


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
        b[:Nx] += -u[0, 1:-1]
        b[::Nx] += -u[1:-1, 0]
        b[Nx-1::Nx] += -u[1:-1, -1]
        b[-Nx:] += -u[-1, 1:-1]
    elif BCtype == 'N':
        b = dx * dx * f[1:-1, 1:-1].reshape(Nx * Ny)
        b[0] += - refP
        A.todense()[range(1, Nx), range(1, Nx)] += 1
        A.todense()[range(-Nx, 0), range(-Nx, 0)] += 1
        A.todense()[range(0, Nx*Ny, Nx), range(0, Nx*Ny, Nx)] += 1
        A.todense()[range(Nx-1, Nx*Ny, Nx), range(Nx-1, Nx*Ny, Nx)] += 1

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

        '''
        if itr_count % 50 == 0:
            print(itr_count, err, aph.max(), aph.min())
        '''

    if itr_count >= int(1e8): 
        print("Poisson solver diverged.")
        sys.exit(2)

    u[1:-1, 1:-1] = x.reshape((Ny, Nx))

    BCHandle[BCtype](u, refP)

    return u, itr_count
