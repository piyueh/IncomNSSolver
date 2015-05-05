# IncomNSSolver
=================

## Current Feature

* 3D incompressible N-S solver
* 2nd-order finite difference method for spatial discretization
* 3rd-order RK method for time marching
* Staggered grid
* two different Poisson solver:
	1. biconjugate gradient method with incomplete LUT preconditioner (single
core code)
	2. Parallel LDLT decomposition using Intel Pardiso library

## Special version

* Cylinder flow using immersed boundary method (hard-coded, use macro
variable `CYLINDER_FLOW` to control it)

## Future Feature

* Convecddtive BC
