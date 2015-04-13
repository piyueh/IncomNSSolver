# include "include/IncomNSSolver.h"

int Mesh::InitMesh(int N[3], double L[3])
{
	Nx = N[0]; Ny = N[1]; Nz = N[2];
	NCells = Nz * Ny * Nz;

	Lx = L[0]; Ly = L[1]; Lz = L[2];

	dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz;

	Nxu = Nx + 1; Nyu = Ny; Nzu = Nz;
	Nxv = Nx; Nyv = Ny + 1; Nzv = Nz;
	Nxw = Nx; Nyw = Ny; Nzw = Nz + 1;

	for(int i=0; i<Nx; ++i) xp.push_back((i+0.5)*dx);
	for(int j=0; j<Ny; ++j) yp.push_back((j+0.5)*dy);
	for(int k=0; k<Nz; ++k) zp.push_back((k+0.5)*dz);

	for(int i=0; i<Nxu; ++i) xu.push_back(i*dx);
	for(int j=0; j<Nyu; ++j) yu.push_back((j+0.5)*dy);
	for(int k=0; k<Nzu; ++k) yu.push_back((k+0.5)*dz);

	for(int i=0; i<Nxv; ++i) xv.push_back((i+0.5)*dx);
	for(int j=0; j<Nyv; ++j) yv.push_back(j*dy);
	for(int k=0; k<Nzv; ++k) yv.push_back((k+0.5)*dz);

	for(int i=0; i<Nxw; ++i) xw.push_back((i+0.5)*dx);
	for(int j=0; j<Nyw; ++j) yw.push_back((j+0.5)*dy);
	for(int k=0; k<Nzw; ++k) yw.push_back(k*dz);

	return 0;
}


int Mesh::addBC(string dir, 
		int pBCtype, double pValue, int vBCtype, double vValue)
{
	BCs[dir] = Boundary(Nx, Ny, Nx, dir, 
			pBCtype, pValue, vBCtype, vValue);
	return 0;
}
