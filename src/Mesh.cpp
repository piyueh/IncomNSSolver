# include "include/class_Mesh.h"

int Mesh::InitMesh(array<int, 3> N, array<double, 3> L)
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
	for(int k=0; k<Nzu; ++k) zu.push_back((k+0.5)*dz);

	for(int i=0; i<Nxv; ++i) xv.push_back((i+0.5)*dx);
	for(int j=0; j<Nyv; ++j) yv.push_back(j*dy);
	for(int k=0; k<Nzv; ++k) zv.push_back((k+0.5)*dz);

	for(int i=0; i<Nxw; ++i) xw.push_back((i+0.5)*dx);
	for(int j=0; j<Nyw; ++j) yw.push_back((j+0.5)*dy);
	for(int k=0; k<Nzw; ++k) zw.push_back(k*dz);

	return 0;
}


int Mesh::addBC(unsigned int dir, int sign, pair<int, double> p, 
		pair<int, double> u, pair<int, double> v, pair<int, double> w)
{
	BCs[dir*sign] = Boundary(Nx, Ny, Nz, dir, sign, p, u, v, w);
	return 0;
}

ostream &operator<<(ostream &os, Mesh &mesh)
{
	os << "Numbers of cells: "
	   << mesh.Nx << " x " << mesh.Ny << " x " << mesh.Nz << endl;
	os << "Shape of u array: "
	   << mesh.Nxu << " x " << mesh.Nyu << " x " << mesh.Nzu << endl;
	os << "Shape of v array: "
	   << mesh.Nxv << " x " << mesh.Nyv << " x " << mesh.Nzv << endl;
	os << "Shape of w array: "
	   << mesh.Nxw << " x " << mesh.Nyw << " x " << mesh.Nzw << endl;
	os << "Domain size: " 
	   << mesh.Lx << " x " << mesh.Ly << " x " << mesh.Lz << endl;
	os << "Cell size: " 
	   << mesh.dx << " x " << mesh.dy << " x " << mesh.dz << endl;
	os << endl;

	os << "Location of p point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xp << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yp << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zp << endl << endl;

	os << "Location of u point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xu << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yu << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zu << endl << endl;

	os << "Location of v point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xv << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yv << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zv << endl << endl;

	os << "Location of w point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xw << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yw << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zw << endl << endl;

	os << "Boundary conditions: " << endl;
	for(auto it=mesh.BCs.begin(); it!=mesh.BCs.end(); ++it)
		cout << it->second << endl;

	return os;
}
