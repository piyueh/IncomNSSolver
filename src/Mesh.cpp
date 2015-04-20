# include "include/IncomNSSolver.h"


Mesh::Mesh(string &fName)
{
	ifstream file(fName);
	string line;

	while(getline(file, line))
	{
		istringstream wholeLine(line);
		string var;

		wholeLine >> var;

		if (var == "Nx") { wholeLine >> Nx; }
		else if (var == "Ny") { wholeLine >> Ny; }
		else if (var == "Nz") { wholeLine >> Nz; }
		else if (var == "Lx") { wholeLine >> Lx; }
		else if (var == "Ly") { wholeLine >> Ly; }
		else if (var == "Lz") { wholeLine >> Lz; }
		else if (var == "BC") {}
		else if (var.empty()) {}
		else
			throw invalid_argument(
					string("Invalid Argument in ") + fName + ": " + var);
	}
	file.close();

	InitMesh({Nx, Ny, Nz}, {Lx, Ly, Lz});

	file.open(fName);
	while(getline(file, line))
	{
		istringstream wholeLine(line);
		string var;

		wholeLine >> var;

		if (var == "BC")
		{
			unsigned int dir;
			wholeLine >> dir;

			int sign;
			wholeLine >> sign;

			pair<int, double> p, u, v, w;
			wholeLine >> p.first >> p.second;
			wholeLine >> u.first >> u.second;
			wholeLine >> v.first >> v.second;
			wholeLine >> w.first >> w.second;

			addBC(dir, sign, p, u, v, w);

		}
	}
	file.close();

}


int Mesh::InitMesh(array<int, 3> N, array<double, 3> L)
{
	Nx = N[0]; Ny = N[1]; Nz = N[2];
	NCells = Nz * Ny * Nz; Nyz = Ny * Nz;

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
	BCs[dir*sign] = Boundary({Nx, Ny, Nz}, dir, sign, p, u, v, w);
	return 0;
}

