/*
 * The class to store the variables of computation 
 * domain, and mesh
 */


# pragma once


class Mesh
{
	friend class NSSolverEuler;
	friend ostream &operator<<(ostream &os, Mesh &mesh);

	public:

		Mesh() = default;
		Mesh(array<int, 3> N, array<double, 3> L){ InitMesh(N, L); };

		int InitMesh(array<int, 3>, array<double, 3>);

		int addBC(unsigned int dir, int sign, pair<int, double> p, 
				pair<int, double> u, pair<int, double> v, pair<int, double> w); 

	private:

		int Nx, Ny, Nz;
		int NCells;

		int Nxu, Nyu, Nzu;
		int Nxv, Nyv, Nzv;
		int Nxw, Nyw, Nzw;

		double Lx, Ly, Lz;
		double dx, dy, dz;

		map<int, Boundary> BCs;	
		
		vector<double> xp, yp, zp;
		vector<double> xu, yu, zu;
		vector<double> xv, yv, zv;
		vector<double> xw, yw, zw;
};
