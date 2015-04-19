/*
 * The class to store the variables of computation 
 * domain, and mesh
 */


# pragma once


class Mesh
{
	friend class NSSolver;
	friend ostream &operator<<(ostream &os, Mesh &mesh);

	public:

		Mesh() = default;
		Mesh(CaryI3 N, CaryD3 L){ InitMesh(N, L); };

		int InitMesh(CaryI3, CaryD3);

		int addBC(CUI dir, CI sign, CPairID p, CPairID u, CPairID v, CPairID w); 

		const map<int, Boundary> & get_BCs() const { return BCs; }

	private:

		int Nx, Ny, Nz;
		int NCells, Nyz;

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
