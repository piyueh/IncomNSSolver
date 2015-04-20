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

		const int get_Nx() const { return Nx; };
		const int get_Ny() const { return Ny; };
		const int get_Nz() const { return Nz; };

		const int get_Nxu() const { return Nxu; };
		const int get_Nyu() const { return Nyu; };
		const int get_Nzu() const { return Nzu; };

		const int get_Nxv() const { return Nxv; };
		const int get_Nyv() const { return Nyv; };
		const int get_Nzv() const { return Nzv; };

		const int get_Nxw() const { return Nxw; };
		const int get_Nyw() const { return Nyw; };
		const int get_Nzw() const { return Nzw; };

		const double get_Lx() const { return Lx; };
		const double get_Ly() const { return Ly; };
		const double get_Lz() const { return Lz; };

		const double get_dx() const { return dx; };
		const double get_dy() const { return dy; };
		const double get_dz() const { return dz; };

		map<int, Boundary> & get_BCs() { return BCs; }

		const vector<double> & get_xp() const { return xp; }
		const vector<double> & get_yp() const { return yp; }
		const vector<double> & get_zp() const { return zp; }

		const vector<double> & get_xu() const { return xu; }
		const vector<double> & get_yu() const { return yu; }
		const vector<double> & get_zu() const { return zu; }

		const vector<double> & get_xv() const { return xv; }
		const vector<double> & get_yv() const { return yv; }
		const vector<double> & get_zv() const { return zv; }

		const vector<double> & get_xw() const { return xw; }
		const vector<double> & get_yw() const { return yw; }
		const vector<double> & get_zw() const { return zw; }

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
