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
		Mesh(string &);
		Mesh(CaryI3 N, CaryD3 L){ InitMesh(N, L); };

		int InitMesh(CaryI3, CaryD3);

		int addBC(CUI dir, CI sign, CPairID p, CPairID u, CPairID v, CPairID w); 

		CI get_Nx() const { return Nx; };
		CI get_Ny() const { return Ny; };
		CI get_Nz() const { return Nz; };

		CI get_Nxu() const { return Nxu; };
		CI get_Nyu() const { return Nyu; };
		CI get_Nzu() const { return Nzu; };

		CI get_Nxv() const { return Nxv; };
		CI get_Nyv() const { return Nyv; };
		CI get_Nzv() const { return Nzv; };

		CI get_Nxw() const { return Nxw; };
		CI get_Nyw() const { return Nyw; };
		CI get_Nzw() const { return Nzw; };

		CD get_Lx() const { return Lx; };
		CD get_Ly() const { return Ly; };
		CD get_Lz() const { return Lz; };

		CD get_dx() const { return dx; };
		CD get_dy() const { return dy; };
		CD get_dz() const { return dz; };

		map<int, Boundary> & get_BCs() { return BCs; }

		const VD & get_xp() const { return xp; }
		const VD & get_yp() const { return yp; }
		const VD & get_zp() const { return zp; }

		const VD & get_xu() const { return xu; }
		const VD & get_yu() const { return yu; }
		const VD & get_zu() const { return zu; }

		const VD & get_xv() const { return xv; }
		const VD & get_yv() const { return yv; }
		const VD & get_zv() const { return zv; }

		const VD & get_xw() const { return xw; }
		const VD & get_yw() const { return yw; }
		const VD & get_zw() const { return zw; }

	private:

		array<array<int, 3>, 4> N;

		int &Nx = N[0][0], &Ny = N[0][1], &Nz = N[0][2];
		int &Nxu = N[1][0], &Nyu = N[1][1], &Nzu = N[1][2];
		int &Nxv = N[2][0], &Nyv = N[2][1], &Nzv = N[2][2];
		int &Nxw = N[3][0], &Nyw = N[3][1], &Nzw = N[3][2];

		int NCells, Nyz;

		double Lx, Ly, Lz;
		double dx, dy, dz;

		map<int, Boundary> BCs;	
		
		vector<double> xp, yp, zp;
		vector<double> xu, yu, zu;
		vector<double> xv, yv, zv;
		vector<double> xw, yw, zw;
};
