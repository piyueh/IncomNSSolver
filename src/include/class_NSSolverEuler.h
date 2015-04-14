/*
 * The Navier-Stoke solver using the basic 1st-order Euler time marching.
 */


# pragma once


class NSSolverEuler
{
	friend ostream &operator<<(ostream &, NSSolverEuler &);
	public:

		NSSolverEuler(Mesh &m, Fluid &f, double t, double Dt): 
			mesh(m), fluid(f), time(t), dt(Dt)
		{
			dx2 = dx * dx; dy2 = dy * dy; dz2 = dz * dz;

			p.initShape(-1, mesh.Nx, -1, mesh.Nx, -1, mesh.Nz);
			u.initShape(-1, mesh.Nxu, -1, mesh.Nyu, -1, mesh.Nzu);
			v.initShape(-1, mesh.Nxv, -1, mesh.Nyv, -1, mesh.Nzv);
			w.initShape(-1, mesh.Nxw, -1, mesh.Nyw, -1, mesh.Nzw);

			u_str.initShape(-1, mesh.Nxu, -1, mesh.Nyu, -1, mesh.Nzu);
			v_str.initShape(-1, mesh.Nxv, -1, mesh.Nyv, -1, mesh.Nzv);
			w_str.initShape(-1, mesh.Nxw, -1, mesh.Nyw, -1, mesh.Nzw);

			p.setZeros();
			u.setZeros(); v.setZeros(); w.setZeros();

			/*
			pSolver.InitLinSys(mesh.Nx, mesh.Ny, mesh.Nz, 
					mesh.dx, mesh.dy, mesh.dz);
			pSolver.setLHS(mesh.BCs);
			pSolver.setRefP(0, 0, 0, 0);
			*/
			
		}

		int updateGhost();
		int test();

		int output(string);
		int output_u();
		int output_v();
		int output_w();
		int output_p();

	private:

		double time, dt;

		Mesh & mesh;
		Fluid & fluid;

		int &Nx=mesh.Nx, &Ny=mesh.Ny, &Nz=mesh.Nz;
		int &Nxu=mesh.Nxu, &Nyu=mesh.Nyu, &Nzu=mesh.Nzu;
		int &Nxv=mesh.Nxv, &Nyv=mesh.Nyv, &Nzv=mesh.Nzv;
		int &Nxw=mesh.Nxw, &Nyw=mesh.Nyw, &Nzw=mesh.Nzw;

		double &dx=mesh.dx, &dy=mesh.dy, &dz=mesh.dz;

		double dx2, dy2, dz2;
		
		// PoissonSolver pSolver;

		Array3D<double> u, v, w, p;
		Array3D<double> u_str, v_str, w_str;

		int PredictStep();
		double ConvectU(int &, int &, int &);
		double ConvectV(int &, int &, int &);
		double ConvectW(int &, int &, int &);
		double DiffusiveU(int &, int &, int &);
		double DiffusiveV(int &, int &, int &);
		double DiffusiveW(int &, int &, int &);
};
