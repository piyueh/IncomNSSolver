/*
 * The Navier-Stoke solver using the basic 1st-order Euler time marching.
 */


# pragma once

class NSSolver
{
	friend ostream &operator<<(ostream &, NSSolver &);

	public:

		NSSolver(Mesh &m, Fluid &f): mesh(m), fluid(f) {};

		int InitSolver(CD t, CD Dt, CaryI3 pIdx, CD pR);

		int solve(int);

		int output(string);
		int output_u(); int output_v(); int output_w(); int output_p();

	private:

		double time, dt;

		Mesh & mesh;
		Fluid & fluid;

		map<int, Boundary> & BCs = mesh.BCs;

		int &Nx=mesh.Nx, &Ny=mesh.Ny, &Nz=mesh.Nz;
		int &Nxu=mesh.Nxu, &Nyu=mesh.Nyu, &Nzu=mesh.Nzu;
		int &Nxv=mesh.Nxv, &Nyv=mesh.Nyv, &Nzv=mesh.Nzv;
		int &Nxw=mesh.Nxw, &Nyw=mesh.Nyw, &Nzw=mesh.Nzw;

		int &Nyz = mesh.Nyz;

		double &dx=mesh.dx, &dy=mesh.dy, &dz=mesh.dz;

		vector<double> &xp=mesh.xp, &yp=mesh.yp, &zp=mesh.zp;
		vector<double> &xu=mesh.xu, &yu=mesh.yu, &zu=mesh.zu;
		vector<double> &xv=mesh.xv, &yv=mesh.yv, &zv=mesh.zv;
		vector<double> &xw=mesh.xw, &yw=mesh.yw, &zw=mesh.zw;

		double dx2, dy2, dz2;
		int uBgIdx, uEdIdx, vBgIdx, vEdIdx, wBgIdx, wEdIdx; 
		
		PoissonSolver pSolver;
		array<int, 3> pRefIdx;
		double pRef;


		Array3D<double> u, v, w;
		Array3D<double> u_str, v_str, w_str;

		VectorXd p, b;

		int updateGhost();

		int PredictStep();
		double ConvectU(int &, int &, int &);
		double ConvectV(int &, int &, int &);
		double ConvectW(int &, int &, int &);
		double DiffusiveU(int &, int &, int &);
		double DiffusiveV(int &, int &, int &);
		double DiffusiveW(int &, int &, int &);

		int updatePoissonSource();

		int updateU();
};
