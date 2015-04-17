/*
 * The Navier-Stoke solver using the basic 1st-order Euler time marching.
 */


# pragma once
# include <cmath>

class NSSolverEuler
{
	friend ostream &operator<<(ostream &, NSSolverEuler &);

	public:

		NSSolverEuler(Mesh &m, Fluid &f): mesh(m), fluid(f) {};

		int InitSolver(double t, double Dt, 
				array<int, 3> pIdx, double pR)
		{
			dt = Dt; time = t;
			pRefIdx = pIdx; pRef = pR;

			dx2 = dx * dx; dy2 = dy * dy; dz2 = dz * dz;

			u.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
			v.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
			w.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

			u_str.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
			v_str.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
			w_str.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

			u.setZeros(); v.setZeros(); w.setZeros();

			for(int i=0; i<Nxu; ++i){
				for(int j=0; j<Nyu; ++j){
					for(int k=0; k<Nzu; ++k){
						u(i, j, k) = - cos(xu[i]) * sin(yu[j]);
					}
				}
			}

			for(int i=0; i<Nxv; ++i){
				for(int j=0; j<Nyv; ++j){
					for(int k=0; k<Nzv; ++k){
						v(i, j, k) = sin(xv[i]) * cos(yv[j]);
					}
				}
			}

			b.resize(Nx * Ny * Nz);
			p.resize(Nx * Ny * Nz);

			pSolver.InitLinSys({Nx, Ny, Nz}, {dx, dy, dz});
			pSolver.setLHS(mesh.get_BCs());
			pSolver.setRefP(pRefIdx, pRef);
			
			return 0;
		}

		int solve(int);

		int output(string);
		int output_u(); int output_v(); int output_w(); int output_p();

	private:

		double time, dt;

		Mesh & mesh;
		Fluid & fluid;

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
