/*
 * The Navier-Stoke solver using the basic 1st-order Euler time marching.
 */


# pragma once

class NSSolver
{
	friend ostream &operator<<(ostream &, NSSolver &);

	public:

		NSSolver(Mesh &m, Fluid &f, Data &d): mesh(m), fluid(f), data(d) {};

		int InitSolver(CD Dt, CI &, CI &, CaryI3 pIdx, CD pR);
		int InitSolver(string &);

		int solve();

	private:

		// computational parameters
		int targetNStep, outputN;
		double dt;
		double pRef;
		array<int, 3> pRefIdx;

		// fluid properties
		Fluid & fluid;
		double &nu = fluid.nu, &rho = fluid.rho;

		// mesh
		Mesh & mesh;
		map<int, Boundary> & BCs = mesh.BCs;

		int &Nx = mesh.Nx, &Ny = mesh.Ny, &Nz = mesh.Nz;
		int &Nxu = mesh.Nxu, &Nyu = mesh.Nyu, &Nzu = mesh.Nzu;
		int &Nxv = mesh.Nxv, &Nyv = mesh.Nyv, &Nzv = mesh.Nzv;
		int &Nxw = mesh.Nxw, &Nyw = mesh.Nyw, &Nzw = mesh.Nzw;

		int &Nyz = mesh.Nyz;

		double &dx = mesh.dx, &dy = mesh.dy, &dz = mesh.dz;

		vector<double> &xp = mesh.xp, &yp = mesh.yp, &zp = mesh.zp;
		vector<double> &xu = mesh.xu, &yu = mesh.yu, &zu = mesh.zu;
		vector<double> &xv = mesh.xv, &yv = mesh.yv, &zv = mesh.zv;
		vector<double> &xw = mesh.xw, &yw = mesh.yw, &zw = mesh.zw;


		// data of primative variables
		Data & data;
		double &time = data.time;
		Array3D<double> &u = data.u, &v = data.v, &w = data.w;
		VectorXd &p = data.p;


		// Poisson solver
		PoissonSolver pSolver;
		VectorXd b;


		// variables that are convinent for simplifying codes
		double dx2, dy2, dz2;
		int uBgIdx, uEdIdx, vBgIdx, vEdIdx, wBgIdx, wEdIdx; 
		Array3D<double> Gu, Gv, Gw;


		// private functions
		int InitLambda();
		int updateGhost();
		int PredictStep(CD & DT);
		int PredictStep(CD & DT, CD & coef);
		int updatePoissonSource(CD &);
		int updateU(CD &);

		function<double(CI &, CI &, CI &)> ConvU, ConvV, ConvW;
		function<double(CI &, CI &, CI &)> DiffU, DiffV, DiffW;
		function<void(CI &, CI &, CI &)> updGu, updGv, updGw;
		function<void(CI &, CI &, CI &, CD &)> updGu2, updGv2, updGw2;
		function<void(CI &, CI &, CI &, CD &)> preU, preV, preW;
		function<int(CI &, CI &, CI &)> calIdx;
		function<void(CI &, CI &)> updUB, updVB, updWB;
		function<void(CI &, CI &, CI &, CD &)> updU, updV, updW;
		function<void(CI &, CI &, CI &)> DivOnPresPt;

};
