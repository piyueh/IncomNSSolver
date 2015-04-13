/*
 * The Navier-Stoke solver using the basic 1st-order Euler time marching.
 */
class NSSolverEuler
{
	public:

		NSSolverEuler(Mesh &m, Fluid &f, double t, double Dt): 
			mesh(m), fluid(f), time(t), dt(Dt)
		{
			p.initShape(-1, mesh.Nx, -1, mesh.Nx, -1, mesh.Nz);
			u.initShape(-1, mesh.Nxu, -1, mesh.Nyu, -1, mesh.Nzu);
			v.initShape(-1, mesh.Nxv, -1, mesh.Nyv, -1, mesh.Nzv);
			w.initShape(-1, mesh.Nxw, -1, mesh.Nyw, -1, mesh.Nzw);

			u_str.initShape(-1, mesh.Nxu, -1, mesh.Nyu, -1, mesh.Nzu);
			v_str.initShape(-1, mesh.Nxv, -1, mesh.Nyv, -1, mesh.Nzv);
			w_str.initShape(-1, mesh.Nxw, -1, mesh.Nyw, -1, mesh.Nzw);

			p.setZeros();
			u.setZeros(); v.setZeros(); w.setZeros();

			pSolver.InitLinSys(mesh.Nx, mesh.Ny, mesh.Nz, 
					mesh.dx, mesh.dy, mesh.dz);
			pSolver.setLHS(mesh.BCs);
			pSolver.setRefP(0, 0, 0, 0);
			
		}

	private:

		double time;
		double dt;
		
		Mesh & mesh;

		Fluid & fluid;

		PoissonSolver pSolver;

		Array3D<double> u, v, w, p;
		Array3D<double> u_str, v_str, w_str;
};
