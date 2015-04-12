# include "include/IncomNSSolver.h"
# include "testFuncs.cpp"
# include <cmath>


double evalRelErr(VectorXd & x, VectorXd & xe)
{
	auto tmp = x - xe;
	return tmp.cwiseAbs().maxCoeff();
}


int main()
{
	int Nx, Ny, Nz;
	int NCell;
	int Nyz;

	double Lx, Ly, Lz;
	double dx, dy, dz;

	ArrayXd x;
	ArrayXd y;
	ArrayXd z;

	VectorXd p;
	vector<Boundary> pBCs;
	PoissonSolver PressureSolver;

	VectorXd u;
	VectorXd v;
	VectorXd w;


	/* Initialize computational domain */
	{ 
		Lx = 1.0; Ly = 1.0; Lz = 1.0;
		Nx = 3; Ny = 3; Nz = 1;
		NCell = Nx * Ny * Nz;
		Nyz = Ny * Nz;
		dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz;

		x.setLinSpaced(Nx, dx/2, Lx-dx/2);
		y.setLinSpaced(Ny, dy/2, Ly-dy/2);
		z.setLinSpaced(Nz, dz/2, Lz-dz/2);
	}

	/* Initialize the primative unknowns */
	{
		p.resize(Nx*Ny);

		u.resize((Nx+1)*Ny*Nz);
		//u.setZero();
		u << 1, 2, 3, 4,
		  	 5, 6, 7, 8,
			 9, 10, 11, 12;

		v.resize(Nx*(Ny+1)*Nz);
		v.setZero();

		w.resize(Nx*Ny*(Nz+1));
		w.setZero();
	}

	/* Set up boundary conditions for the pressure and 
	 * initialize the pressure */
	{
		pBCs = genBCs(Nx, Ny, Nz);
		PressureSolver.InitLinSys(Nx, Ny, Nz, dx, dy, dz);
		PressureSolver.setLHS(pBCs);	
		PressureSolver.setRefP(0, 0, 0, 0);
	}

	VectorXd Fcu((Nx+1)*Ny);
	/* update u* */
	{
		int k = 0;
		int idx;
		double Fcu1, Fcu2, Fcu3;

		for(int i=1; i<Nx; ++i) {
			for(int j=1; j<Ny-1; ++j) {

				idx = i * Nyz + j * Nz + k;

				cout << i << ", " << j << ", " << idx << endl;

				cout << u[idx-Nyz] << ", " << u[idx] << ", " << u[idx+Nyz] << endl;

				Fcu1 = 0.25 * (u[idx-Nyz] + 2 * u[idx] + u[idx+Nyz]) * 
							  (u[idx+Nyz] - u[idx-Nyz]) / dx;


				Fcu2 = 0.25 * ((u[idx+Nz] + u[idx]) * (v[idx-Nyz+Nz] + v[idx+Nz]) - 
							   (u[idx] + u[idx-Nz]) * (v[idx-Nyz] + v[idx])) / dy;

				Fcu[idx] = Fcu1 + Fcu2;
				cout << Fcu[idx] << endl;
			}
		}
	}
	cout << Fcu << endl;


	/* Solve pressure field */
	{
		PressureSolver.setRHS(sourceTerm(Nx, Ny, Nz, x, y, z, 7)); 
		PressureSolver.Solve(p);
	}


	
	return 0;
}


