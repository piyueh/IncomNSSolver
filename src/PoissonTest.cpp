# include "include/IncomNSSolver.h"
# include "testFuncs.cpp"
# include <cmath>


double evalRelErr(VectorXd & x, VectorXd & xe)
{
	auto tmp = x - xe;
	//return tmp.array().pow(2).sum() / xe.array().pow(2).sum();
	return tmp.cwiseAbs().maxCoeff();
}


int main()
{
	int Nx, Ny, Nz;
	int NCell;

	double Lx, Ly, Lz;
	double dx, dy, dz;

	ArrayXd x;
	ArrayXd y;
	ArrayXd z;

	VectorXd p;
	VectorXd p_exact;

	vector<Boundary> BCs;
	PoissonSolver test;


	Lx = 1.0; Ly = 1.0; Lz = 1.0;

	Nx = 100; Ny = 100; Nz = 1;
	NCell = Nx * Ny * Nz;

	dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz;

	x.setLinSpaced(Nx, dx/2, Lx-dx/2);
	y.setLinSpaced(Ny, dy/2, Ly-dy/2);
	z.setLinSpaced(Nz, dz/2, Lz-dz/2);

	p_exact = exactSoln(Nx, Ny, Nz, x, y, z, 7);

	BCs = genBCs(Nx, Ny, Nz);

	test.InitLinSys(Nx, Ny, Nz, dx, dy, dz);

	test.setLHS(BCs);	

	test.setRefP(0, 0, 0, p_exact[0]);
	test.setRHS(sourceTerm(Nx, Ny, Nz, x, y, z, 7)); 


	test.Solve(p);

	cout << p << endl;

	double err = evalRelErr(p, p_exact);

	cout << err << endl;
	cout << sqrt(err) << endl;

	
	return 0;
}


