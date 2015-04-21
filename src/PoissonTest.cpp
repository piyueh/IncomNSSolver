# include "include/IncomNSSolver.h"
# include "testFuncs.cpp"


double evalRelErr(VectorXd & x, VectorXd & xe)
{
	auto tmp = x - xe;
	//return tmp.array().pow(2).sum() / xe.array().pow(2).sum();
	return tmp.cwiseAbs().maxCoeff();
}


int main()
{
	int Nx, Ny, Nz;
	double Lx, Ly, Lz;
	double dx, dy, dz;

	ArrayXd x, y, z;

	VectorXd p, f;
	VectorXd p_exact;

	Mesh mesh;
	PoissonSolver solver;


	Lx = 1.0; Ly = 1.0; Lz = 1.0;
	Nx = 100; Ny = 100; Nz = 1;
	dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz; 


	x.setLinSpaced(Nx, dx/2, Lx-dx/2);
	y.setLinSpaced(Ny, dy/2, Ly-dy/2);
	z.setLinSpaced(Nz, dz/2, Lz-dz/2);

	p_exact = exactSoln(Nx, Ny, Nz, x, y, z, 25);

	f = sourceTerm(Nx, Ny, Nz, x, y, z, 25);

	mesh.InitMesh({Nx, Ny, Nz}, {Lx, Ly, Lz});

	mesh.addBC(1, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(2, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(3, 1, {0, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(1, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(2, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(3, -1, {0, 0}, {0, 0}, {0, 0}, {0, 0});


	
	solver.InitLinSys({Nx, Ny, Nz}, {dx, dy, dz});
	solver.setLHS(mesh.get_BCs());
	solver.setRefP({0, 0, 0}, p_exact[0]);

	solver.Solve(f, p);

	double err = evalRelErr(p, p_exact);

	cout << err << endl;

	
	return 0;
}


