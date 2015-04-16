# include "include/IncomNSSolver.h"
# include "testFuncs.cpp"


double evalRelErr(VectorXd & x, VectorXd & xe)
{
	auto tmp = x - xe;
	return tmp.cwiseAbs().maxCoeff();
}


int main()
{
	int Nx, Ny, Nz;
	double Lx, Ly, Lz;

	Fluid fluid(1, 1);
	Mesh mesh;
	NSSolverEuler solver(mesh, fluid);

	Nx = 100; Ny = 100; Nz = 1;
	Lx = M_PI * 2; Ly = M_PI * 2; Lz = 1;

	mesh.InitMesh({Nx, Ny, Nz}, {Lx, Ly, Lz});

	mesh.addBC(1, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(2, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(3, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(1, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(2, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(3, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});

	solver.InitSolver(0., 0.001, {Nx-1, 0, 0}, 0.);
	solver.solve(2000);
	
	solver.output("Data.txt");

	return 0;
}


