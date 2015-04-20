# include "include/IncomNSSolver.h"


int main()
{
	int Nx, Ny, Nz;
	double Lx, Ly, Lz;

	Fluid fluid(0.01, 1.);
	//Fluid fluid(1., 1.);
	Mesh mesh;
	Data data;
	NSSolver solver(mesh, fluid, data);

	Nx = 100; Ny = 100; Nz = 1;
	Lx = 1; Ly = 1; Lz = 0.01;
	//Lx = 2*M_PI; Ly = 2*M_PI; Lz = 0.1;

	mesh.InitMesh({Nx, Ny, Nz}, {Lx, Ly, Lz});

	/*
	mesh.addBC(1, 1, {-1, 0}, {-1, 0}, {1, 0}, {0, 0});
	mesh.addBC(2, 1, {-1, 0}, {1, 0}, {-1, 0}, {0, 0});
	mesh.addBC(3, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(1, -1, {-1, 0}, {-1, 0}, {1, 0}, {0, 0});
	mesh.addBC(2, -1, {-1, 0}, {1, 0}, {-1, 0}, {0, 0});
	mesh.addBC(3, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	*/
	mesh.addBC(1, 1, {-1, 0}, {1, 0}, {1, 0}, {0, 0});
	mesh.addBC(2, 1, {-1, 0}, {1, 1}, {1, 0}, {0, 0});
	mesh.addBC(3, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	mesh.addBC(1, -1, {-1, 0}, {1, 0}, {1, 0}, {0, 0});
	mesh.addBC(2, -1, {-1, 0}, {1, 0}, {1, 0}, {0, 0});
	mesh.addBC(3, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});

	data.InitData(mesh);	

	solver.InitSolver(0.001, {Nx-1, 0, 0}, 0.);

	solver.solve(150000, 50);
	//solver.solve(1);
	
	data.output("Data.txt");

	return 0;
}


