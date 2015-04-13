# include "include/IncomNSSolver.h"


int main()
{

	Mesh testMesh;
	Fluid testFluid(1., 1.);

	testMesh.InitMesh({4, 4, 1}, {1, 1, 1});

	testMesh.addBC(1, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	testMesh.addBC(2, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	testMesh.addBC(3, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	testMesh.addBC(1, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	testMesh.addBC(2, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	testMesh.addBC(3, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	

	NSSolverEuler testSolve(testMesh, testFluid, 0., 0.1);

	testSolve.test();

	testSolve.output("A.txt");

	testSolve.updateGhost();

	testSolve.output("B.txt");

	return 0;
}
