# include <iostream>
# include <vector>
# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>

using namespace std;
using namespace Eigen;


# include "LaplaceSolver.h"


int main()
{
	int Nx, Ny, Nz;
	int NCell;

	double Lx, Ly, Lz;
	double dx, dy, dz;

	Lx = 1.0; Ly = 1.0; Lz = 1.0;

	Nx = 3; Ny = 3; Nz = 1;
	NCell = Nx * Ny * Nz;

	dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz;

	auto A = CreateA(Nx, Ny, Nz, dx, dy, dz); 

	cout << A << endl;
	cout << "Row = " << A.rows() << endl;
	cout << "Col = " << A.cols() << endl;
	cout << "NonZeros = " << A.nonZeros() << endl;
	cout << A.isCompressed() << endl;
	cout << sizeof(A) << endl;



	return 0;
}


