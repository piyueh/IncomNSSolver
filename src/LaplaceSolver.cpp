
# include <vector>
# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>


using namespace std;
using namespace Eigen;


SparseMatrix<double> CreateA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz)
{
	SparseMatrix<double> A;

	A = InitializeA(Nx, Ny, Nz, dx, dy, dz);

	return A;
}


SparseMatrix<double> InitializeA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz)
{
	int NCell = Nx * Ny * Nz;

	double dxInv2, dyInv2, dzInv2;
	double dgValue;

	vector<Triplet<double>> CoeffMatrix;
	SparseMatrix<double> A;

	dxInv2 = 1.0 / (dx * dx);
	dyInv2 = 1.0 / (dy * dy);
	dzInv2 = 1.0 / (dz * dz);

	dgValue = - 2.0 * (dxInv2 + dyInv2 + dzInv2);

	for(int idx=0; idx<NCell; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx, dgValue));

	for(int idx=0; idx<NCell-1; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx+1, dxInv2));

	for(int idx=1; idx<NCell; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx-1, dxInv2));

	for(int idx=0; idx<NCell-Nx; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx+Nx, dyInv2));

	for(int idx=Nx; idx<NCell; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx-Nx, dyInv2));

	for(int idx=0; idx<NCell-Nx*Ny; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx+Nx*Ny, dzInv2));

	for(int idx=Nx*Ny; idx<NCell; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx-Nx*Ny, dzInv2));

	A.resize(NCell, NCell);
	A.setFromTriplets(CoeffMatrix.begin(), CoeffMatrix.end());

	return A;
}


SparseMatrix<double> BCCorrectA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz, vector<int> & BCType)
{
	SparseMatrix<double> BC;

	return BC;
}
