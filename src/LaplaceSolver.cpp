
# include <vector>
# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>


using namespace std;
using namespace Eigen;


SparseMatrix<double> CreateA(int *N, double *dL)
{
	SparseMatrix<double> A;

	A = InitializeA(N, dL);

	return A;
}


/*
 * Initialize the matrix A in the linear system Ax=b for the 3D Poisson problems.
 * The boundary conditions will not be implemented in the returned matrix A.
 */
SparseMatrix<double> InitializeA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz)
{
	// the numbers of diagonal sub-diagonal elements
	int NDG0 = Nx * Ny * Nz;
	int NDG1 = NDG0 - 1;
	int NDG2 = NDG0 - Nz;
	int NDG3 = NDG0 - Ny * Nz;

	double dxInv2, dyInv2, dzInv2;
	double dgValue;

	// the container used to initialize the matrix A
	vector<Triplet<double>> CoeffMatrix;

	// the matrix A in matrix equation Ax=b for the Laplace equation
	SparseMatrix<double> A;

	// values of non-diagonal elements
	dxInv2 = 1.0 / (dx * dx);
	dyInv2 = 1.0 / (dy * dy);
	dzInv2 = 1.0 / (dy * dy);

	// the value of diagonal elements
	dgValue = - 2.0 * (dxInv2 + dyInv2 + dzInv2);

	// assign the values to the 0th diagonal elements
	for(int idx=0; idx<NDG0; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx, dgValue));

	// assign the values to the 1st diagonal elements
	for(int idx=0; idx<NDG1; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx+1, dzInv2));

	// assign the values to the -1st diagonal elements
	for(int idx=1; idx<NDG0; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx-1, dzInv2));

	// assign the values to the 2nd diagonal elements (A[i, i+Nz])
	for(int idx=0; idx<NDG2; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx+Nz, dyInv2));

	// assign the values to the -2nd diagonal elements (A[i, i-Nz])
	for(int idx=Nz; idx<NDG0; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx-Nz, dyInv2));

	// assign the values to the 3nd diagonal elements (A[i, i+Ny*Nz])
	for(int idx=0; idx<NDG3; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx+Ny*Nz, dxInv2));

	// assign the values to the -3nd diagonal elements (A[i, i-Ny*Nz])
	for(int idx=Ny*Nz; idx<NDG0; ++idx)
		CoeffMatrix.push_back(Triplet<double>(idx, idx-Ny*Nz, dxInv2));

	// allocate the space for the sparse matrix A
	A.resize(NDG0, NDG0);

	// initialize the sparse matrix A
	A.setFromTriplets(CoeffMatrix.begin(), CoeffMatrix.end());

	return A;
}


/*
 * Create the correction matrix for the matrix A according to the boundary condidtion.
 * The sparse matrix form is used.
 */
SparseMatrix<double> BCCorrectA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz, string dir, int BCType, double BCValue)
{
	SparseMatrix<double> BC;

	return BC;
}
