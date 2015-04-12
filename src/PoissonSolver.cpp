# include "include/IncomNSSolver.h"

PoissonSolver::PoissonSolver(int N1, int N2, int N3, double d1, double d2, double d3)
{
	Nx = N1; Ny = N2; Nz = N3;
	dx = d1; dy = d2; dz = d3;

	N = Nx * Ny * Nz;

	A.resize(N, N);
	A.setZero();

	b.resize(N);
	b.setZero();
}	


int PoissonSolver::InitLinSys(int N1, int N2, int N3, double d1, double d2, double d3)
{
	Nx = N1; Ny = N2; Nz = N3;
	dx = d1; dy = d2; dz = d3;

	N = Nx * Ny * Nz;

	A.resize(N, N);
	A.setZero();

	b.resize(N);
	b.setZero();

	return 0;
}	



int PoissonSolver::setLHS(vector<Boundary> & BC)
{
	int err;

	err = InitA();

	for(auto &it: BC)
		err = BCCorrectA(it);

	cgSolver.compute(A);

	return 0;
}


int PoissonSolver::setRHS(VectorXd f)
{
	assert(f.size()==Nx*Ny*Nz);

	b += f;
	return 0;
}


int PoissonSolver::setRefP(int idx[3], double value)
{
	int I = idx[0] * Ny * Nz + idx[1] * Nz + idx[2];

	refIdx[0] = idx[0];
	refIdx[1] = idx[1];
	refIdx[2] = idx[2];
	refP = value;

	b.coeffRef(I) -= value * A.coeffRef(I, I);
	A.coeffRef(I, I) = 0.;
	return 0;
}


int PoissonSolver::setRefP(int i, int j, int k, double value)
{
	int I = i * Ny * Nz + j * Nz + k;

	refIdx[0] = i;
	refIdx[1] = j;
	refIdx[2] = k;
	refP = value;

	b.coeffRef(I) -= value * A.coeffRef(I, I);
	A.coeffRef(I, I) = 0.;
	return 0;
}


int PoissonSolver::Solve(VectorXd &soln)
{
	soln = cgSolver.solve(b.matrix());

	return 0;
}

/*
 * Initialize the matrix A in the linear system Ax=b for the 3D Poisson problems.
 * The boundary conditions will not be implemented in the returned matrix A.
 */
int PoissonSolver::InitA()
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


	// values of non-diagonal elements
	dxInv2 = 1.0 / (dx * dx);
	dyInv2 = 1.0 / (dy * dy);
	dzInv2 = 1.0 / (dz * dz);

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
	A.setFromTriplets(CoeffMatrix.cbegin(), CoeffMatrix.cend());
	A.makeCompressed();


	return 0;
}


/*
 * Create the correction matrix for the matrix A according to the boundary condidtion.
 * The sparse matrix form is used.
 */
int PoissonSolver::BCCorrectA(Boundary & Surf)
{
	int j0;
	int type = Surf.get_pType();

	double v = Surf.get_pBCvalue();
	double adjValue;

	auto bg = Surf.bgCell();
	auto ed = Surf.edCell();


	switch (Surf.get_Direction()[1])
	{
		case 'x': j0 = Ny * Nz; adjValue = 1 / (dx * dx); break;
		case 'y': j0 = Nz; adjValue = 1 / (dy * dy); break;
		case 'z': j0 = 1; adjValue = 1 / (dz * dz); break;
		default:
			throw invalid_argument("The direction of the boundary is wrong!");
			break;
	}

	switch (Surf.get_Direction()[0])
	{
		case '+': break;
		case '-': j0 *= (-1); break;
		default:
			throw invalid_argument("The sign of the boundary is wrong!");
			break;
	}
		

	switch (type)
	{
		case 0:

			for(auto it=bg, it2=Surf.bgOppCell(); it<ed; ++it, ++it2){
				A.coeffRef(*it, *it2) += adjValue;
				if (*it + j0 >= 0 && *it + j0 < A.cols())
					A.coeffRef(*it, *it+j0) -= adjValue;
			}
			break;

		case 1:

			for(auto it=bg; it<ed; ++it){
				A.coeffRef(*it, *it) -= adjValue;
				if (*it + j0 >= 0 && *it + j0 < A.cols())
					A.coeffRef(*it, *it+j0) -= adjValue;
			}
			break;

		case -1:

			for(auto it=bg; it<ed; ++it){
				A.coeffRef(*it, *it) += adjValue;
				if (*it + j0 >= 0 && *it + j0 < A.cols())
					A.coeffRef(*it, *it+j0) -= adjValue;
			}
			break;

		default:

			throw invalid_argument("The sign of the boundary is wrong!");
			break;
	}


	return 0;
}


void PoissonSolver::printA()
{
	cout << "LHS, Matrix A: " << endl;
	cout << A << endl;
}


void PoissonSolver::printb()
{
	cout << "RHS, Vector b: " << endl;
	cout << b << endl;
}


