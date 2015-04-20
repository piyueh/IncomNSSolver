# include "include/IncomNSSolver.h"


int PoissonSolver::InitLinSys(const array<int, 3> n, const array<double, 3> dl)
{
	Nx = n[0]; Ny = n[1]; Nz = n[2]; 
	dx = dl[0]; dy = dl[1]; dz = dl[2]; 

	NCells = Nx * Ny * Nz; Nyz = Ny * Nz;

	A.resize(NCells, NCells);
	A.setZero();

	return 0;
}	



int PoissonSolver::setLHS(const map<int, Boundary> & BC)
{

	InitA();

	for(auto &it: BC)
		BCCorrectA(it.second.get_Dir(), it.second.get_Sign(), 
				it.second.get_pType(), it.second.get_pBCvalue());

	cgSolver.compute(A);

	return 0;
}


int PoissonSolver::setRefP(const array<int, 3> & Idx, const double & value)
{
	refLoc = getIdx(Idx);
	refP = value;

	A.coeffRef(refLoc, refLoc) = 1;
	
	if ((refLoc - Nyz) >= 0) A.coeffRef(refLoc, refLoc - Nyz) = 0;
	if ((refLoc + Nyz) < NCells) A.coeffRef(refLoc, refLoc + Nyz) = 0;
	if ((refLoc - Nz) >= 0) A.coeffRef(refLoc, refLoc - Nz) = 0;
	if ((refLoc + Nz) < NCells) A.coeffRef(refLoc, refLoc + Nz) = 0;
	if ((refLoc - 1) >= 0) A.coeffRef(refLoc, refLoc - 1) = 0;
	if ((refLoc + 1) < NCells) A.coeffRef(refLoc, refLoc + 1) = 0;

	A.makeCompressed();

	return 0;
}


pair<int, double> PoissonSolver::Solve(VectorXd & f, VectorXd & soln)
{
	assert(f.size() == NCells);

	f[refLoc] = refP;

	soln = cgSolver.solveWithGuess(f, soln);

	return {cgSolver.iterations(), cgSolver.error()};
}


/*
 * Initialize the matrix A in the linear system Ax=b for the 3D Poisson problems.
 * The boundary conditions will not be implemented in the returned matrix A.
 */
int PoissonSolver::InitA()
{
	// the numbers of diagonal sub-diagonal elements
	int & NDG0 = NCells;
	int NDG1 = NDG0 - 1;
	int NDG2 = NDG0 - Nz;
	int NDG3 = NDG0 - Nyz;

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
int PoissonSolver::BCCorrectA(const unsigned int & dir, const int & sign,
		const int & type, const double & value)
{
	int jTarget;
	double adjValue;

	const auto Cells = ObtainCells(dir, sign);
	auto bg_i = Cells.first.cbegin(), ed_i = Cells.first.cend();
	decltype(bg_i) bg_j;

	switch (dir)
	{
		case 1: jTarget = Nyz; adjValue = 1 / (dx * dx); break;
		case 2: jTarget = Nz; adjValue = 1 / (dy * dy); break;
		case 3: jTarget = 1; adjValue = 1 / (dz * dz); break;
		default:
			throw invalid_argument("The direction of the boundary is wrong!");
			break;
	}

	jTarget *= sign;

	if (type == 1) adjValue *= -1;

	switch (type)
	{
		case 0: bg_j = Cells.second.cbegin(); break;
		case 1: bg_j = Cells.first.cbegin(); adjValue *= -1; break;
		case -1: bg_j = Cells.first.cbegin(); break;
		default:
			throw invalid_argument("The sign of the boundary is wrong!");
			break;
	}

	for(auto it=bg_i, it2=bg_j; it<ed_i; ++it, ++it2){
		A.coeffRef(*it, *it2) += adjValue;
		if (*it + jTarget >= 0 && *it + jTarget < A.cols())
			A.coeffRef(*it, *it+jTarget) -= adjValue;
	}

	return 0;
}


inline pair<vector<int>, vector<int>>  
PoissonSolver::ObtainCells(const unsigned int & dir, const int & sign)
{
	int ia = 0, ib = Nx;
	int ja = 0, jb = Ny;
	int ka = 0, kb = Nz;
	int adjIdx, idx;

	pair<vector<int>, vector<int>> Cells;

	switch (dir)
	{
		case 1: ib = 1; adjIdx = (Nx - 1) * Nyz; break;
		case 2: jb = 1; adjIdx = (Ny - 1) * Nz; break;
		case 3: kb = 1; adjIdx = Nz - 1; break;
	}

	for(int i=ia; i<ib; ++i) {
		for(int j=ja; j<jb; ++j) {
			for(int k=ka; k<kb; ++k) {
				idx = getIdx(i, j, k);
				Cells.first.push_back(idx);
				Cells.second.push_back(idx + adjIdx);
			}
		}
	}

	if (sign == 1) Cells.first.swap(Cells.second);

	return Cells;
}


inline int PoissonSolver::getIdx(const int & i, const int & j, const int & k) 
{ return i * Nyz + j * Nz + k; }

inline int PoissonSolver::getIdx(const array<int, 3> & idx) 
{ return idx[0] * Nyz + idx[1] * Nz + idx[2]; }


void PoissonSolver::printA()
{
	cout << "LHS, Matrix A: " << endl;
	cout << A << endl;
}



