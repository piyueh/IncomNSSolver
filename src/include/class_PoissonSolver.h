/*
 * Class of the Poisson solver.
 */


# pragma once


class PoissonSolver
{
	public:
	
		PoissonSolver() = default;
		PoissonSolver(const array<int, 3> n, const array<double, 3> d) 
		{InitLinSys(n, d);}

		int InitLinSys(const array<int, 3>, const array<double, 3>);

		int setRefP(const array<int, 3> &, const double &);

		int setLHS(map<int, Boundary> &);

# ifdef BICGSTAB
		int setTolerance(CD &);
# endif

		pair<int, double> Solve(Map<VectorXd> &, Map<VectorXd> &);

		void printA();
	
	private:

		int Nx, Ny, Nz;
		int NCells, Nyz;

		double dx, dy, dz;

		int refLoc;
		double refP=0;


		SparseMatrix<double> A;


# ifdef BICGSTAB
		BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> Solver;
# endif

# ifdef PARDISO_LDLT
		PardisoLDLT<SparseMatrix<double>, Upper> Solver;
# endif

		int InitA();
		int getIdx(const int & i, const int & j, const int & k); 
		int getIdx(const array<int, 3> & idx);
		pair<vector<int>, vector<int>> ObtainCells(const unsigned int &, const int &);
		int BCCorrectA(const unsigned int &, const int &, const int &, const double &);
};


