/*
 * Class of the Poisson solver.
 */


# pragma once


class PoissonSolver
{
	public:
	
		PoissonSolver() = default;
		PoissonSolver(const array<int, 3> n, const array<double, 3> d) {InitLinSys(n, d);}

		int InitLinSys(const array<int, 3>, const array<double, 3>);

		int setRefP(const array<int, 3> &, const double &);

		int setLHS(const map<int, Boundary> &);

		pair<int, double> Solve(VectorXd &, VectorXd &);

		void printA();
	
	private:

		int Nx, Ny, Nz;
		int NCells, Nyz;

		double dx, dy, dz;

		int refLoc;
		double refP=0;


		SparseMatrix<double> A;
		//ConjugateGradient<SparseMatrix<double>> cgSolver;
		BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> cgSolver;


		int InitA();
		int getIdx(const int & i, const int & j, const int & k); 
		int getIdx(const array<int, 3> & idx);
		pair<vector<int>, vector<int>> ObtainCells(const unsigned int &, const int &);
		int BCCorrectA(const unsigned int &, const int &, const int &, const double &);
};


