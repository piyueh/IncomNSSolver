/*
 * Class of the Poisson solver.
 */


# pragma once


class PoissonSolver
{
	public:
	
		PoissonSolver() = default;
		PoissonSolver(int a, int b, int c, double d, double e, double f)
			{InitLinSys(a, b, c, d, e, f);}

		int InitLinSys(int, int, int, double, double, double);

		int setLHS(vector<Boundary> &);
		int setRHS(VectorXd);

		int setRefP(int [3], double);
		int setRefP(int, int, int, double);

		int Solve(VectorXd &);

		void printA();
		void printb();
		void printx();
	
	private:

		int N;
		int Nx, Ny, Nz;
		double dx, dy, dz;

		int refIdx[3]={0, 0, 0};
		double refP=0;

		SparseMatrix<double> A;
		VectorXd b;

		ConjugateGradient<SparseMatrix<double>> cgSolver;


		int InitA();
		int BCCorrectA(Boundary &);
};


