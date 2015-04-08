
class PoissonSolver
{
	public:
	
		PoissonSolver() = default;

		int InitLHS(int, int, int, double, double, double, vector<Boundary> &);
		int InitRHS(Matrix<double, 1, Dynamic>);

		int setRefP(int [3], double);
		int setRefP(int , int, int, double);

		int Solve(VectorXd &);

		void printA();
		void printb();
		void printx();
	
	private:

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






























