
class PoissonSolver
{
	public:
	
		PoissonSolver() = default;

		int InitLinearSys(int, int, int, double, double, double, vector<Boundary> &);

		void printA();
	
	private:

		int Nx, Ny, Nz;
		double dx, dy, dz;

		SparseMatrix<double> A;
		Matrix<double, Dynamic, Dynamic> b;
		Matrix<double, Dynamic, Dynamic> x;


		int InitA();
		int BCCorrectA(Boundary &);
};






























