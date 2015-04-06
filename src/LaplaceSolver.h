SparseMatrix<double> InitializeA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz);
SparseMatrix<double> CreateA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz);
SparseMatrix<double> BCCorrectA(int Nx, int Ny, int Nz, 
		double dx, double dy, double dz, vector<int> & BCType);
