# include "include/IncomNSSolver.h"

int main()
{
	int Nx, Ny, Nz;
	int NCell;

	double Lx, Ly, Lz;
	double dx, dy, dz;

	vector<Boundary> BCs;
	PoissonSolver test;



	Lx = 1.0; Ly = 1.0; Lz = 1.0;

	Nx = 3; Ny = 3; Nz = 1;
	NCell = Nx * Ny * Nz;

	dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz;



	{
		vector<int> temp; 	
		int i = 0;
		for(int j=0; j<Ny; ++j)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Ny*Nz, "-x", temp, 1, 0.0));
	}

	{
		vector<int> temp; 	
		int i = Nx - 1;
		for(int j=0; j<Ny; ++j)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Ny*Nz, "+x", temp, 1, 0.0));
	}

	{
		vector<int> temp; 	
		int j = 0;
		for(int i=0; i<Nx; ++i)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Nx*Nz, "-y", temp, 1, 0.0));
	}

	{
		vector<int> temp; 	
		int j = Ny - 1;
		for(int i=0; i<Nx; ++i)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Nx*Nz, "+y", temp, 1, 0.0));
	}

	/*
	{
		vector<int> temp; 	
		int k = 0;
		for(int i=0; i<Nx; ++i)
		{
			for(int j=0; j<Ny; ++j)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Nx*Ny, "-z", temp, 1, 0.0));
	}


	{
		vector<int> temp; 	
		int k = Nz - 1;
		for(int i=0; i<Nx; ++i)
		{
			for(int j=0; j<Ny; ++j)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Nx*Ny, "+z", temp, 1, 0.0));
	}
	*/


	test.InitLinearSys(Nx, Ny, Nz, dx, dy, dz, BCs);	
	test.printA();


	return 0;
}


