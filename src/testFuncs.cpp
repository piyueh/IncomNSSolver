
vector<Boundary> genBCs(int Nx, int Ny, int Nz)
{
	vector<Boundary> BCs;
	{
		vector<int> temp; 	
		int i = 0;
		for(int j=0; j<Ny; ++j)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Ny*Nz, "-x", -1, 0.0, temp));
	}

	{
		vector<int> temp; 	
		int i = Nx - 1;
		for(int j=0; j<Ny; ++j)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Ny*Nz, "+x", -1, 0.0, temp));
	}

	{
		vector<int> temp; 	
		int j = 0;
		for(int i=0; i<Nx; ++i)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Nx*Nz, "-y", -1, 0.0, temp));
	}

	{
		vector<int> temp; 	
		int j = Ny - 1;
		for(int i=0; i<Nx; ++i)
		{
			for(int k=0; k<Nz; ++k)
				temp.push_back(i*Ny*Nz+j*Nz+k);
		}
		BCs.push_back(Boundary(Nx*Nz, "+y", -1, 0.0, temp));
	}

	{
		vector<int> temp, temp2; 	
		int k0 = 0, k1 = Nz - 1;
		for(int i=0; i<Nx; ++i)
		{
			for(int j=0; j<Ny; ++j)
			{
				temp.push_back(i*Ny*Nz+j*Nz+k0);
				temp2.push_back(i*Ny*Nz+j*Nz+k1);
			}
		}
		BCs.push_back(Boundary(Nx*Ny, "-z", 0, 0.0, temp, temp2));
	}


	{
		vector<int> temp, temp2; 	
		int k0 = Nz-1, k1 = 0;
		for(int i=0; i<Nx; ++i)
		{
			for(int j=0; j<Ny; ++j)
			{
				temp.push_back(i*Ny*Nz+j*Nz+k0);
				temp2.push_back(i*Ny*Nz+j*Nz+k1);
			}
		}
		BCs.push_back(Boundary(Nx*Ny, "+z", 0, 0.0, temp, temp2));
	}

	return BCs;
}


Array<double, 1, Dynamic> sourceTerm(int Nx, int Ny, int Nz,
		Array<double, 1, Dynamic> &X, 
		Array<double, 1, Dynamic> &Y,
	   	Array<double, 1, Dynamic> &Z, int n)
{
	int Nyz = Ny * Nz;
	double M_PI2 = M_PI * M_PI;

	auto xcos = (2 * M_PI * n * X).cos();
	auto ycos = (2 * M_PI * n * Y).cos();
	Array<double, 1, Dynamic> f(Nx*Ny*Nz);


	for(int i=0; i<Nx; ++i){
		for(int j=0; j<Ny; ++j){
			for(int k=0; k<Nz; ++k)
				f(i*Nyz+j*Nz+k) = - 8. * n * n * M_PI2 * xcos(i) * ycos(j); 		
		}
	}

	return f;
}


Array<double, 1, Dynamic> exactSoln(int Nx, int Ny, int Nz,
		Array<double, 1, Dynamic> &X, 
		Array<double, 1, Dynamic> &Y,
	   	Array<double, 1, Dynamic> &Z, int n)
{
	int Nyz = Ny * Nz;
	double M_PI2 = M_PI * M_PI;

	auto xcos = (2 * M_PI * n * X).cos();
	auto ycos = (2 * M_PI * n * Y).cos();
	Array<double, 1, Dynamic> P(Nx*Ny*Nz);


	for(int i=0; i<Nx; ++i){
		for(int j=0; j<Ny; ++j){
			for(int k=0; k<Nz; ++k)
				P(i*Nyz+j*Nz+k) = xcos(i) * ycos(j); 		
		}
	}

	return P;
}
