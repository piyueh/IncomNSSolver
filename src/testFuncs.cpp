VectorXd sourceTerm(int Nx, int Ny, int Nz,
		ArrayXd &X, ArrayXd &Y, ArrayXd &Z, int n)
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

	return f.matrix();
}


VectorXd exactSoln(int Nx, int Ny, int Nz,
		ArrayXd &X, ArrayXd &Y, ArrayXd &Z, int n)
{
	int Nyz = Ny * Nz;

	auto xcos = (2 * M_PI * n * X).cos();
	auto ycos = (2 * M_PI * n * Y).cos();
	Array<double, 1, Dynamic> P(Nx*Ny*Nz);


	for(int i=0; i<Nx; ++i){
		for(int j=0; j<Ny; ++j){
			for(int k=0; k<Nz; ++k)
				P(i*Nyz+j*Nz+k) = xcos(i) * ycos(j); 		
		}
	}

	return P.matrix();
}
