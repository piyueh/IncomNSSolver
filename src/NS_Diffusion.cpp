//double NSSolver

double NSSolver::DiffusiveU(int &i, int &j, int &k)
{
	double Fdu;
	Fdu = fluid.nu * (
			(u(i+1, j, k) - 2 * u(i, j, k) + u(i-1, j, k)) / dx2 +
			(u(i, j+1, k) - 2 * u(i, j, k) + u(i, j-1, k)) / dy2 +
			(u(i, j, k+1) - 2 * u(i, j, k) + u(i, j, k-1)) / dz2 );
	return Fdu;
}


double NSSolver::DiffusiveV(int &i, int &j, int &k)
{
	double Fdv;
	Fdv = fluid.nu * (
			(v(i+1, j, k) - 2 * v(i, j, k) + v(i-1, j, k)) / dx2 +
			(v(i, j+1, k) - 2 * v(i, j, k) + v(i, j-1, k)) / dy2 +
			(v(i, j, k+1) - 2 * v(i, j, k) + v(i, j, k-1)) / dz2 );
	return Fdv;
}


double NSSolver::DiffusiveW(int &i, int &j, int &k)
{
	double Fdw;
	Fdw = fluid.nu * (
			(w(i+1, j, k) - 2 * w(i, j, k) + w(i-1, j, k)) / dx2 +
			(w(i, j+1, k) - 2 * w(i, j, k) + w(i, j-1, k)) / dy2 +
			(w(i, j, k+1) - 2 * w(i, j, k) + w(i, j, k-1)) / dz2 );
	return Fdw;
}
