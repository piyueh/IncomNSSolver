
double NSSolver::ConvectU(int &i, int &j, int &k)
{
	double Fcu;
	Fcu =
		0.25 * ((u(i+1, j, k) + u(i, j, k)) * (u(i+1, j, k) + u(i, j, k)) -
				(u(i, j, k) + u(i-1, j, k)) * (u(i, j, k) + u(i-1, j, k))) / dx + 

		0.25 * ((u(i, j+1, k) + u(i, j, k)) * (v(i, j+1, k) + v(i-1, j+1, k)) -
				(u(i, j, k) + u(i, j-1, k)) * (v(i, j, k) + v(i-1, j, k))) / dy +

		0.25 * ((u(i, j, k+1) + u(i, j, k)) * (w(i, j, k+1) + w(i-1, j, k+1)) -
				(u(i, j, k) + u(i, j, k-1)) * (w(i, j, k) + w(i-1, j, k))) / dz; 
	return Fcu;
}


double NSSolver::ConvectV(int &i, int &j, int &k)
{
	double Fcv;
	Fcv =
		0.25 * ((v(i+1, j, k) + v(i, j, k)) * (u(i+1, j, k) + u(i+1, j-1, k)) -
				(v(i, j, k) + v(i-1, j, k)) * (u(i, j, k) + u(i, j-1, k))) / dx + 

		0.25 * ((v(i, j+1, k) + v(i, j, k)) * (v(i, j+1, k) + v(i, j, k)) -
				(v(i, j, k) + v(i, j-1, k)) * (v(i, j, k) + v(i, j-1, k))) / dy +

		0.25 * ((v(i, j, k+1) + v(i, j, k)) * (w(i, j, k+1) + w(i, j-1, k+1)) -
				(v(i, j, k) + v(i, j, k-1)) * (w(i, j, k) + w(i, j-1, k))) / dz; 
	return Fcv;
}


double NSSolver::ConvectW(int &i, int &j, int &k)
{
	double Fcw;
	Fcw =
		0.25 * ((w(i+1, j, k) + v(i, j, k)) * (u(i+1, j, k) + u(i+1, j, k-1)) -
				(w(i, j, k) + v(i-1, j, k)) * (u(i, j, k) + u(i, j, k-1))) / dx +

		0.25 * ((w(i, j+1, k) + w(i, j, k)) * (v(i, j+1, k) + v(i, j+1, k-1)) -
				(w(i, j, k) + w(i, j-1, k)) * (v(i, j, k) + v(i, j, k-1))) / dy + 

		0.25 * ((w(i, j, k+1) + w(i, j, k)) * (w(i, j, k+1) + w(i, j, k)) -
				(w(i, j, k) + w(i, j, k-1)) * (w(i, j, k) + w(i, j, k-1))) / dz;

	return Fcw;
}

