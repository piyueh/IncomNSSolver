# include "include/IncomNSSolver.h"


int updXGhost(Array3D<double> &, int &, int &, 
		int &, int &, int &, double, double &);
int updYGhost(Array3D<double> &, int &, int &, 
		int &, int &, int &, double, double &);
int updZGhost(Array3D<double> &, int &, int &, 
		int &, int &, int &, double, double &);


int NSSolverEuler::updateGhost()
{
	for(auto &bcPair: mesh.BCs)
	{
		auto & bc = bcPair.second;
		
		switch (bc.dir) {
			case 1: // x

				updXGhost(u, mesh.Nyu, mesh.Nzu, 
						bc.uType, bc.uBCIdx, bc.uBCcorIdx, 2*mesh.dx, bc.uBCvalue);

				updXGhost(v, mesh.Nyv, mesh.Nzv, 
						bc.vType, bc.vBCIdx, bc.vBCcorIdx, mesh.dx, bc.vBCvalue);

				updXGhost(w, mesh.Nyw, mesh.Nzw, 
						bc.wType, bc.wBCIdx, bc.wBCcorIdx, mesh.dx, bc.wBCvalue);

				updXGhost(p, mesh.Ny, mesh.Nz, 
						bc.pType, bc.pBCIdx, bc.pBCcorIdx, mesh.dx, bc.pBCvalue);

				break;

			case 2: // y

				updYGhost(u, mesh.Nxu, mesh.Nzu, 
						bc.uType, bc.uBCIdx, bc.uBCcorIdx, mesh.dy, bc.uBCvalue);

				updYGhost(v, mesh.Nxv, mesh.Nzv, 
						bc.vType, bc.vBCIdx, bc.vBCcorIdx, 2*mesh.dy, bc.vBCvalue);

				updYGhost(w, mesh.Nxw, mesh.Nzw, 
						bc.wType, bc.wBCIdx, bc.wBCcorIdx, mesh.dy, bc.wBCvalue);

				updYGhost(p, mesh.Nx, mesh.Nz, 
						bc.pType, bc.pBCIdx, bc.pBCcorIdx, mesh.dy, bc.pBCvalue);

				break;

			case 3: // y

				updZGhost(u, mesh.Nxu, mesh.Nyu, 
						bc.uType, bc.uBCIdx, bc.uBCcorIdx, mesh.dy, bc.uBCvalue);

				updZGhost(v, mesh.Nxv, mesh.Nyv, 
						bc.vType, bc.vBCIdx, bc.vBCcorIdx, mesh.dy, bc.vBCvalue);

				updZGhost(w, mesh.Nxw, mesh.Nyw, 
						bc.wType, bc.wBCIdx, bc.wBCcorIdx, 2*mesh.dy, bc.wBCvalue);

				updZGhost(p, mesh.Nx, mesh.Ny, 
						bc.pType, bc.pBCIdx, bc.pBCcorIdx, mesh.dy, bc.pBCvalue);

				break;
		}
	}
	return 0;
}


int NSSolverEuler::PredictStep()
{

	for(int i=0; i<Nxu; ++i){
		for(int j=0; j<Nyu; ++j){
			for(int k=0; k<Nzu; ++k){
				u_str(i, j, k) = u(i, j, k) + 
					dt * (DiffusiveU(i, j, k) - ConvectU(i, j, k));
			}
		}
	}

	for(int i=0; i<Nxv; ++i){
		for(int j=0; j<Nyv; ++j){
			for(int k=0; k<Nzv; ++k){
				v_str(i, j, k) = v(i, j, k) + 
					dt * (DiffusiveV(i, j, k) - ConvectV(i, j, k));
			}
		}
	}

	for(int i=0; i<Nxw; ++i){
		for(int j=0; j<Nyw; ++j){
			for(int k=0; k<Nzw; ++k){
				w_str(i, j, k) = w(i, j, k) + 
					dt * (DiffusiveW(i, j, k) - ConvectW(i, j, k));
			}
		}
	}

	return 0;
}


double NSSolverEuler::ConvectU(int &i, int &j, int &k)
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


double NSSolverEuler::ConvectV(int &i, int &j, int &k)
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


double NSSolverEuler::ConvectW(int &i, int &j, int &k)
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


double NSSolverEuler::DiffusiveU(int &i, int &j, int &k)
{
	double Fdu;
	Fdu = fluid.nu * (
			(u(i+1, j, k) - 2 * u(i, j, k) + u(i-1, j, k)) / dx2 +
			(u(i, j+1, k) - 2 * u(i, j, k) + u(i, j-1, k)) / dy2 +
			(u(i, j, k+1) - 2 * u(i, j, k) + u(i, j, k-1)) / dz2 );
	return Fdu;
}


double NSSolverEuler::DiffusiveV(int &i, int &j, int &k)
{
	double Fdv;
	Fdv = fluid.nu * (
			(v(i+1, j, k) - 2 * v(i, j, k) + v(i-1, j, k)) / dx2 +
			(v(i, j+1, k) - 2 * v(i, j, k) + v(i, j-1, k)) / dy2 +
			(v(i, j, k+1) - 2 * v(i, j, k) + v(i, j, k-1)) / dz2 );
	return Fdv;
}


double NSSolverEuler::DiffusiveW(int &i, int &j, int &k)
{
	double Fdw;
	Fdw = fluid.nu * (
			(w(i+1, j, k) - 2 * w(i, j, k) + w(i-1, j, k)) / dx2 +
			(w(i, j+1, k) - 2 * w(i, j, k) + w(i, j-1, k)) / dy2 +
			(w(i, j, k+1) - 2 * w(i, j, k) + w(i, j, k-1)) / dz2 );
	return Fdw;
}



int updXGhost(Array3D<double> &u, int &Ny, int &Nz, 
		int &type, int &ghIdx, int &corIdx, double coeff, double &BCvalue)
{
	switch (type)
	{
		case 0:
			for(int j=0; j<Ny; ++j){
				for(int k=0; k<Nz; ++k){
					u(ghIdx, j, k) = u(corIdx, j, k); } }

			break;

		case 1: case -1:
			for(int j=0; j<Ny; ++j){
				for(int k=0; k<Nz; ++k){
					u(ghIdx, j, k) = 
						- type * u(corIdx, j, k) + coeff * BCvalue; } }

			break;

		default:
			throw invalid_argument("In EulerNSSsolver.cpp->updXGhost");
			break;
	}
	return 0;
}


int updYGhost(Array3D<double> &u, int &Nx, int &Nz, 
		int &type, int &ghIdx, int &corIdx, double coeff, double &BCvalue)
{
	switch (type)
	{
		case 0:
			for(int i=0; i<Nx; ++i){
				for(int k=0; k<Nz; ++k){
					u(i, ghIdx, k) = u(i, corIdx, k); } }

			break;

		case 1: case -1:
			for(int i=0; i<Nx; ++i){
				for(int k=0; k<Nz; ++k){
					u(i, ghIdx, k) = 
						- type * u(i, corIdx, k) + coeff * BCvalue; } }

			break;

		default:
			throw invalid_argument("In EulerNSSsolver.cpp->updXGhost");
			break;
	}
	return 0;
}


int updZGhost(Array3D<double> &u, int &Nx, int &Ny, 
		int &type, int &ghIdx, int &corIdx, double coeff, double &BCvalue)
{
	switch (type)
	{
		case 0:
			for(int i=0; i<Nx; ++i){
				for(int j=0; j<Ny; ++j){
					u(i, j, ghIdx) = u(i, j, corIdx); } }

			break;

		case 1: case -1:
			for(int i=0; i<Nx; ++i){
				for(int j=0; j<Ny; ++j){
					u(i, j, ghIdx) = 
						- type * u(i, j, corIdx) + coeff * BCvalue; } }

			break;

		default:
			throw invalid_argument("In EulerNSSsolver.cpp->updXGhost");
			break;
	}
	return 0;
}


int NSSolverEuler::test()
{
	int c = 1;	
	for(int i=0; i<mesh.Nxu; ++i){
		for(int j=0; j<mesh.Nyu; ++j){
			for(int k=0; k<mesh.Nzu; ++k){
				u(i, j, k) = c;
				++c;
			}
		}
	}

	c = 1;	
	for(int i=0; i<mesh.Nxv; ++i){
		for(int j=0; j<mesh.Nyv; ++j){
			for(int k=0; k<mesh.Nzv; ++k){
				v(i, j, k) = c;
				++c;
			}
		}
	}

	c = 1;	
	for(int i=0; i<mesh.Nxw; ++i){
		for(int j=0; j<mesh.Nyw; ++j){
			for(int k=0; k<mesh.Nzw; ++k){
				w(i, j, k) = c;
				++c;
			}
		}
	}

	return 0;
}




