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

	return 0;
}




