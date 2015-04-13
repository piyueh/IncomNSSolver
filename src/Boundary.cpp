# include "include/class_Boundary.h"


void set_correspondIdx(int &type, int &corIdx, int v1, int v0);


/*
 * Nx, Ny, Nz: numbers of cells in the three directions
 * N: number of cells on this boundary
 * dir: normal vector of the boundary, +-1, +-2, +-3 represent +-x, +-y, +-z
 * p, u, v, w: c++ pairs of boundary types and values
 * 			   the first element in the pair is type:
 * 		 			1: Dirichlet
 * 	 	 			0: Periodic
 * 					-1: Neumann
 * 			   the second element is value
 */
Boundary::Boundary(int Nx, int Ny, int Nz, unsigned int Dir, int Sign,
		pair<int, double> p, 
		pair<int, double> u, pair<int, double> v, pair<int, double> w)
{

	dir = Dir; sign = Sign;

	uType = u.first; vType = v.first; wType = w.first; pType = p.first;

	uBCvalue = u.second; vBCvalue = v.second; wBCvalue = w.second; pBCvalue = p.second;

	switch (dir)
	{
		case 1: // x
			switch (sign){
				case -1: // -x
					pBCIdx = -1; uBCIdx = -1; vBCIdx = -1; wBCIdx = -1;
					set_correspondIdx(pType, pBCcorIdx, 0, Nx-1);
					set_correspondIdx(uType, uBCcorIdx, 1, Nx-1);
					set_correspondIdx(vType, vBCcorIdx, 0, Nx-1);
					set_correspondIdx(wType, wBCcorIdx, 0, Nx-1);
					break;
				case 1: // +x
					pBCIdx = Nx; uBCIdx = Nx + 1; vBCIdx = Nx; wBCIdx = Nx;
					set_correspondIdx(pType, pBCcorIdx, Nx-1, 0);
					set_correspondIdx(uType, uBCcorIdx, Nx-1, 1);
					set_correspondIdx(vType, vBCcorIdx, Nx-1, 0);
					set_correspondIdx(wType, wBCcorIdx, Nx-1, 0);
					break;
			}
			break;

		case 2: // y
			switch (sign){
				case -1: // -y
					pBCIdx = -1; uBCIdx = -1; vBCIdx = -1; wBCIdx = -1;
					set_correspondIdx(pType, pBCcorIdx, 0, Ny-1);
					set_correspondIdx(uType, uBCcorIdx, 0, Ny-1);
					set_correspondIdx(vType, vBCcorIdx, 1, Ny-1);
					set_correspondIdx(wType, wBCcorIdx, 0, Ny-1);
					break;
				case 1: // +y
					pBCIdx = Ny; uBCIdx = Ny; vBCIdx = Ny + 1; wBCIdx = Ny;
					set_correspondIdx(pType, pBCcorIdx, Ny-1, 0);
					set_correspondIdx(uType, uBCcorIdx, Ny-1, 0);
					set_correspondIdx(vType, vBCcorIdx, Ny-1, 1);
					set_correspondIdx(wType, wBCcorIdx, Ny-1, 0);
					break;
			}
			break;

		case 3: // z
			switch (sign){
				case -1: // -z
					pBCIdx = -1; uBCIdx = -1; vBCIdx = -1; wBCIdx = -1;
					set_correspondIdx(pType, pBCcorIdx, 0, Nz-1);
					set_correspondIdx(uType, uBCcorIdx, 0, Nz-1);
					set_correspondIdx(vType, vBCcorIdx, 0, Nz-1);
					set_correspondIdx(wType, wBCcorIdx, 1, Nz-1);
					break;
				case 1: // +z
					pBCIdx = Nz; uBCIdx = Nz; vBCIdx = Nz; wBCIdx = Nz + 1;
					set_correspondIdx(pType, pBCcorIdx, Nz-1, 0);
					set_correspondIdx(uType, uBCcorIdx, Nz-1, 0);
					set_correspondIdx(vType, vBCcorIdx, Nz-1, 0);
					set_correspondIdx(wType, wBCcorIdx, Nz-1, 1);
					break;
			}
			break;

		default:
			throw invalid_argument("Invalid boundary direction");
			break;
	}
}


void set_correspondIdx(int &type, int &corIdx, int v1, int v0)
{
	switch (type) {
		case -1: case  1: corIdx = v1; break; // Neumann and Dirichlet
		case 0: corIdx = v0; break;      // Periodic
		default:
			throw invalid_argument("Error in Boundary.cpp");
			break;
	}
}


ostream &operator<<(ostream &os, Boundary &BC)
{
	os << "Direction: " << BC.dir*BC.sign << endl;

	os << "Type of pressure BC: " << BC.pType << endl;
	os << "Type of u BC: " << BC.uType << endl;
	os << "Type of v BC: " << BC.uType << endl;
	os << "Type of w BC: " << BC.uType << endl;

	os << "Value of pressure BC: " << BC.pBCvalue << endl;
	os << "Value of u BC: " << BC.uBCvalue << endl;
	os << "Value of v BC: " << BC.vBCvalue << endl;
	os << "Value of w BC: " << BC.wBCvalue << endl;

	os << "Index of ghost cell of pressure: " << BC.pBCIdx << endl;
	os << "Index of ghost cell of u: " << BC.uBCIdx << endl;
	os << "Index of ghost cell of v: " << BC.vBCIdx << endl;
	os << "Index of ghost cell of w: " << BC.wBCIdx << endl;

	os << "Index of the corresponding cell of pressure: "
	   << BC.pBCcorIdx << endl;
	os << "Index of the corresponding cell of u: "
	   << BC.uBCcorIdx << endl;
	os << "Index of the corresponding cell of v: "
	   << BC.vBCcorIdx << endl;
	os << "Index of the corresponding cell of w: "
	   << BC.wBCcorIdx << endl;
	return os;
}

