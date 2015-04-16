# include "include/IncomNSSolver.h"


int set_correspondIdx(int &type, int &corIdx, int v1, int v0);


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
Boundary::Boundary(const array<int, 3> N, const unsigned int Dir, 
		const int Sign, const pair<int, double> p, 
		const pair<int, double> u, const pair<int, double> v, const pair<int, double> w)
{
	auto & Nx = N[0], & Ny = N[1], & Nz = N[2];

	dir = Dir; sign = Sign;

	uType = u.first; vType = v.first; wType = w.first; pType = p.first;

	uBCvalue = u.second; vBCvalue = v.second; 
	wBCvalue = w.second; pBCvalue = p.second;

	if ((uType == -1) && (sign == -1)) uBCvalue *= -1;
	if ((vType == -1) && (sign == -1)) vBCvalue *= -1;
	if ((wType == -1) && (sign == -1)) wBCvalue *= -1;
	if ((pType == -1) && (sign == -1)) pBCvalue *= -1;

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


int set_correspondIdx(int &type, int &corIdx, int v1, int v0)
{
	switch (type) {
		case -1: case  1: corIdx = v1; break; // Neumann and Dirichlet
		case 0: corIdx = v0; break;      // Periodic
		default:
			throw invalid_argument("Error in Boundary.cpp");
			break;
	}
	return 0;
}


