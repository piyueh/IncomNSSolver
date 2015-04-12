# include "include/IncomNSSolver.h"


/*
 * Nx, Ny, Nz: numbers of cells in the three directions
 * N: number of cells on this boundary
 * dir: normal vector of the boundary, toward outside of the  domain
 * pt: BC type of pressure
 * 		 1: Dirichlet
 * 	 	 0: Periodic
 * 		-1: Neumann
 * pv: the value of the boundary condition
 * vt: BC type of velocity, using 1, 0, -1, like pt
 * vv: value of velocity BC
 */
Boundary::Boundary(int Nx, int Ny, int Nz, 
		string dir, int pt, double pv, int vt, double vv)
{

	direction = dir;
	pType = pt;
	vType = vt;
	pBCvalue = pv;
	vBCvalue = vv;

	switch (dir[1])
	{
		case 'x':
			Ncells = Ny * Nz;
			for(int j=0; j<Ny; ++j){
				for(int k=0; k<Nz; ++k){
					Cell.push_back(j*Nz+k);
					OppCell.push_back((Nx-1)*Ny*Nz+j*Nz+k);
				}
			}
			break;
		
		case 'y':
			Ncells = Nx * Nz;
			for(int i=0; i<Nx; ++i){
				for(int k=0; k<Nz; ++k){
					Cell.push_back(i*Ny*Nz+k);
					OppCell.push_back(i*Ny*Nz+(Ny-1)*Nz+k);
				}
			}
			break;

		case 'z':
			Ncells = Nx * Ny;
			for(int i=0; i<Nx; ++i){
				for(int j=0; j<Ny; ++j){
					Cell.push_back(i*Ny*Nz+j*Nz);
					OppCell.push_back(i*Ny*Nz+j*Nz+Nz-1);
				}
			}
			break;

		default:
			throw invalid_argument("Invalid boundary direction");
			break;
	}

	switch (dir[0])
	{

		case '-': break;
		case '+':
			Cell.swap(OppCell);
			break;

		default:
			throw invalid_argument("Invalid boundary direction");
			break;
	}
}


