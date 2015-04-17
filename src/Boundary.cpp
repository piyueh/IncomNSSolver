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

# include "include/IncomNSSolver.h"


Boundary::Boundary(const array<int, 3> N, const unsigned int Dir, 
	const int Sign, const pair<int, double> p, 
	const pair<int, double> u, const pair<int, double> v, const pair<int, double> w)
{
	auto & Nx = N[0], & Ny = N[1], & Nz = N[2];

	auto Idx = [&N, &Dir, &Sign] (const int & i) -> int
		{ return (Sign == -1) ? -1 : (i == Dir) ? N[Dir-1]+1 : N[Dir-1]; };	

	auto corIdx = [&N, &Dir, &Sign] (const int & i, const int & type) -> int
		{ return ((Sign == 1) ? (type != 0) : (type == 0)) ? 
								N[Dir-1] - 1 : (i == Dir) ? 1 : 0; };	
	
	dir = Dir; sign = Sign;

	uType = u.first; vType = v.first; wType = w.first; pType = p.first;

	uBCvalue = ((sign == -1) && (uType == -1)) ? - u.second : u.second;
	vBCvalue = ((sign == -1) && (vType == -1)) ? - v.second : v.second;
	wBCvalue = ((sign == -1) && (wType == -1)) ? - w.second : w.second;
	pBCvalue = ((sign == -1) && (pType == -1)) ? - p.second : p.second;

	for(int i=0; i<4; ++i)
	{
		BCIdx[i] = Idx(i);
		BCcorIdx[i] = corIdx(i, BCtype[i]); 
	}

}


Boundary& Boundary::operator=(const Boundary & A)
{
	dir = A.dir; sign = A.sign;
	BCtype = A.BCtype;
	BCIdx = A.BCIdx;
	BCcorIdx = A.BCcorIdx;
	BCvalues = A.BCvalues;

	return *this;
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

