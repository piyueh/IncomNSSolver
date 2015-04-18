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


Boundary::Boundary(CaryI3 N, CUI Dir, CI Sign, 
		CPairID p, CPairID u, CPairID v, CPairID w)
{
	auto & Nx = N[0], & Ny = N[1], & Nz = N[2];

	auto Idx = [&N, &Dir, &Sign] (CI & i) -> int
		{ return (Sign == -1) ? -1 : (i == Dir) ? N[Dir-1]+1 : N[Dir-1]; };	

	auto corIdx = [&N, &Dir, &Sign] (CI & i, CI & type) -> int
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

	InitUpdGh();
}


int Boundary::InitUpdGh()
{

	Array3D<function<void(CI &, CI &, A3Dd &, CD &)>> f;
	f.initShape(1, 3, -1, 1, 0, 3);


	f(1, -1, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(pBCIdx, i, j) =  A(pBCcorIdx, i, j) + dL * pBCvalue; };
	f(1, -1, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(uBCIdx, i, j) =  A(uBCcorIdx, i, j) + 2 * dL * uBCvalue; };
	f(1, -1, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(vBCIdx, i, j) =  A(vBCcorIdx, i, j) + dL * vBCvalue; };
	f(1, -1, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(wBCIdx, i, j) =  A(wBCcorIdx, i, j) + dL * wBCvalue; };


	f(1, 0, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(pBCIdx, i, j) = A(pBCcorIdx, i, j); };
	f(1, 0, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(uBCIdx, i, j) = A(uBCcorIdx, i, j); };
	f(1, 0, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(vBCIdx, i, j) = A(vBCcorIdx, i, j); };
	f(1, 0, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(wBCIdx, i, j) = A(wBCcorIdx, i, j); };


	f(1, 1, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(pBCIdx, i, j) =  - A(pBCcorIdx, i, j) + 2 * pBCvalue; };
	f(1, 1, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(uBCIdx, i, j) =  - A(uBCcorIdx, i, j) + 2 * uBCvalue; };
	f(1, 1, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(vBCIdx, i, j) =  - A(vBCcorIdx, i, j) + 2 * vBCvalue; };
	f(1, 1, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(wBCIdx, i, j) =  - A(wBCcorIdx, i, j) + 2 * wBCvalue; };



	f(2, -1, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, pBCIdx, j) =  A(i, pBCcorIdx, j) + dL * pBCvalue; };
	f(2, -1, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, uBCIdx, j) =  A(i, uBCcorIdx, j) + dL * uBCvalue; };
	f(2, -1, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, vBCIdx, j) =  A(i, vBCcorIdx, j) + 2 * dL * vBCvalue; };
	f(2, -1, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, wBCIdx, j) =  A(i, wBCcorIdx, j) + dL * wBCvalue; };


	f(2, 0, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, pBCIdx, j) = A(i, pBCcorIdx, j); };
	f(2, 0, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, uBCIdx, j) = A(i, uBCcorIdx, j); };
	f(2, 0, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, vBCIdx, j) = A(i, vBCcorIdx, j); };
	f(2, 0, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, wBCIdx, j) = A(i, wBCcorIdx, j); };


	f(2, 1, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, pBCIdx, j) =  - A(i, pBCcorIdx, j) + 2 * pBCvalue; };
	f(2, 1, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, uBCIdx, j) =  - A(i, uBCcorIdx, j) + 2 * uBCvalue; };
	f(2, 1, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, vBCIdx, j) =  - A(i, vBCcorIdx, j) + 2 * vBCvalue; };
	f(2, 1, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, wBCIdx, j) =  - A(i, wBCcorIdx, j) + 2 * wBCvalue; };



	f(3, -1, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, pBCIdx) =  A(i, j, pBCcorIdx) + dL * pBCvalue; };
	f(3, -1, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, uBCIdx) =  A(i, j, uBCcorIdx) + dL * uBCvalue; };
	f(3, -1, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, vBCIdx) =  A(i, j, vBCcorIdx) + dL * vBCvalue; };
	f(3, -1, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, wBCIdx) =  A(i, j, wBCcorIdx) + 2 * dL * wBCvalue; };


	f(3, 0, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, pBCIdx) = A(i, j, pBCcorIdx); };
	f(3, 0, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, uBCIdx) = A(i, j, uBCcorIdx); };
	f(3, 0, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, vBCIdx) = A(i, j, vBCcorIdx); };
	f(3, 0, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, wBCIdx) = A(i, j, wBCcorIdx); };


	f(3, 1, 0) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, pBCIdx) =  - A(i, j, pBCcorIdx) + 2 * pBCvalue; };
	f(3, 1, 1) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, uBCIdx) =  - A(i, j, uBCcorIdx) + 2 * uBCvalue; };
	f(3, 1, 2) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, vBCIdx) =  - A(i, j, vBCcorIdx) + 2 * vBCvalue; };
	f(3, 1, 3) = [this] (CI & i, CI & j, A3Dd & A, CD & dL) -> void
					{ A(i, j, wBCIdx) =  - A(i, j, wBCcorIdx) + 2 * wBCvalue; };


	for(int i=0; i<4; ++i) updOneGh[i] = f(dir, BCtype[i], i);	

	return 0;
}


int Boundary::updGhost(CI & N1, CI & N2, CI & field, A3Dd & A, CD & dL)
{
	for(int i=0; i<N1; ++i)
	{
		for(int j=0; j<N2; ++j)
		{
			updOneGh[field](i, j, A, dL);
		}
	}

	return 0;
}


Boundary& Boundary::operator=(const Boundary & A)
{
	dir = A.dir; sign = A.sign;
	BCtype = A.BCtype;
	BCIdx = A.BCIdx;
	BCcorIdx = A.BCcorIdx;
	BCvalues = A.BCvalues;

	InitUpdGh();

	return *this;
}


