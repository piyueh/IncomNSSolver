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


Boundary::Boundary(array<array<int, 3>, 4> & N, CUI Dir, CI Sign, 
		CPairID p, CPairID u, CPairID v, CPairID w)
{
	dir = Dir; sign = Sign;

	switch (dir)
	{
		case 1: corDir[0] = 1; corDir[1] = 2; break;
		case 2: corDir[0] = 0; corDir[1] = 2; break;
		case 3: corDir[0] = 0; corDir[1] = 1; break;
	}


	uType = u.first; vType = v.first; wType = w.first; pType = p.first;

	pBCvalue = ((sign == -1) && (pType == -1)) ? - p.second : p.second;
	uBCvalue = ((sign == -1) && (uType == -1)) ? - u.second : u.second;
	vBCvalue = ((sign == -1) && (vType == -1)) ? - v.second : v.second;
	wBCvalue = ((sign == -1) && (wType == -1)) ? - w.second : w.second;


	/************************************************************************
	* this function determine the index of the ghost cell of primative variables
	* there are no ghost cells for the velocity normal to the boundary
	************************************************************************/
	auto Idx = [&N, &Dir, &Sign] (CI & i) -> int
	{
		return (Sign == -1) ? -1 : N[i][Dir-1];	
	};


	/************************************************************************
	* this function determine the index of corresponding cells, which will be 
	* used during the process of updating ghost cells, of the ghost cells
	* on the current boundary
	************************************************************************/
	auto corIdx = [&N, &Dir, &Sign] (CI & i, CI & type) -> int
	{
		return (i == Dir) ? 
			(((Sign == 1) ? (type == 0) : (type != 0)) ? 1 : N[i][Dir-1]-2) : 
			(((Sign == 1) ? (type == 0) : (type != 0)) ? 0 : N[i][Dir-1]-1) ;
	};

	
	/************************************************************************
	* determine the index and coresponding index
	************************************************************************/
	for(int i=0; i<4; ++i)
	{
		BCIdx[i] = Idx(i);
		BCcorIdx[i] = corIdx(i, BCtype[i]); 
	}

	InitUpdGh();
}


int Boundary::InitUpdGh()
{

	Array3D<function<void(CI &, CI &, CI &, aA3Dd_ptr_4 &, CaryD3 &)>> f;
	f.initShape(1, 3, -1, 1, 0, 0);


	/************************************************************************
	* create a pool of functions that update the ghost cells according to
	* different boundary conditions
	* for the normal velocity, it is not ghost cell, instead, it updates the 
	* velocity nodes which are located right on the boundary
	************************************************************************/

	// ---------------------------x direction------------------------------//
	// x direction, Neumann
	f(1, -1, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(BCIdx[id], i, j) = 
				(*F[id])(BCcorIdx[id], i, j) + dL[dir-1] * BCvalues[id]; 
		};


	// x direction, Periodic
	f(1, 0, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(BCIdx[id], i, j) = (*F[id])(BCcorIdx[id], i, j); 
		};


	// x direction, Dirichlet
	f(1, 1, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(BCIdx[id], i, j) = 
				- (*F[id])(BCcorIdx[id], i, j) + 2 * BCvalues[id];
		};

	



	// ---------------------------y direction------------------------------//
	// y direction, Neumann
	f(2, -1, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(i, BCIdx[id], j) = 
				(*F[id])(i, BCcorIdx[id], j) + dL[dir-1] * BCvalues[id]; 
		};


	// y direction, Periodic
	f(2, 0, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(i, BCIdx[id], j) = (*F[id])(i, BCcorIdx[id], j); 
		};


	// y direction, Dirichlet
	f(2, 1, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(i, BCIdx[id], j) = 
				- (*F[id])(i, BCcorIdx[id], j) + 2 * BCvalues[id];
		};




	// ---------------------------z direction------------------------------//
	// z direction, Neumann
	f(3, -1, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(i, j, BCIdx[id]) = 
				(*F[id])(i, j, BCcorIdx[id]) + dL[dir-1] * BCvalues[id]; 
		};

	// z direction, Periodic
	f(3, 0, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(i, j, BCIdx[id]) = (*F[id])(i, j, BCcorIdx[id]); 
		};


	// z direction, Dirichlet
	f(3, 1, 0) = 
		[this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		{ 
			(*F[id])(i, j, BCIdx[id]) = 
				- (*F[id])(i, j, BCcorIdx[id]) + 2 * BCvalues[id];
		};




	// -----------pool for the case the normal velocity is Neumann--------//
	vector<function<void(CI &, CI &, CI &, aA3Dd_ptr_4 &, CaryD3 &)>> h;
	h.resize(3);

	// fix the u component => directly assign BC value to the boundary nodes
	h[0] = [this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		   { 
			   (*F[1])(BCIdx[1], i, j) = 
					(*F[1])(BCcorIdx[1], i, j) + 2 * dL[0] * BCvalues[id]; 
		   };

	// fix the v component => directly assign BC value to the boundary nodes
	h[1] = [this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		   { 
			   (*F[2])(i, BCIdx[2], j) = 
					(*F[2])(i, BCcorIdx[2], j) + 2 * dL[1] * BCvalues[id]; 
		   };

	// fix the w component => directly assign BC value to the boundary nodes
	h[2] = [this] (CI & i, CI & j, CI & id, aA3Dd_ptr_4 & F, CaryD3 & dL) -> void
		   { 
			   (*F[3])(i, j, BCIdx[3]) = 
					(*F[3])(i, j, BCcorIdx[3]) + 2 * dL[2] * BCvalues[id]; 
		   };

	// finally, determine the corresponding functions
	for(int i=0; i<4; ++i) updOneGh[i] = f(dir, BCtype[i], 0);	

	if (BCtype[dir] == -1) updOneGh[dir] = h[dir-1];



	return 0;
}


int Boundary::updGhost(array<array<int, 3>, 4> & N, aA3Dd_ptr_4 & F, CaryD3 & dL)
{


	for(int id=0; id<4; ++id)
	{	
		for(int i=0; i<N[id][corDir[0]]; ++i)
		{
			for(int j=0; j<N[id][corDir[1]]; ++j)
			{
				updOneGh[id](i, j, id, F, dL);
			}
		}
	}

	return 0;
}


int Boundary::updDirichBC(array<array<int, 3>, 4> & N, aA3Dd_ptr_4 & F)
{

	// -----------pool for the case the normal velocity is dirichlet--------//
	map<int, function<void(CI &, CI &, A3Dd &)>> g;

	// fix the u component => directly assign BC value to the boundary nodes
	g[1] = [this] (CI & i, CI & j, A3Dd & A) -> void
		   { A(0, i, j) = BCvalues[1]; };

	// fix the v component => directly assign BC value to the boundary nodes
	g[2] = [this] (CI & i, CI & j, A3Dd & A) -> void
		   { A(i, 0, j) = BCvalues[2]; };

	// fix the w component => directly assign BC value to the boundary nodes
	g[3] = [this] (CI & i, CI & j, A3Dd & A) -> void
		   { A(i, j, 0) = BCvalues[3]; };

	// fix the u component => directly assign BC value to the boundary nodes
	g[-1] = [this, & N] (CI & i, CI & j, A3Dd & A) -> void
		   { A(N[1][0]-1, i, j) = BCvalues[1]; };

	// fix the v component => directly assign BC value to the boundary nodes
	g[-2] = [this, & N] (CI & i, CI & j, A3Dd & A) -> void
		   { A(i, N[2][1]-1, j) = BCvalues[2]; };

	// fix the w component => directly assign BC value to the boundary nodes
	g[-3] = [this, & N] (CI & i, CI & j, A3Dd & A) -> void
		   { A(i, j, N[3][2]-1) = BCvalues[3]; };


	if (BCtype[dir] == 1)
	{
		int id = dir * sign;

		for(int i=0; i<N[dir][corDir[0]]; ++i)
		{
			for(int j=0; j<N[dir][corDir[1]]; ++j) g[id](i, j, *F[dir]);
		}
	}
			

	return 0;
}


Boundary& Boundary::operator=(const Boundary & A)
{

	dir = A.dir; sign = A.sign;
	corDir = A.corDir;
	BCtype = A.BCtype;
	BCIdx = A.BCIdx;
	BCcorIdx = A.BCcorIdx;
	BCvalues = A.BCvalues;

	InitUpdGh();

	return *this;
}


