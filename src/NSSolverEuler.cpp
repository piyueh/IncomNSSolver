# include "include/IncomNSSolver.h"
# include "NS_Convection.cpp"
# include "NS_Diffusion.cpp"


int NSSolver::InitSolver(CD t, CD Dt, CaryI3 pIdx, CD pR)
{
	dt = Dt; time = t;
	pRefIdx = pIdx; pRef = pR;

	dx2 = dx * dx; dy2 = dy * dy; dz2 = dz * dz;

	u.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	v.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	w.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

	u_str.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	v_str.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	w_str.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

	u.setZeros(); v.setZeros(); w.setZeros();

	b.resize(Nx * Ny * Nz);
	p.resize(Nx * Ny * Nz);

	pSolver.InitLinSys({Nx, Ny, Nz}, {dx, dy, dz});
	pSolver.setLHS(mesh.get_BCs());
	pSolver.setRefP(pRefIdx, pRef);

	uBgIdx = (BCs[-1].get_uType() == 1) ? 1 : 0;
	uEdIdx = (BCs[1].get_uType() == 1) ? Nxu-1 : Nxu;
	vBgIdx = (BCs[-2].get_vType() == 1) ? 1 : 0;
	vEdIdx = (BCs[2].get_vType() == 1) ? Nyv-1 : Nyv;
	wBgIdx = (BCs[-3].get_wType() == 1) ? 1 : 0;
	wEdIdx = (BCs[3].get_wType() == 1) ? Nzw-1 : Nzw;

	if (BCs[-1].get_uType() == 1)
	{
		for(int i=0; i<Nyu; ++i){
			for(int j=0; j<Nzu; ++j){
				u(0, i, j) = BCs[-1].get_uBCvalue();
			}
		}
	}

	if (BCs[-2].get_uType() == 1)
	{
		for(int i=0; i<Nxv; ++i){
			for(int j=0; j<Nzv; ++j){
				v(i, 0, j) = BCs[-2].get_uBCvalue();
			}
		}
	}
	
	if (BCs[-3].get_uType() == 1)
	{
		for(int i=0; i<Nxw; ++i){
			for(int j=0; j<Nyw; ++j){
				w(i, j, 0) = BCs[-3].get_uBCvalue();
			}
		}
	}

	if (BCs[1].get_uType() == 1)
	{
		for(int i=0; i<Nyu; ++i){
			for(int j=0; j<Nzu; ++j){
				u(Nxu-1, i, j) = BCs[1].get_uBCvalue();
			}
		}
	}

	if (BCs[2].get_uType() == 1)
	{
		for(int i=0; i<Nxv; ++i){
			for(int j=0; j<Nzv; ++j){
				v(i, Nyv-1, j) = BCs[2].get_uBCvalue();
			}
		}
	}
	
	if (BCs[3].get_uType() == 1)
	{
		for(int i=0; i<Nxw; ++i){
			for(int j=0; j<Nyw; ++j){
				w(i, j, Nzw-1) = BCs[3].get_uBCvalue();
			}
		}
	}
	return 0;
}


int NSSolver::solve(int targetNStep)
{

	for(int n=0; n<targetNStep; ++n)
	{

		updateGhost();

		PredictStep();

		updatePoissonSource();

		pSolver.Solve(b, p);

		updateU();

		time += dt;

		if (n % 50 == 0)
		{
			output(to_string(n)+".txt");
			cout << "n=" << n+1 << " ";
			cout << "time = " << time << endl;
		}
	}

	return 0;
}


int NSSolver::updateGhost()
{
	//BCs[-3].updGhost(Nx, Ny, 0, p, dz);
	BCs[-3].updGhost(Nxu, Nyu, 1, u, dz);
	BCs[-3].updGhost(Nxv, Nyv, 2, v, dz);
	BCs[-3].updGhost(Nxw, Nyw, 3, w, dz);

	//BCs[-2].updGhost(Nx, Nz, 0, p, dy);
	BCs[-2].updGhost(Nxu, Nzu, 1, u, dy);
	BCs[-2].updGhost(Nxv, Nzv, 2, v, dy);
	BCs[-2].updGhost(Nxw, Nzw, 3, w, dy);

	//BCs[-1].updGhost(Ny, Nz, 0, p, dx);
	BCs[-1].updGhost(Nyu, Nzu, 1, u, dx);
	BCs[-1].updGhost(Nyv, Nzv, 2, v, dx);
	BCs[-1].updGhost(Nyw, Nzw, 3, w, dx);


	//BCs[1].updGhost(Ny, Nz, 0, p, dx);
	BCs[1].updGhost(Nyu, Nzu, 1, u, dx);
	BCs[1].updGhost(Nyv, Nzv, 2, v, dx);
	BCs[1].updGhost(Nyw, Nzw, 3, w, dx);

	//BCs[2].updGhost(Nx, Nz, 0, p, dy);
	BCs[2].updGhost(Nxu, Nzu, 1, u, dy);
	BCs[2].updGhost(Nxv, Nzv, 2, v, dy);
	BCs[2].updGhost(Nxw, Nzw, 3, w, dy);

	//BCs[3].updGhost(Nx, Ny, 0, p, dz);
	BCs[3].updGhost(Nxu, Nyu, 1, u, dz);
	BCs[3].updGhost(Nxv, Nyv, 2, v, dz);
	BCs[3].updGhost(Nxw, Nyw, 3, w, dz);


	return 0;
}


int NSSolver::PredictStep()
{

	for(int i=uBgIdx; i<uEdIdx; ++i){
		for(int j=0; j<Nyu; ++j){
			for(int k=0; k<Nzu; ++k){
				u_str(i, j, k) = u(i, j, k) + 
					dt * (DiffusiveU(i, j, k) - ConvectU(i, j, k));
			}
		}
	}

	for(int i=0; i<Nxv; ++i){
		for(int j=vBgIdx; j<vEdIdx; ++j){
			for(int k=0; k<Nzv; ++k){
				v_str(i, j, k) = v(i, j, k) + 
					dt * (DiffusiveV(i, j, k) - ConvectV(i, j, k));
			}
		}
	}

	for(int i=0; i<Nxw; ++i){
		for(int j=0; j<Nyw; ++j){
			for(int k=wBgIdx; k<wEdIdx; ++k){
				w_str(i, j, k) = w(i, j, k) + 
					dt * (DiffusiveW(i, j, k) - ConvectW(i, j, k));
			}
		}
	}

	return 0;
}


int NSSolver::updatePoissonSource()
{
	for(int i=0; i<Nx; ++i){
		for(int j=0; j<Ny; ++j){
			for(int k=0; k<Nz; ++k){
				
				b(i * Nyz + j * Nz + k) = 
					((u_str(i+1, j, k) - u_str(i, j, k)) / dx + 
					 (v_str(i, j+1, k) - v_str(i, j, k)) / dy + 
					 (w_str(i, j, k+1) - w_str(i, j, k)) / dz) / dt;
			}
		}
	}

	return 0;
}


int NSSolver::updateU()
{
	auto calIdx = [this] (CI & i, CI & j, CI & k) -> int
				  { return i * Nyz + j * Nz + k; };

	
	// update u
	auto fuB = [this] (CI & i, CI & j) -> void 
			   { u(0, i, j) = u_str(0, i, j); 
				 u(Nxu-1, i, j) = u_str(Nxu-1, i, j); };

	auto fu = [&calIdx, this] (CI & i, CI & j, CI & k) -> void
			  { int tmp1 = calIdx(i, j, k); int tmp2 = tmp1 - Nyz; 
				u(i, j, k) = u_str(i, j, k) - dt * (p(tmp1) - p(tmp2)) / dx; };

	dualLoop(0, Nyu, 0, Nzu, fuB);

	tripleLoop(1, Nxu-1, 0, Nyu, 0, Nzu, fu); 



	// update v
	auto fvB = [this] (CI & i, CI & j) -> void 
			   { v(i, 0, j) = v_str(i, 0, j); 
				 v(i, Nyv-1, j) = v_str(i, Nyv-1, j); };

	auto fv = [&calIdx, this] (CI & i, CI & j, CI & k) -> void
			  { int tmp1 = calIdx(i, j, k); int tmp2 = tmp1 - Nz; 
				v(i, j, k) = v_str(i, j, k) - dt * (p(tmp1) - p(tmp2)) / dy; };

	dualLoop(0, Nxv, 0, Nzv, fvB);

	tripleLoop(0, Nxv, 1, Nyv-1, 0, Nzv, fv); 


	//update w
	auto fwB = [this] (CI & i, CI & j) -> void 
			   { w(i, j, 0) = w_str(i, j, 0); 
				 w(i, j, Nzw-1) = w_str(i, j, Nzw-1); };

	auto fw = [&calIdx, this] (CI & i, CI & j, CI & k) -> void
			  { int tmp1 = calIdx(i, j, k); int tmp2 = tmp1 - 1; 
				w(i, j, k) = w_str(i, j, k) - dt * (p(tmp1) - p(tmp2)) / dz; };

	dualLoop(0, Nxw, 0, Nyw, fwB);

	tripleLoop(0, Nxw, 0, Nyw, 1, Nzw-1, fw); 

	return 0;
}




