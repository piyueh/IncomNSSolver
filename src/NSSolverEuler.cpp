# include "include/IncomNSSolver.h"
# include "NS_ElementryFuncs.cpp"
# include <cmath>


int NSSolver::InitSolver(CD t, CD Dt, CaryI3 pIdx, CD pR)
{
	dt = Dt; time = t;
	pRefIdx = pIdx; pRef = pR;

	dx2 = dx * dx; dy2 = dy * dy; dz2 = dz * dz;

	u.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	v.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	w.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

	Gu.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	Gv.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	Gw.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

	u.setZeros(); v.setZeros(); w.setZeros();

	TGVortex(Nxu, Nyu, Nzu, u, xu, yu,
			 Nxv, Nyv, Nzv, v, xv, yv, time);

	b.resize(Nx * Ny * Nz);
	p.resize(Nx * Ny * Nz);

	output("Init.txt");

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

	if (BCs[-2].get_vType() == 1)
	{
		for(int i=0; i<Nxv; ++i){
			for(int j=0; j<Nzv; ++j){
				v(i, 0, j) = BCs[-2].get_vBCvalue();
			}
		}
	}
	
	if (BCs[-3].get_wType() == 1)
	{
		for(int i=0; i<Nxw; ++i){
			for(int j=0; j<Nyw; ++j){
				w(i, j, 0) = BCs[-3].get_wBCvalue();
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

	if (BCs[2].get_vType() == 1)
	{
		for(int i=0; i<Nxv; ++i){
			for(int j=0; j<Nzv; ++j){
				v(i, Nyv-1, j) = BCs[2].get_vBCvalue();
			}
		}
	}
	
	if (BCs[3].get_wType() == 1)
	{
		for(int i=0; i<Nxw; ++i){
			for(int j=0; j<Nyw; ++j){
				w(i, j, Nzw-1) = BCs[3].get_wBCvalue();
			}
		}
	}


	
	InitLambda();

	return 0;
}


int NSSolver::solve(int targetNStep)
{
	pair<int, double> Itr;

	for(int n=0; n<targetNStep; ++n)
	{

		updateGhost();
		PredictStep(dt/3);
		updatePoissonSource(dt/3);
		Itr = pSolver.Solve(b, p);
		updateU(dt/3);

		updateGhost();
		PredictStep(15*dt/16, -5./9.);
		updatePoissonSource(5*dt/12);
		Itr = pSolver.Solve(b, p);
		updateU(5*dt/12);

		updateGhost();
		PredictStep(8*dt/15, -153./128.);
		updatePoissonSource(dt/4);
		Itr = pSolver.Solve(b, p);
		updateU(dt/4);

		p *= rho;

		time += dt;

		if (n % 50 == 0)
		{
			output(to_string(n)+".txt");
			cout << "n=" << n+1 << " ";
			cout << "time = " << time << " ";
			cout << Itr.first << " " << Itr.second << endl;
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


int NSSolver::PredictStep(CD & DT)
{
	tripleLoop(uBgIdx, uEdIdx, 0, Nyu, 0, Nzu, updGu);
	tripleLoop(0, Nxv, vBgIdx, vEdIdx, 0, Nzv, updGv);
	tripleLoop(0, Nxw, 0, Nyw, wBgIdx, wEdIdx, updGw);

	tripleLoop(uBgIdx, uEdIdx, 0, Nyu, 0, Nzu, DT, preU);
	tripleLoop(0, Nxv, vBgIdx, vEdIdx, 0, Nzv, DT, preV);
	tripleLoop(0, Nxw, 0, Nyw, wBgIdx, wEdIdx, DT, preW);

	return 0;
}


int NSSolver::PredictStep(CD & DT, CD & coef)
{
	tripleLoop(uBgIdx, uEdIdx, 0, Nyu, 0, Nzu, coef, updGu2);
	tripleLoop(0, Nxv, vBgIdx, vEdIdx, 0, Nzv, coef, updGv2);
	tripleLoop(0, Nxw, 0, Nyw, wBgIdx, wEdIdx, coef, updGw2);

	tripleLoop(uBgIdx, uEdIdx, 0, Nyu, 0, Nzu, DT, preU);
	tripleLoop(0, Nxv, vBgIdx, vEdIdx, 0, Nzv, DT, preV);
	tripleLoop(0, Nxw, 0, Nyw, wBgIdx, wEdIdx, DT, preW);

	return 0;
}


int NSSolver::updatePoissonSource(CD & DT)
{
	tripleLoop(0, Nx, 0, Ny, 0, Nz, DivOnPresPt);

	b /= DT;

	return 0;
}


int NSSolver::updateU(CD & DT)
{
	// update u
	//dualLoop(0, Nyu, 0, Nzu, updUB);
	tripleLoop(1, Nxu-1, 0, Nyu, 0, Nzu, DT, updU); 

	// update v
	//dualLoop(0, Nxv, 0, Nzv, updVB);
	tripleLoop(0, Nxv, 1, Nyv-1, 0, Nzv, DT, updV); 

	//update w
	//dualLoop(0, Nxw, 0, Nyw, updWB);
	tripleLoop(0, Nxw, 0, Nyw, 1, Nzw-1, DT, updW); 

	return 0;
}




