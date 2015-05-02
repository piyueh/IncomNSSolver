# include "include/IncomNSSolver.h"
# include "NS_ElementryFuncs.cpp"


NSSolver::NSSolver(Mesh &m, Fluid &f, Data &d, Solid &s, string &fName): 
	mesh(m), fluid(f), data(d), cyln(s)
{
	ifstream file(fName);
	string line;

	int tNStep = -1, ON = 100, SN = 1;
	double DT = -1; double pR = 0;
	array<int, 3> pRIdx = {0, 0, 0};

	while(getline(file, line))
	{
		istringstream OneLine(line);
		string var;

		OneLine >> var;

		if (var == "DT") { OneLine >> DT; }
		else if (var == "TargetNStep") { OneLine >> tNStep; }
		else if (var == "OutputNStep") { OneLine >> ON; }
		else if (var == "ScreenNStep") { OneLine >> SN; }
		else if (var == "RefPIdx") 
		{ 
			OneLine >> pRIdx[0] >> pRIdx[1] >> pRIdx[2]; 
		}
		else if (var == "RefP") { OneLine >> pR; }
		else if (var == "PTOL") { OneLine >> ptol; }
		else
		{
			throw invalid_argument(
					string("Invalid Argument in ") + fName + ": " + var);
		}
	}

	file.close();

	if (tNStep == -1) 
		throw invalid_argument(
				string("No Setting of total running time ") + 
				"steps is detected! Use the keyworld: TargetNStep");

	if (DT == -1) 
		throw invalid_argument(
				string("No Setting of time step is detected! ") +
				"Use the keyworld: DT");

	InitSolver(DT, tNStep, ON, SN, pRIdx, pR);

}


int NSSolver::InitSolver(CD & Dt, CI & tNStep, CI & ON, CI & SN, CaryI3 & pIdx, CD & pR)
{
	dt = Dt;
	targetNStep = tNStep; outputN = ON; screenN = SN;
	pRefIdx = pIdx; pRef = pR;

	dx2 = dx * dx; dy2 = dy * dy; dz2 = dz * dz;

	Gu.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	Gv.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	Gw.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

	dp.initShape(Nx, Ny, Nz); dp.setZeros();
	b.initShape(Nx, Ny, Nz); b.setZeros();

	pSolver.InitLinSys({Nx, Ny, Nz}, {dx, dy, dz});
	pSolver.setLHS(mesh.get_BCs());
	pSolver.setRefP(pRefIdx, pRef);
	pSolver.setTolerance(ptol);

	uBgIdx = (BCs[-1].get_uType() == 1) ? 1 : 0;
	uEdIdx = (BCs[1].get_uType() == 1) ? Nxu-1 : Nxu;
	vBgIdx = (BCs[-2].get_vType() == 1) ? 1 : 0;
	vEdIdx = (BCs[2].get_vType() == 1) ? Nyv-1 : Nyv;
	wBgIdx = (BCs[-3].get_wType() == 1) ? 1 : 0;
	wEdIdx = (BCs[3].get_wType() == 1) ? Nzw-1 : Nzw;
	
	InitLambda();

	return 0;
}


int NSSolver::solve()
{
	pair<int, double> Itr;

	clock_t t0, t;

	Map<VectorXd> dp_Eigen(dp.data(), dp.size());
	Map<VectorXd> b_Eigen(b.data(), b.size());

	for(int n=0; n<targetNStep; ++n)
	{
		t0 = clock();

		updateGhost();
		PredictStep(dt/3.);
		cyln.updVelocity(u, v, w);
		updatePoissonSource(1);
		Itr = pSolver.Solve(b_Eigen, dp_Eigen);
		updateField(dt/3.);

		updateGhost();
		PredictStep(15.*dt/16., -5./9.);
		cyln.updVelocity(u, v, w);
		updatePoissonSource(1);
		Itr = pSolver.Solve(b_Eigen, dp_Eigen);
		updateField(15.*dt/16);

		updateGhost();
		PredictStep(8.*dt/15., -153./128.);
		cyln.updVelocity(u, v, w);
		updatePoissonSource(1);
		Itr = pSolver.Solve(b_Eigen, dp_Eigen);
		updateField(8.*dt/15.);
		
		t = clock() - t0;

		time += dt;

		if (n % screenN == 0) 
		{
			cout << "n=" << n+1 << " ";
			cout << "time = " << time << " ";
			cout << Itr.first << " " << Itr.second << " ";
			cout << ((float)t)/CLOCKS_PER_SEC << endl;
		}

		if (n % outputN == 0) data.output(to_string(n)+".txt");
	}

	return 0;
}


int NSSolver::updateGhost()
{
	BCs[-3].updGhost(Nx, Ny, 0, p, dz);
	BCs[-3].updGhost(Nxu, Nyu, 1, u, dz);
	BCs[-3].updGhost(Nxv, Nyv, 2, v, dz);
	BCs[-3].updGhost(Nxw, Nyw, 3, w, dz);

	BCs[-2].updGhost(Nx, Nz, 0, p, dy);
	BCs[-2].updGhost(Nxu, Nzu, 1, u, dy);
	BCs[-2].updGhost(Nxv, Nzv, 2, v, dy);
	BCs[-2].updGhost(Nxw, Nzw, 3, w, dy);

	BCs[-1].updGhost(Ny, Nz, 0, p, dx);
	BCs[-1].updGhost(Nyu, Nzu, 1, u, dx);
	BCs[-1].updGhost(Nyv, Nzv, 2, v, dx);
	BCs[-1].updGhost(Nyw, Nzw, 3, w, dx);


	BCs[1].updGhost(Ny, Nz, 0, p, dx);
	BCs[1].updGhost(Nyu, Nzu, 1, u, dx);
	BCs[1].updGhost(Nyv, Nzv, 2, v, dx);
	BCs[1].updGhost(Nyw, Nzw, 3, w, dx);

	BCs[2].updGhost(Nx, Nz, 0, p, dy);
	BCs[2].updGhost(Nxu, Nzu, 1, u, dy);
	BCs[2].updGhost(Nxv, Nzv, 2, v, dy);
	BCs[2].updGhost(Nxw, Nzw, 3, w, dy);

	BCs[3].updGhost(Nx, Ny, 0, p, dz);
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

	return 0;
}


int NSSolver::updateField(CD & coeff)
{
	// update u
	tripleLoop(1, Nxu-1, 0, Nyu, 0, Nzu, updU); 

	// update v
	tripleLoop(0, Nxv, 1, Nyv-1, 0, Nzv, updV); 

	//update w
	tripleLoop(0, Nxw, 0, Nyw, 1, Nzw-1, updW); 

	// update p
	tripleLoop(0, Nx, 0, Ny, 0, Nz, coeff, updP);

	return 0;
}




