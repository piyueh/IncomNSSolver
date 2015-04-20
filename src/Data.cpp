# include "include/IncomNSSolver.h"

int Data::InitData(Mesh & mesh)
{

	map<int, Boundary> & BCs = mesh.get_BCs();
	const int &Nxu=mesh.get_Nxu(), &Nyu=mesh.get_Nyu(), &Nzu=mesh.get_Nzu();
	const int &Nxv=mesh.get_Nxv(), &Nyv=mesh.get_Nyv(), &Nzv=mesh.get_Nzv();
	const int &Nxw=mesh.get_Nxw(), &Nyw=mesh.get_Nyw(), &Nzw=mesh.get_Nzw();

	time = 0;

	Nx = mesh.get_Nx(); Ny = mesh.get_Ny(); Nz = mesh.get_Nz();

	u.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	v.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	w.initShape(-1, Nxw, -1, Nyw, -1, Nzw);

	u.setZeros(); v.setZeros(); w.setZeros();

	p.resize(Nx * Ny * Nz); p.setZero();

	function<void(CI &, CI &)> f;

	f = [&, this] (CI & i, CI & j) 
	{ u(0, i, j) = BCs[-1].get_uBCvalue();} ;

	if (BCs[-1].get_uType() == 1) dualLoop(0, Nyu, 0, Nzu, f);

	f = [&, this] (CI & i, CI & j) 
	{ v(i, 0, i) = BCs[-2].get_vBCvalue();} ;

	if (BCs[-2].get_vType() == 1) dualLoop(0, Nxv, 0, Nzv, f);
	
	f = [&, this] (CI & i, CI & j) 
	{ w(i, j, 0) = BCs[-3].get_wBCvalue();} ;

	if (BCs[-3].get_wType() == 1) dualLoop(0, Nxw, 0, Nyw, f);

	f = [&, this] (CI & i, CI & j) 
	{ u(Nxu-1, i, j) = BCs[1].get_uBCvalue();} ;

	if (BCs[1].get_uType() == 1) dualLoop(0, Nyu, 0, Nzu, f);

	f = [&, this] (CI & i, CI & j) 
	{ v(i, Nyv-1, i) = BCs[2].get_vBCvalue();} ;

	if (BCs[2].get_vType() == 1) dualLoop(0, Nxv, 0, Nzv, f);
	
	f = [&, this] (CI & i, CI & j) 
	{ w(i, j, Nzw-1) = BCs[3].get_wBCvalue();} ;

	if (BCs[3].get_wType() == 1) dualLoop(0, Nxw, 0, Nyw, f);

	return 0;
}

