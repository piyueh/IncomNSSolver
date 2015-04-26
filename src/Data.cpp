# include "include/IncomNSSolver.h"

Data::Data(string & fName)
{
	ifstream file(fName);
	string line;

	while (getline(file, line))
	{
		istringstream OneLine(line);
		string var;
		
		int N1, N2, N3;

		OneLine >> var;
		if (var == "TIME") {OneLine >> time;}
		else if (var == "Nu") 
		{
			OneLine >> N1 >> N2 >> N3;
			u.initShape(-1, N1-2, -1, N2-2, -1, N3-2);
		}
		else if (var == "Nv") 
		{
			OneLine >> N1 >> N2 >> N3;
			v.initShape(-1, N1-2, -1, N2-2, -1, N3-2);
		}
		else if (var == "Nw") 
		{
			OneLine >> N1 >> N2 >> N3;
			w.initShape(-1, N1-2, -1, N2-2, -1, N3-2);
		}
		else if (var == "Np") 
		{
			OneLine >> Nx >> Ny >> Nz;
			p.initShape(-1, N1-2, -1, N2-2, -1, N3-2);
		}
		else if ((var.empty()) || (var == "u") || 
				(var == "v") || (var == "w") || (var == "p")) {}
		else
			throw invalid_argument(
					string("Invalid Argument in ") + fName + ": " + var);
	}
	file.close();

	file.open(fName);
	while (getline(file, line))
	{
		istringstream OneLine(line);
		string var;
		double value;
		
		OneLine >> var;
		if (var == "u") 
		{ 
			Array3D<double> tmp;
			while (OneLine >> value) tmp.push_back(value);
			if (tmp.size() == u.size()) { u = tmp; }
			else if (tmp.size() == 1) { u.setConstant(tmp[0]); }
			else 
			{ 
				throw invalid_argument(string("Size of input u is wrong! ") +
						string("tmp.size=") + to_string(tmp.size()) + " while " + 
						string("u.size=") + to_string(u.size())); 
			}
		}
		else if (var == "v") 
		{ 
			Array3D<double> tmp;
			while (OneLine >> value) tmp.push_back(value);
			if (tmp.size() == v.size()) { v = tmp; }
			else if (tmp.size() == 1) { v.setConstant(tmp[0]); }
			else
			{ 
				throw invalid_argument(string("Size of input v is wrong! ") +
						string("tmp.size=") + to_string(tmp.size()) + " while " + 
						string("v.size=") + to_string(v.size())); 
			}
		}
		else if (var == "w") 
		{ 
			Array3D<double> tmp;
			while (OneLine >> value) tmp.push_back(value);
			if (tmp.size() == w.size()) { w = tmp; }
			else if (tmp.size() == 1) { w.setConstant(tmp[0]); }
			else
			{ 
				throw invalid_argument(string("Size of input w is wrong! ") +
						string("tmp.size=") + to_string(tmp.size()) + " while " + 
						string("w.size=") + to_string(w.size())); 
			}
		}
		else if (var == "p") 
		{ 
			Array3D<double> tmp;
			while (OneLine >> value) tmp.push_back(value);
			if (tmp.size() == p.size()) { p = tmp; }
			else if (tmp.size() == 1) { p.setConstant(tmp[0]); }
			else
			{ 
				throw invalid_argument(string("Size of input p is wrong! ") +
						string("tmp.size=") + to_string(tmp.size()) + " while " + 
						string("p.size=") + to_string(w.size())); 
			}
		}
	}

	file.close();

}


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
	p.initShape(-1, Nx, -1, Ny, -1, Nz);

	u.setZeros(); v.setZeros(); w.setZeros(); p.setZeros();


	/*********************************************************
	 * Setting the value on the boundary for Dirichlet BCs.
	 * Note: pressure is always Neumann BC in this solver.
	 *********************************************************/
	function<void(CI &, CI &)> f;

	// -x direction
	f = [&, this] (CI & i, CI & j) 
	{ u(0, i, j) = BCs[-1].get_uBCvalue();} ;

	if (BCs[-1].get_uType() == 1) dualLoop(0, Nyu, 0, Nzu, f);

	// -y direction
	f = [&, this] (CI & i, CI & j) 
	{ v(i, 0, i) = BCs[-2].get_vBCvalue();} ;

	if (BCs[-2].get_vType() == 1) dualLoop(0, Nxv, 0, Nzv, f);
	
	// -z direction
	f = [&, this] (CI & i, CI & j) 
	{ w(i, j, 0) = BCs[-3].get_wBCvalue();} ;

	if (BCs[-3].get_wType() == 1) dualLoop(0, Nxw, 0, Nyw, f);

	// +x direction
	f = [&, this] (CI & i, CI & j) 
	{ u(Nxu-1, i, j) = BCs[1].get_uBCvalue();} ;

	if (BCs[1].get_uType() == 1) dualLoop(0, Nyu, 0, Nzu, f);

	// +y direction
	f = [&, this] (CI & i, CI & j) 
	{ v(i, Nyv-1, i) = BCs[2].get_vBCvalue();} ;

	if (BCs[2].get_vType() == 1) dualLoop(0, Nxv, 0, Nzv, f);
	
	// +z direction
	f = [&, this] (CI & i, CI & j) 
	{ w(i, j, Nzw-1) = BCs[3].get_wBCvalue();} ;

	if (BCs[3].get_wType() == 1) dualLoop(0, Nxw, 0, Nyw, f);

	return 0;
}

