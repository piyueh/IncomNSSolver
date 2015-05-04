# include "include/IncomNSSolver.h"

Data::Data(const string & fName, Mesh & mesh)
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
			N[1][0] = N1 - 2; N[1][1] = N2 - 2; N[1][2] = N3 - 2;
			u.initShape(-1, N[1][0], -1, N[1][1], -1, N[1][2]);
		}
		else if (var == "Nv") 
		{
			OneLine >> N1 >> N2 >> N3;
			N[2][0] = N1 - 2; N[2][1] = N2 - 2; N[2][2] = N3 - 2;
			v.initShape(-1, N[2][0], -1, N[2][1], -1, N[2][2]);
		}
		else if (var == "Nw") 
		{
			OneLine >> N1 >> N2 >> N3;
			N[3][0] = N1 - 2; N[3][1] = N2 - 2; N[3][2] = N3 - 2;
			w.initShape(-1, N[3][0], -1, N[3][1], -1, N[3][2]);
		}
		else if (var == "Np") 
		{
			OneLine >> N1 >> N2 >> N3;
			N[0][0] = N1 - 2; N[0][1] = N2 - 2; N[0][2] = N3 - 2;
			p.initShape(-1, N[0][0], -1, N[0][1], -1, N[0][2]);
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


	if (N[1] != array<int, 3>({mesh.get_Nxu(), mesh.get_Nyu(), mesh.get_Nzu()}))
	{
		throw invalid_argument(string("The mesh setting and ") +
				"the initial value data differ!! => u");
	}

	if (N[2] != array<int, 3>({mesh.get_Nxv(), mesh.get_Nyv(), mesh.get_Nzv()}))
	{
		throw invalid_argument(string("The mesh setting and ") +
				"the initial value data differ!! => v");
	}
	
	if (N[3] != array<int, 3>({mesh.get_Nxw(), mesh.get_Nyw(), mesh.get_Nzw()}))
	{
		cout << N[3][0] << " " << N[3][1] << " " << N[3][2] << endl;
		cout << mesh.get_Nxw() << " " << mesh.get_Nyw() << " " << mesh.get_Nzw() << endl;
		throw invalid_argument(string("The mesh setting and ") +
				"the initial value data differ!! => w");
	}
	
	if (N[0] != array<int, 3>({mesh.get_Nx(), mesh.get_Ny(), mesh.get_Nz()}))
	{
		throw invalid_argument(string("The mesh setting and ") +
				"the initial value data differ!! => p");
	}

}


int Data::InitData(Mesh & mesh)
{

	map<int, Boundary> & BCs = mesh.get_BCs();

	time = 0;

	N[0][0] = mesh.get_Nx(); N[0][1] = mesh.get_Ny(); N[0][2] = mesh.get_Nz();
	N[1][0] = mesh.get_Nxu(), N[1][1] = mesh.get_Nyu(), N[1][2] = mesh.get_Nzu();
	N[2][0] = mesh.get_Nxv(), N[2][1] = mesh.get_Nyv(), N[2][2] = mesh.get_Nzv();
	N[3][0] = mesh.get_Nxw(), N[3][1] = mesh.get_Nyw(), N[3][2] = mesh.get_Nzw();

	p.initShape(-1, N[0][0], -1, N[0][1], -1, N[0][2]);
	u.initShape(-1, N[1][0], -1, N[1][1], -1, N[1][2]);
	v.initShape(-1, N[2][0], -1, N[2][1], -1, N[2][2]);
	w.initShape(-1, N[3][0], -1, N[3][1], -1, N[3][2]);

	p.setZeros(); u.setZeros(); v.setZeros(); w.setZeros(); 

	return 0;
}


