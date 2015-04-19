/*
 * This file defines overloaded iostream functionx.
 */


# include "include/IncomNSSolver.h"


template<typename T>
ostream & operator<<(ostream &os, vector<T> x)
{
	for(auto &i: x) os << i << " ";
	return os;
}



template<typename T>
ostream & operator<<(ostream &os, Array3D<T> &A)
{
	os << A.shape() << endl;
	
	for(auto &i: A) os << i << " ";

	return os;
}


ostream &operator<<(ostream &os, NSSolver &solver)
{
	os << "Current time: " << solver.time << endl;
	os << "Time step: " << solver.dt << endl;

	return os;
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


ostream &operator<<(ostream &os, Mesh &mesh)
{
	os << "Numbers of cells: "
	   << mesh.Nx << " x " << mesh.Ny << " x " << mesh.Nz << endl;
	os << "Shape of u array: "
	   << mesh.Nxu << " x " << mesh.Nyu << " x " << mesh.Nzu << endl;
	os << "Shape of v array: "
	   << mesh.Nxv << " x " << mesh.Nyv << " x " << mesh.Nzv << endl;
	os << "Shape of w array: "
	   << mesh.Nxw << " x " << mesh.Nyw << " x " << mesh.Nzw << endl;
	os << "Domain size: " 
	   << mesh.Lx << " x " << mesh.Ly << " x " << mesh.Lz << endl;
	os << "Cell size: " 
	   << mesh.dx << " x " << mesh.dy << " x " << mesh.dz << endl;
	os << endl;

	os << "Location of p point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xp << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yp << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zp << endl << endl;

	os << "Location of u point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xu << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yu << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zu << endl << endl;

	os << "Location of v point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xv << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yv << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zv << endl << endl;

	os << "Location of w point: " << endl;
	os << "\t x: " << endl;
	os << "\t\t" << mesh.xw << endl << endl;
	os << "\t y: " << endl;
	os << "\t\t" << mesh.yw << endl << endl;
	os << "\t z: " << endl;
	os << "\t\t" << mesh.zw << endl << endl;

	os << "Boundary conditions: " << endl;
	for(auto it=mesh.BCs.begin(); it!=mesh.BCs.end(); ++it)
		os << it->second << endl;

	return os;
}


int NSSolver::output(string fileName)
{
	ofstream file(fileName);

	file << u << endl;
	file << v << endl;
	file << w << endl;
	file << p << endl;

	return 0;
}

int NSSolver::output_u() { cout << u << endl; return 0; }
int NSSolver::output_v() { cout << v << endl; return 0; }
int NSSolver::output_w() { cout << w << endl; return 0; }
int NSSolver::output_p() { cout << p << endl; return 0; }
