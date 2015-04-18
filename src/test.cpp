# include "include/IncomNSSolver.h"
# include "Misc.cpp"
# include "io.cpp"

int main()
{
	array<int ,3> N, Nu, Nv, Nw;

	map<int, Boundary> BCs;

	N = {3, 3, 1};
	Nu = N; Nv = N; Nw = N;
	Nu[0] += 1; Nv[1] += 1; Nw[2] += 1;

	BCs[1] = Boundary(N, 1, 1, {-1, 0}, {1, 1}, {0, 0}, {0, 0});
	BCs[2] = Boundary(N, 2, 1, {-1, 0}, {1, 1}, {0, 0}, {0, 0});
	BCs[3] = Boundary(N, 3, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});
	BCs[-1] = Boundary(N, 1, -1, {-1, 0}, {1, 1}, {0, 0}, {0, 0});
	BCs[-2] = Boundary(N, 2, -1, {-1, 0}, {1, 1}, {0, 0}, {0, 0});
	BCs[-3] = Boundary(N, 3, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});

	Array3D<double> u, v, w, p;

	u.initShape(-1, Nu[0], -1, Nu[1], -1, Nu[2]);
	v.initShape(-1, Nv[0], -1, Nv[1], -1, Nv[2]);
	w.initShape(-1, Nw[0], -1, Nw[1], -1, Nw[2]);
	p.initShape(-1, N[0], -1, N[0], -1, N[0]);

	u.setZeros(); v.setZeros(); w.setZeros(); p.setZeros(); 

	for(int i=0; i<Nu[0]; ++i){
		for(int j=0; j<Nu[1]; ++j){
			for(int k=0; k<Nu[2]; ++k)
				u(i, j, k) = (i+1) * 10 + j+1;
		}
	}

	for(int i=0; i<Nv[0]; ++i){
		for(int j=0; j<Nv[1]; ++j){
			for(int k=0; k<Nv[2]; ++k)
				v(i, j, k) = (i+1) * 10 + j+1;
		}
	}

	for(int i=0; i<Nw[0]; ++i){
		for(int j=0; j<Nw[1]; ++j){
			for(int k=0; k<Nw[2]; ++k)
				w(i, j, k) = (i+1) * 10 + j+1;
		}
	}

	ofstream file("Data0.txt");
	file << u << endl;	
	file << v << endl;	
	file << w << endl;	

	BCs[1].updGhost(Nu[1], Nu[2], 1, u, 1./N[0]); 
	BCs[2].updGhost(Nu[0], Nu[2], 1, u, 1./N[1]); 
	BCs[3].updGhost(Nu[0], Nu[1], 1, u, 1./N[2]); 
	BCs[-1].updGhost(Nu[1], Nu[2], 1, u, 1./N[0]); 
	BCs[-2].updGhost(Nu[0], Nu[2], 1, u, 1./N[1]); 
	BCs[-3].updGhost(Nu[0], Nu[1], 1, u, 1./N[2]); 

	BCs[1].updGhost(Nv[1], Nv[2], 2, v, 1./N[0]); 
	BCs[2].updGhost(Nv[0], Nv[2], 2, v, 1./N[1]); 
	BCs[3].updGhost(Nv[0], Nv[1], 2, v, 1./N[2]); 
	BCs[-1].updGhost(Nv[1], Nv[2], 2, v, 1./N[0]); 
	BCs[-2].updGhost(Nv[0], Nv[2], 2, v, 1./N[1]); 
	BCs[-3].updGhost(Nv[0], Nv[1], 2, v, 1./N[2]); 

	ofstream file2("Data1.txt");
	file2 << u << endl;	
	file2 << v << endl;	
	file2 << w << endl;	

	return 0;
}

