

int updXGhost(Array3D<double> &u, const int &Ny, const int &Nz, 
		const int &type, const int &ghIdx, const int &corIdx, 
		const double coeff, const double &BCvalue)
{
	switch (type)
	{
		case 0:
			for(int j=0; j<Ny; ++j){
				for(int k=0; k<Nz; ++k){
					u(ghIdx, j, k) = u(corIdx, j, k); } }

			break;

		case 1: case -1:
			for(int j=0; j<Ny; ++j){
				for(int k=0; k<Nz; ++k){
					u(ghIdx, j, k) = 
						- type * u(corIdx, j, k) + coeff * BCvalue; } }

			break;

		default:
			throw invalid_argument("In EulerNSSsolver.cpp->updXGhost");
			break;
	}
	return 0;
}


int updYGhost(Array3D<double> &u, const int &Nx, const int &Nz, 
		const int &type, const int &ghIdx, const int &corIdx, 
		const double coeff, const double &BCvalue)
{
	switch (type)
	{
		case 0:
			for(int i=0; i<Nx; ++i){
				for(int k=0; k<Nz; ++k){
					u(i, ghIdx, k) = u(i, corIdx, k); } }

			break;

		case 1: case -1:
			for(int i=0; i<Nx; ++i){
				for(int k=0; k<Nz; ++k){
					u(i, ghIdx, k) = 
						- type * u(i, corIdx, k) + coeff * BCvalue; } }

			break;

		default:
			throw invalid_argument("In EulerNSSsolver.cpp->updXGhost");
			break;
	}
	return 0;
}


int updZGhost(Array3D<double> &u, const int &Nx, const int &Ny, 
	const int &type, const int &ghIdx, const int &corIdx, 
	const double coeff, const double &BCvalue)
{
	switch (type)
	{
		case 0:
			for(int i=0; i<Nx; ++i){
				for(int j=0; j<Ny; ++j){
					u(i, j, ghIdx) = u(i, j, corIdx); } }

			break;

		case 1: case -1:
			for(int i=0; i<Nx; ++i){
				for(int j=0; j<Ny; ++j){
					u(i, j, ghIdx) = 
						- type * u(i, j, corIdx) + coeff * BCvalue; } }

			break;

		default:
			throw invalid_argument("In EulerNSSsolver.cpp->updXGhost");
			break;
	}
	return 0;
}


int NSSolverEuler::updateGhost()
{
	for(auto &bcPair: mesh.get_BCs())
	{
		auto & bc = bcPair.second;
		
		switch (bc.get_Dir()) {
			case 1: // x

				updXGhost(u, Nyu, Nzu, 
						bc.get_uType(), bc.get_uBCIdx(), bc.get_uBCcorIdx(), 
						2*dx, bc.get_uBCvalue());

				updXGhost(v, Nyv, Nzv, 
						bc.get_vType(), bc.get_vBCIdx(), bc.get_vBCcorIdx(), 
						dx, bc.get_vBCvalue());

				updXGhost(w, Nyw, Nzw, 
						bc.get_wType(), bc.get_wBCIdx(), bc.get_wBCcorIdx(), 
						dx, bc.get_wBCvalue());

				break;

			case 2: // y

				updYGhost(u, Nxu, Nzu, 
						bc.get_uType(), bc.get_uBCIdx(), bc.get_uBCcorIdx(), 
						dy, bc.get_uBCvalue());

				updYGhost(v, Nxv, Nzv, 
						bc.get_vType(), bc.get_vBCIdx(), bc.get_vBCcorIdx(), 
						2*dy, bc.get_vBCvalue());

				updYGhost(w, Nxw, Nzw, 
						bc.get_wType(), bc.get_wBCIdx(), bc.get_wBCcorIdx(), 
						dy, bc.get_wBCvalue());
				break;

			case 3: // y

				updZGhost(u, Nxu, Nyu, 
						bc.get_uType(), bc.get_uBCIdx(), bc.get_uBCcorIdx(), 
						dy, bc.get_uBCvalue());

				updZGhost(v, Nxv, Nyv, 
						bc.get_vType(), bc.get_vBCIdx(), bc.get_vBCcorIdx(), 
						dy, bc.get_vBCvalue());

				updZGhost(w, Nxw, Nyw, 
						bc.get_wType(), bc.get_wBCIdx(), bc.get_wBCcorIdx(), 
						2*dy, bc.get_wBCvalue());
				break;
		}
	}
	return 0;
}
