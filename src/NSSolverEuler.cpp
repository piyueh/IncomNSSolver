# include "include/IncomNSSolver.h"
# include "NS_updateGhost.cpp"
# include "NS_Convection.cpp"
# include "NS_Diffusion.cpp"


int NSSolverEuler::solve(int targetNStep)
{

	for(int n=0; n<targetNStep; ++n)
	{
		cout << "n=" << n << " ";

		updateGhost();

		PredictStep();

		updatePoissonSource();

		pSolver.Solve(b, p);

		updateU();

		time += dt;

		cout << "time = " << time << endl;
	}

	return 0;
}




int NSSolverEuler::PredictStep()
{

	for(int i=0; i<Nxu; ++i){
		for(int j=0; j<Nyu; ++j){
			for(int k=0; k<Nzu; ++k){
				u_str(i, j, k) = u(i, j, k) + 
					dt * (DiffusiveU(i, j, k) - ConvectU(i, j, k));
			}
		}
	}

	for(int i=0; i<Nxv; ++i){
		for(int j=0; j<Nyv; ++j){
			for(int k=0; k<Nzv; ++k){
				v_str(i, j, k) = v(i, j, k) + 
					dt * (DiffusiveV(i, j, k) - ConvectV(i, j, k));
			}
		}
	}

	for(int i=0; i<Nxw; ++i){
		for(int j=0; j<Nyw; ++j){
			for(int k=0; k<Nzw; ++k){
				w_str(i, j, k) = w(i, j, k) + 
					dt * (DiffusiveW(i, j, k) - ConvectW(i, j, k));
			}
		}
	}

	return 0;
}


int NSSolverEuler::updatePoissonSource()
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


int NSSolverEuler::updateU()
{
	int tmp1, tmp2;
	
	// update u
	for(int j=0; j<Nyu; ++j){
		for(int k=0; k<Nzu; ++k){
			u(0, j, k) = u_str(0, j, k);
			u(Nxu-1, j, k) = u_str(Nxu-1, j, k);
		}
	}

	for(int i=1; i<Nxu-1; ++i){
		for(int j=0; j<Nyu; ++j){
			for(int k=0; k<Nzu; ++k){
				tmp1 = i * Nyz + j * Nz + k; tmp2 = tmp1 - Nyz;
				u(i, j, k) = u_str(i, j, k) -
					dt * (p(tmp1) - p(tmp2)) / dx;	
			}
		}
	}

	// update v
	for(int i=0; i<Nxv; ++i){
		for(int k=0; k<Nzv; ++k){
			v(i, 0, k) = v_str(i, 0, k);
			v(i, Nyv-1, k) = u_str(i, Nyv-1, k);
		}
	}

	for(int i=0; i<Nxv; ++i){
		for(int j=1; j<Nyv-1; ++j){
			for(int k=0; k<Nzv; ++k){
				tmp1 = i * Nyz + j * Nz + k; tmp2 = tmp1 - Nyz;
				v(i, j, k) = v_str(i, j, k) -
					dt * (p(tmp1) - p(tmp2)) / dy;	
			}
		}
	}

	//update w
	for(int i=0; i<Nxw; ++i){
		for(int j=0; j<Nyw; ++j){
			w(i, j, 0) = w_str(i, j, 0);
			v(i, j, Nzw-1) = u_str(i, j, Nzw-1);
		}
	}

	for(int i=0; i<Nxw; ++i){
		for(int j=0; j<Nyw; ++j){
			for(int k=1; k<Nzw-1; ++k){
				tmp1 = i * Nyz + j * Nz + k; tmp2 = tmp1 - Nyz;
				w(i, j, k) = w_str(i, j, k) -
					dt * (p(tmp1) - p(tmp2)) / dz;	
			}
		}
	}
	   	
	return 0;
}




