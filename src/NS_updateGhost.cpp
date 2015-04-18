int NSSolverEuler::updateGhost()
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
