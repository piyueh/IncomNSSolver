
NSSolverEuler::solve()
{
	// predict step
	predict_U();
	predict_V();
	predict_W();	
	// set LHS of Poisson solver
	
	// solve pressure
	pSolver.solve(p);

	// update velocity
	
}


NSSolverEuler::predict_u()
{
	// update the ghost cells
	update_ghostU();

	// add convective term to u_star
	convection_U();

	// add diffusion term to u_star
	diffusion_U();
}


NSSolverEuler::update_ghostU();
{
	u_str.setZeros();

	for(int j=0; j<Nuy; ++j)
	{
		for(int k=0; k<Nuz; ++k)
		{
			u_str(-1, j, k) = u(Nux-2, j, k);
			u_str(Nux, j, k) = u(1, j, k);
		}
	}
			
	for(int i=0; i<Nux; ++i)
	{
		for(int k=0; k<Nuz; ++k)
		{
			u_str(i, -1, k) = u(i, Nuy-1, k);
			u_str(i, Nuy, k) = u(i, 0, k);
		}
	}
			
	for(int i=0; i<Nux; ++i)
	{
		for(int j=0; j<Nuy; ++j)
		{
			u_str(i, j, -1) = u(i, j, Nuz-1);
			u_str(i, j, Nuz) = u(i, j, 0);
		}
	}
}


NSSolverEuler::update_ghostV();
{
	v_str.setZeros();

	for(int j=0; j<Nvy; ++j)
	{
		for(int k=0; k<Nvz; ++k)
		{
			v_str(-1, j, k) = v(Nvx-1, j, k);
			v_str(Nvx, j, k) = v(0, j, k);
		}
	}
			
	for(int i=0; i<Nvx; ++i)
	{
		for(int k=0; k<Nvz; ++k)
		{
			v_str(i, -1, k) = v(i, Nuy-1, k);
			v_str(i, Nuy, k) = v(i, 0, k);
		}
	}
			
	for(int i=0; i<Nux; ++i)
	{
		for(int j=0; j<Nuy; ++j)
		{
			u_str(i, j, -1) = u(i, j, Nuz-1);
			u_str(i, j, Nuz) = u(i, j, 0);
		}
	}
}
