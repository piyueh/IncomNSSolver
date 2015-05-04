int NSSolver::InitLambda()
{
	calIdx = [this] (CI & i, CI & j, CI & k) -> int 
	{ 
		return i * Nyz + j * Nz + k; 
	};
	

	/*****************************************************************
	 * Updating the field
	 *****************************************************************/
	updU = [this] (CI & i, CI & j, CI & k) -> void 
	{ 
		u(i, j, k) -= (dp(i, j, k) - dp(i-1, j, k)) / dx; 
	};

	updV = [this] (CI & i, CI & j, CI & k) -> void 
	{ 
		v(i, j, k) -= (dp(i, j, k) - dp(i, j-1, k)) / dy; 
	};

	updW = [this] (CI & i, CI & j, CI & k) -> void 
	{ 
		w(i, j, k) -= (dp(i, j, k) - dp(i, j, k-1)) / dz; 
	};

	updP = [this] (CI & i, CI & j, CI & k, CD & coeff) -> void 
	{ 
		p(i, j, k) += dp(i, j, k) / coeff; 
	};


	/*****************************************************************
	 * 
	 *****************************************************************/
	DivOnPresPt = [this] (CI & i, CI & j, CI & k) -> void
	{
		b(i, j, k) = 
			(u(i+1, j, k) - u(i, j, k)) / dx + 
			(v(i, j+1, k) - v(i, j, k)) / dy + 
			(w(i, j, k+1) - w(i, j, k)) / dz;

	};


	/*****************************************************************
	 * Convective terms start from here
	 *****************************************************************/
	ConvU = [this] (CI & i, CI & j, CI & k) -> double
	{
		return
			- 0.25 * ((u(i+1, j, k) + u(i, j, k)) * (u(i+1, j, k) + u(i, j, k)) -
					  (u(i, j, k) + u(i-1, j, k)) * (u(i, j, k) + u(i-1, j, k))) / dx - 

			  0.25 * ((u(i, j+1, k) + u(i, j, k)) * (v(i, j+1, k) + v(i-1, j+1, k)) -
					  (u(i, j, k) + u(i, j-1, k)) * (v(i, j, k) + v(i-1, j, k))) / dy -

			  0.25 * ((u(i, j, k+1) + u(i, j, k)) * (w(i, j, k+1) + w(i-1, j, k+1)) -
					  (u(i, j, k) + u(i, j, k-1)) * (w(i, j, k) + w(i-1, j, k))) / dz; 
	};

	ConvV = [this] (CI & i, CI & j, CI & k) -> double
	{
		return
			- 0.25 * ((v(i+1, j, k) + v(i, j, k)) * (u(i+1, j, k) + u(i+1, j-1, k)) -
					  (v(i, j, k) + v(i-1, j, k)) * (u(i, j, k) + u(i, j-1, k))) / dx - 

			  0.25 * ((v(i, j+1, k) + v(i, j, k)) * (v(i, j+1, k) + v(i, j, k)) -
					  (v(i, j, k) + v(i, j-1, k)) * (v(i, j, k) + v(i, j-1, k))) / dy -

			  0.25 * ((v(i, j, k+1) + v(i, j, k)) * (w(i, j, k+1) + w(i, j-1, k+1)) -
					  (v(i, j, k) + v(i, j, k-1)) * (w(i, j, k) + w(i, j-1, k))) / dz; 
	};


	ConvW = [this] (CI & i, CI & j, CI & k) -> double
	{
		return
			- 0.25 * ((w(i+1, j, k) + w(i, j, k)) * (u(i+1, j, k) + u(i+1, j, k-1)) -
					  (w(i, j, k) + w(i-1, j, k)) * (u(i, j, k) + u(i, j, k-1))) / dx -

			  0.25 * ((w(i, j+1, k) + w(i, j, k)) * (v(i, j+1, k) + v(i, j+1, k-1)) -
					  (w(i, j, k) + w(i, j-1, k)) * (v(i, j, k) + v(i, j, k-1))) / dy - 

			  0.25 * ((w(i, j, k+1) + w(i, j, k)) * (w(i, j, k+1) + w(i, j, k)) -
					  (w(i, j, k) + w(i, j, k-1)) * (w(i, j, k) + w(i, j, k-1))) / dz;

	};


	/*****************************************************************
	 * Diffusive terms start from here
	 *****************************************************************/
	DiffU = [this] (CI & i, CI & j, CI & k) -> double
	{
		return nu * (
				(u(i+1, j, k) - 2 * u(i, j, k) + u(i-1, j, k)) / dx2 +
				(u(i, j+1, k) - 2 * u(i, j, k) + u(i, j-1, k)) / dy2 +
				(u(i, j, k+1) - 2 * u(i, j, k) + u(i, j, k-1)) / dz2 );
	};


	DiffV = [this] (CI & i, CI & j, CI & k) -> double
	{
		return nu * (
				(v(i+1, j, k) - 2 * v(i, j, k) + v(i-1, j, k)) / dx2 +
				(v(i, j+1, k) - 2 * v(i, j, k) + v(i, j-1, k)) / dy2 +
				(v(i, j, k+1) - 2 * v(i, j, k) + v(i, j, k-1)) / dz2 );
	};


	DiffW = [this] (CI & i, CI & j, CI & k) -> double
	{
		return nu * (
				(w(i+1, j, k) - 2 * w(i, j, k) + w(i-1, j, k)) / dx2 +
				(w(i, j+1, k) - 2 * w(i, j, k) + w(i, j-1, k)) / dy2 +
				(w(i, j, k+1) - 2 * w(i, j, k) + w(i, j, k-1)) / dz2 );
	};


	/*****************************************************************
	 * Pressure terms start from here
	 *****************************************************************/
	PresU = [this] ( CI & i, CI & j, CI & k) -> double
	{
		return - (p(i, j, k) - p(i-1, j, k)) / dx; 
	};

	PresV = [this] ( CI & i, CI & j, CI & k) -> double
	{
		return - (p(i, j, k) - p(i, j-1, k)) / dy; 
	};
		
	PresW = [this] ( CI & i, CI & j, CI & k) -> double
	{
		return - (p(i, j, k) - p(i, j, k-1)) / dz; 
	};


	/*****************************************************************
	 * Functions to calculate Gu, Gv, Gw in the interior area
	 *****************************************************************/
	updGu = [this] (CI & i, CI & j, CI & k) -> void
	{
		Gu(i, j, k) = DiffU(i, j, k) + ConvU(i, j, k) + PresU(i, j, k);
	};


	updGv = [this] (CI & i, CI & j, CI & k) -> void
	{
		Gv(i, j, k) = DiffV(i, j, k) + ConvV(i, j, k) + PresV(i, j, k);
	};


	updGw = [this] (CI & i, CI & j, CI & k) -> void
	{
		Gw(i, j, k) = DiffW(i, j, k) + ConvW(i, j, k) + PresW(i, j, k);
	};


	/*****************************************************************
	 * Overloading version of the functions calculating Gu, Gv, Gw
	 *****************************************************************/
	updGu2 = [this] (CI & i, CI & j, CI & k, CD & coef) -> void
	{
		Gu(i, j, k) = coef * Gu(i, j, k) + 
			DiffU(i, j, k) + ConvU(i, j, k) + PresU(i, j, k);
	};


	updGv2 = [this] (CI & i, CI & j, CI & k, CD & coef) -> void
	{
		Gv(i, j, k) = coef * Gv(i, j, k) + 
			DiffV(i, j, k) + ConvV(i, j, k) + PresV(i, j, k);
	};


	updGw2 = [this] (CI & i, CI & j, CI & k, CD & coef) -> void
	{
		Gw(i, j, k) = coef * Gw(i, j, k) + 
			DiffW(i, j, k) + ConvW(i, j, k) + PresW(i, j, k);
	};


	/*****************************************************************
	 * Updating velocity in prediction step
	 *****************************************************************/
	preU = [this] (CI & i, CI & j, CI & k, CD & DT) -> void
	{
		u(i, j, k) = u(i, j, k) + DT * Gu(i, j, k);
	};


	preV = [this] (CI & i, CI & j, CI & k, CD & DT) -> void
	{
		v(i, j, k) = v(i, j, k) + DT * Gv(i, j, k);
	};


	preW = [this] (CI & i, CI & j, CI & k, CD & DT) -> void
	{
		w(i, j, k) = w(i, j, k) + DT * Gw(i, j, k);
	};

	return 0;
}

