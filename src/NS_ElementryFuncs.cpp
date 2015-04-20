int NSSolver::InitLambda()
{
	calIdx = [this] (CI & i, CI & j, CI & k) -> int 
	{ 
		return i * Nyz + j * Nz + k; 
	};
	

	/*
	updUB = [this] (CI & i, CI & j) -> void 
	{ 
		u(0, i, j) = u_str(0, i, j); 
		u(Nxu-1, i, j) = u_str(Nxu-1, i, j); 
	};

	updVB = [this] (CI & i, CI & j) -> void 
	{ 
		v(i, 0, j) = v_str(i, 0, j); 
		v(i, Nyv-1, j) = v_str(i, Nyv-1, j); 
	};

	updWB = [this] (CI & i, CI & j) -> void 
	{ 
		w(i, j, 0) = w_str(i, j, 0); 
		w(i, j, Nzw-1) = w_str(i, j, Nzw-1); 
	}; 
	*/

	updU = [this] (CI & i, CI & j, CI & k, CD & DT) -> void 
	{ 
		int tmp1 = calIdx(i, j, k); int tmp2 = tmp1 - Nyz; 
		u(i, j, k) = u(i, j, k) - DT * (p(tmp1) - p(tmp2)) / dx; 
	};

	updV = [this] (CI & i, CI & j, CI & k, CD & DT) -> void 
	{ 
		int tmp1 = calIdx(i, j, k); int tmp2 = tmp1 - Nz; 
		v(i, j, k) = v(i, j, k) - DT * (p(tmp1) - p(tmp2)) / dy; 
	};

	updW = [this] (CI & i, CI & j, CI & k, CD & DT) -> void 
	{ 
		int tmp1 = calIdx(i, j, k); int tmp2 = tmp1 - 1; 
		w(i, j, k) = w(i, j, k) - DT * (p(tmp1) - p(tmp2)) / dz; 
	};


	DivOnPresPt = [this] (CI & i, CI & j, CI & k) -> void
	{
		int tmp = calIdx(i, j, k);
		b(tmp) = 
			(u(i+1, j, k) - u(i, j, k)) / dx + 
			(v(i, j+1, k) - v(i, j, k)) / dy + 
			(w(i, j, k+1) - w(i, j, k)) / dz;

	};


	ConvU = [this] (CI & i, CI & j, CI & k) -> double
	{
		return
			0.25 * ((u(i+1, j, k) + u(i, j, k)) * (u(i+1, j, k) + u(i, j, k)) -
					(u(i, j, k) + u(i-1, j, k)) * (u(i, j, k) + u(i-1, j, k))) / dx + 

			0.25 * ((u(i, j+1, k) + u(i, j, k)) * (v(i, j+1, k) + v(i-1, j+1, k)) -
					(u(i, j, k) + u(i, j-1, k)) * (v(i, j, k) + v(i-1, j, k))) / dy +

			0.25 * ((u(i, j, k+1) + u(i, j, k)) * (w(i, j, k+1) + w(i-1, j, k+1)) -
					(u(i, j, k) + u(i, j, k-1)) * (w(i, j, k) + w(i-1, j, k))) / dz; 
	};

	ConvV = [this] (CI & i, CI & j, CI & k) -> double
	{
		return
			0.25 * ((v(i+1, j, k) + v(i, j, k)) * (u(i+1, j, k) + u(i+1, j-1, k)) -
					(v(i, j, k) + v(i-1, j, k)) * (u(i, j, k) + u(i, j-1, k))) / dx + 

			0.25 * ((v(i, j+1, k) + v(i, j, k)) * (v(i, j+1, k) + v(i, j, k)) -
					(v(i, j, k) + v(i, j-1, k)) * (v(i, j, k) + v(i, j-1, k))) / dy +

			0.25 * ((v(i, j, k+1) + v(i, j, k)) * (w(i, j, k+1) + w(i, j-1, k+1)) -
					(v(i, j, k) + v(i, j, k-1)) * (w(i, j, k) + w(i, j-1, k))) / dz; 
	};


	ConvW = [this] (CI & i, CI & j, CI & k) -> double
	{
		return
			0.25 * ((w(i+1, j, k) + v(i, j, k)) * (u(i+1, j, k) + u(i+1, j, k-1)) -
					(w(i, j, k) + v(i-1, j, k)) * (u(i, j, k) + u(i, j, k-1))) / dx +

			0.25 * ((w(i, j+1, k) + w(i, j, k)) * (v(i, j+1, k) + v(i, j+1, k-1)) -
					(w(i, j, k) + w(i, j-1, k)) * (v(i, j, k) + v(i, j, k-1))) / dy + 

			0.25 * ((w(i, j, k+1) + w(i, j, k)) * (w(i, j, k+1) + w(i, j, k)) -
					(w(i, j, k) + w(i, j, k-1)) * (w(i, j, k) + w(i, j, k-1))) / dz;

	};


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


	updGu = [this] (CI & i, CI & j, CI & k) -> void
	{
		Gu(i, j, k) = DiffU(i, j, k) - ConvU(i, j, k);
	};


	updGv = [this] (CI & i, CI & j, CI & k) -> void
	{
		Gv(i, j, k) = DiffV(i, j, k) - ConvV(i, j, k);
	};


	updGw = [this] (CI & i, CI & j, CI & k) -> void
	{
		Gw(i, j, k) = DiffW(i, j, k) - ConvW(i, j, k);
	};


	updGu2 = [this] (CI & i, CI & j, CI & k, CD & coef) -> void
	{
		Gu(i, j, k) = coef * Gu(i, j, k) + DiffU(i, j, k) - ConvU(i, j, k);
	};


	updGv2 = [this] (CI & i, CI & j, CI & k, CD & coef) -> void
	{
		Gv(i, j, k) = coef * Gv(i, j, k) +  DiffV(i, j, k) - ConvV(i, j, k);
	};


	updGw2 = [this] (CI & i, CI & j, CI & k, CD & coef) -> void
	{
		Gw(i, j, k) = coef * Gw(i, j, k) +  DiffW(i, j, k) - ConvW(i, j, k);
	};


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

