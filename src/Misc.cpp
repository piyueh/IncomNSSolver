# include "include/IncomNSSolver.h"

int tripleLoop(CI & ia, CI & ib, CI & ja, CI & jb, CI & ka, CI & kb, 
		function<void(CI &, CI &, CI &)> f)
{
	for(int i=ia; i<ib; ++i)
	{
		for(int j=ja; j<jb; ++j)
		{
			for(int k=ka; k<kb; ++k)
			   
				f(i, j, k);
		}
	}

	return 0;
}


int tripleLoop(CI & ia, CI & ib, CI & ja, CI & jb, CI & ka, CI & kb, 
		CD & coef, function<void(CI &, CI &, CI &, CD &)> f)
{
	for(int i=ia; i<ib; ++i)
	{
		for(int j=ja; j<jb; ++j)
		{
			for(int k=ka; k<kb; ++k)
			   
				f(i, j, k, coef);
		}
	}

	return 0;
}


int dualLoop(CI & ia, CI & ib, CI & ja, CI & jb,
	   	function<void(CI &, CI &)> f)
{
	for(int i=ia; i<ib; ++i)
	{
		for(int j=ja; j<jb; ++j)

			f(i, j);
	}

	return 0;
}


int TGVortex(CI & Nxu, CI & Nyu, CI & Nzu, A3Dd & u, VD & xu, VD & yu,
		CI & Nxv, CI & Nyv, CI & Nzv, A3Dd & v, VD & xv, VD & yv, CD & t)
{
	auto u_ex = [&] (CI & i, CI & j, CI & k) -> void
	{
		u(i, j, k) = - exp(- 2 * t) * cos(xu[i]) * sin(yu[j]);
	};

	auto v_ex = [&] (CI & i, CI & j, CI & k) -> void
	{
		v(i, j, k) = exp(- 2 * t) * sin(xv[i]) * cos(yv[j]);
	};

	tripleLoop(0, Nxu, 0, Nyu, 0, Nzu, u_ex);
	tripleLoop(0, Nxv, 0, Nyv, 0, Nzv, v_ex);

	return 0;
}
