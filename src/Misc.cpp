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
