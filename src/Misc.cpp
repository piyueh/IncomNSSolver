# include "include/IncomNSSolver.h"

int tripleLoop(const int &ia, const int &ib, const int &ja, const int &jb, 
		const int &ka, const int &kb, function<void(int &, int &, int &)> f)
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


int dualLoop(const int &ia, const int &ib, const int &ja, const int &jb, 
		function<void(int &, int &)> f)
{
	for(int i=ia; i<ib; ++i)
	{
		for(int j=ja; j<jb; ++j)

			f(i, j);
	}

	return 0;
}
