int tripleLoop(int &ia, int &ib, int &ja, int &jb, int &ka, int &kb,  
		function<void(int &, int &, int &)> f)
{
	for(int i=ia; i<ib, ++i)
	{
		for(int j=ja; j<jb; ++j)
		{
			for(int k=ka; k<kb; ++k)
			   
				f(i, j, k);
		}
	}
}


int dualLoop(int &ia, int &ib, int &ja, int &jb, 
		function<void(int &, int &)> f)
{
	for(int i=ia; i<ib, ++i)
	{
		for(int j=ja; j<jb; ++j)

			f(i, j);
	}
}
