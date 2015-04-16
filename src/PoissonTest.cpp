# include "include/IncomNSSolver.h"
# include "testFuncs.cpp"


double evalRelErr(VectorXd & x, VectorXd & xe)
{
	auto tmp = x - xe;
	return tmp.cwiseAbs().maxCoeff();
}


int main()
{

	
	return 0;
}


