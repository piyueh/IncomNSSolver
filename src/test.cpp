# include <iostream>
# include <string>
# include <vector>
# include "include/IncomNSSolver.h"

using namespace std;

int main()
{
	int N = 9;
	vector<int> idxs;
	string dir = "+x";
	int type = 1;
	double value=0.;

	for(int j=0; j<3; ++j)
	{
		for(int k=0; k<3; ++k)
			idxs.push_back(j*3+k);
	}

	Boundary test;
	test = {N, dir, idxs, type, value};
	test.print();
	
	return 0;
}
