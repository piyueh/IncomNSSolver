# include <iostream>
# include <string>
# include <vector>
# include <array>
# include <map>
using namespace std;

# include "include/IncomNSSolver.h"
# include "include/io.h"



int main()
{
	map<string, Boundary> BCs;

	BCs["+x"] = Boundary(3, 3, 1, "+x", 1, 0, 1, 0) ;	
	BCs["-x"] = Boundary(3, 3, 1, "-x", 1, 0, 1, 0) ;	
	BCs["+y"] = Boundary(3, 3, 1, "+y", 1, 0, 1, 0) ;	
	BCs["-y"] = Boundary(3, 3, 1, "-y", 1, 0, 1, 0) ;	
	BCs["+z"] = Boundary(3, 3, 1, "+z", 1, 0, 1, 0) ;	
	BCs["-z"] = Boundary(3, 3, 1, "-z", 1, 0, 1, 0) ;	

	Array3D<double> A;

	A.initShape(-1, 3, 0, 2, 0, 2);

	int i=0;
	for(auto it=A.begin(); it!=A.end(); ++it)
	{
		*it = i;
		++i;
	}

	for(auto &i: A)
		cout << i << ", ";
	cout << endl << endl;

	for(i=-1; i<4; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			for(int k=0; k<3; ++k)
				cout << A(i, j, k) << ", ";
		}
	}
	cout << endl;
	cout << A << endl;
	cout << A.shape() << endl;
	/*
	for(auto it=BCs.begin(); it!=BCs.end(); ++it)
		cout << it->second << endl;
		*/


	/*
	int A[3][3];

	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j)
			A[i][j] = i * 10 + j;
	}

	for(auto iti=begin(A); iti<end(A); ++iti){
		for(auto itj=begin(*iti); itj<end(*iti); ++itj)
			cout << *itj << endl;
	}

	cout << A << endl;
	cout << begin(A) << endl;
	cout << &A[2][2] << endl;
	cout << end(A) << endl;
	*/

	return 0;
}
