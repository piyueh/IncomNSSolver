# include <iostream>
# include <string>
# include <vector>
# include "class_Boundary.h"


Boundary::Boundary(int N, string dir, vector<int> idx, int t, double v)
{
	Ncells = N;
	direction = dir;
	Cell = idx;
	Type = t;
	value = v;
}


void Boundary::print()
{
	cout << "# of cells in this boundary: " << Ncells << endl;
	cout << "The normal direction: " << direction << endl;
	cout << "The type of BC: " << Type << endl;
	cout << "The value of BC: " << value << endl;

	cout << "The index of cells in this boundary (in global index): " << endl;
	for(auto it=Cell.cbegin(); it<Cell.cend(); ++it)
		cout << *it << " ";
	cout << endl;
}
