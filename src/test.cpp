# include <string>
# include <vector>
# include <array>
# include <map>
# include <iostream>
# include <functional>
# include <utility>

using namespace std;

# include "include/class_Boundary.h"
# include "Boundary.cpp"

ostream &operator<<(ostream &os, Boundary &BC);

int main()
{
	Boundary BC;
	BC = Boundary({3, 3, 3}, 2, 1, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0});

	cout << BC << endl;

	return 0;
}
