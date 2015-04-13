# include <iostream>
# include <string>
# include <vector>
# include <array>
# include <map>
using namespace std;

# include "include/class_Mesh.h"



int main()
{

	Mesh test;

	test.InitMesh({4, 4, 1}, {1, 1, 1});

	test.addBC(1, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	test.addBC(2, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	test.addBC(3, 1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	test.addBC(1, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	test.addBC(2, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	
	test.addBC(3, -1, {-1, 0}, {0, 0}, {0, 0}, {0, 0});	

	cout << test << endl;

	return 0;
}
