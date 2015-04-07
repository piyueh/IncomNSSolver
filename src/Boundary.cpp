# include "include/IncomNSSolver.h"


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
	cout << "# of cells in this boundary: " << getNcells() << endl;
	cout << "The normal direction: " << direction << endl;
	cout << "The type of BC: " << getType() << endl;
	cout << "The value of BC: " << value << endl;

	cout << "The index of cells in this boundary (in global index): " << endl;
	for(auto it=Cell.cbegin(); it<Cell.cend(); ++it)
		cout << *it << " ";
	cout << endl;
}


vector<int>::const_iterator Boundary::bgCell()
{
	return Cell.cbegin();
}


vector<int>::const_iterator Boundary::edCell()
{
	return Cell.cend();
}


double Boundary::getBCvalue()
{
	return value;
}


int Boundary::getNcells()
{
	return Ncells;
}


string Boundary::getDirection()
{
	return direction;
}


int Boundary::getType()
{
	return Type;
}

