# include "include/IncomNSSolver.h"


Boundary::Boundary(int N, string dir, int t, double v, 
		vector<int> idx)
{
	Ncells = N;
	direction = dir;
	Type = t;
	value = v;
	Cell = idx;
}


Boundary::Boundary(int N, string dir, int t, double v, 
		vector<int> idx, vector<int> idx2)
{
	Ncells = N;
	direction = dir;
	Type = t;
	value = v;
	Cell = idx;
	OppCell = idx2;
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

	if (getType() == 0)
	{
		cout << "This is a periodic BC." 
			"The corresponding elements on the opposite surface are: " << endl;
		for(auto it=OppCell.cbegin(); it<OppCell.cend(); ++it)
			cout << *it << " ";
		cout << endl;
	}
}


int Boundary::getCell(int idx)
{
	return Cell[idx];
}



int Boundary::getOppCell(int idx)
{
	return OppCell[idx];
}


vector<int>::const_iterator Boundary::bgCell()
{
	return Cell.cbegin();
}


vector<int>::const_iterator Boundary::edCell()
{
	return Cell.cend();
}



vector<int>::const_iterator Boundary::bgOppCell()
{
	return OppCell.cbegin();
}


vector<int>::const_iterator Boundary::edOppCell()
{
	return OppCell.cend();
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

