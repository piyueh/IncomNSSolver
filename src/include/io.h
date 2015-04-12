
ostream & operator<<(ostream &os, vector<int> x)
{
	for(auto &i: x) cout << i << " ";
	cout << endl;
	return os;
}


template<typename T>
ostream & operator<<(ostream &os, Array3D<T> &A)
{
	os << "Nx Ny Nz: "
	   << A.shape() << endl;
	
	os << "Values" << endl;
	for(auto &i: A) cout << i << " ";
	cout << endl;


	return os;
}


ostream &operator<<(ostream &os, Boundary &BC)
{
	os << "# of cells in this boundary: " << BC.get_Ncells() << endl;
	os << "The normal direction: " << BC.get_Direction() << endl;
	os << "The type of pressure BC: " << BC.get_pType() << endl;
	os << "The type of velocity BC: " << BC.get_vType() << endl;
	os << "The value of pressure BC: " << BC.get_pBCvalue() << endl;
	os << "The value of velocity BC: " << BC.get_vBCvalue() << endl;

	os << "The index of cells in this boundary (in 1D index): " << endl;
	for(auto it=BC.bgCell(); it<BC.edCell(); ++it)
		os << *it << " ";
	os << endl;

	os << "The index of cells in opposite boundary (in 1D index): " << endl;
	for(auto it=BC.bgOppCell(); it<BC.edOppCell(); ++it)
		os << *it << " ";
	os << endl;
	return os;
}
