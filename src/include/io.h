
template<typename T>
ostream & operator<<(ostream &os, vector<T> x)
{
	for(auto &i: x) cout << i << " ";
	cout << endl;
	return os;
}




