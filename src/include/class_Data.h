
class Data
{
	public:

		Data() = default;
		Data(const string &, Mesh & mesh);

		int InitData(Mesh & mesh);

		// output
		int output(string);

		// data included in the class
		double time;

		array<array<int, 3>, 4> N;

		Array3D<double> u, v, w, p;

		array<Array3D<double> *, 4> PVar = {&p, &u, &v, &w};
	
	private:

};
