
class Data
{
	public:

		Data() = default;
		Data(string &);

		int InitData(Mesh & mesh);

		// output
		int output(string);

		// data included in the class
		double time;
		int Nx, Ny, Nz;
		Array3D<double> u, v, w, p;
};
