
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
		int Nx, Ny, Nz;
		Array3D<double> u, v, w, p;
	
	private:

		int SetBCvalues(Mesh & mesh);
};
