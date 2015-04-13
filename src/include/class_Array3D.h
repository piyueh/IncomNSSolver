/*
 * This defines a class that use 1D vector to represent 3D array.
 */

template<typename T>
class Array3D: public vector<T>
{
	public:

		Array3D()=default;


		T & operator()(int i, int j, int k)
		{
			i -= xbgIdx; j -= ybgIdx; k -= zbgIdx;	
			return *(this->_M_impl._M_start + i*Nyz + j*Nz + k); 
		}


		void initShape(int N1, int N2, int N3) 
		{
			Nx = N1; Ny = N2; Nz = N3; Nyz = N2 * N3;
			xedIdx = Nx - 1; yedIdx = Ny - 1; zedIdx = Nz - 1;

			this->resize(N1 * N2 * N3); 
		}


		void initShape(int Idx1a, int Idx1b,
				int Idx2a, int Idx2b,int Idx3a, int Idx3b) 
		{
			xbgIdx = Idx1a; xedIdx = Idx1b;
			ybgIdx = Idx2a; yedIdx = Idx2b;
			zbgIdx = Idx3a; zedIdx = Idx3b;

			Nx = xedIdx - xbgIdx + 1; 
			Ny = yedIdx - ybgIdx + 1; 
			Nz = zedIdx - zbgIdx + 1; 
			Nyz = Ny * Nz;

			this->resize(Nx * Ny * Nz); 
		}


		void setConstant(T value){
			for(auto &i: *this) i = value;}


		void setZeros() {
			for(auto &i: *this) i = 0;}


		vector<int> shape(){
			return vector<int> {Nx, Ny, Nz};}
		

	private:

		int Nx, Ny, Nz;
		int Nyz;

		int xbgIdx = 0, xedIdx;
		int ybgIdx = 0, yedIdx;
		int zbgIdx = 0, zedIdx;
};


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
