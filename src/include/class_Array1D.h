/*
 * This defines a class that use 1D vector to represent 1D array.
 */

template<typename T>
class Array1D: public vector<T>
{
	public:

		Array1D()=default;


		T & operator()(int i, int j, int k)
		{
			i -= bgIdx;	
			return *(this->_M_impl._M_start + i); 
		}


		void initShape(int N1) 
		{
			N = N1; edIdx = N - 1;
			this->resize(N); 
		}


		void initShape(int Idx1a, int Idx1b,
				int Idx2a, int Idx2b,int Idx3a, int Idx3b) 
		{
			bgIdx = Idx1a; edIdx = Idx1b;
			N = edIdx - bgIdx + 1; 
			this->resize(N); 
		}


		void setConstant(T value){
			for(auto &i: *this) i = value;}


		void setZeros() {
			for(auto &i: *this) i = 0;}

	private:

		int N;
		int bgIdx = 0, edIdx;
};


