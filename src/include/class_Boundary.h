# include <iostream>
# include <string>
# include <stdexcept>
# include <utility>

using namespace std;

/*
 * The class to store data of each boundary.
 */
class Boundary
{
	friend class Mesh;
	friend ostream &operator<<(ostream &, Boundary &);

	public:

		Boundary() = default;

		/*
		 * Input:
		 * Nx, Ny, Nz, Direction, P BC, P BC Value, V BC, V BC Value
		 */
		Boundary(int Nx, int Ny, int Nz, unsigned int Dir, int Sign, 
				pair<int, double> p, pair<int, double> u, 
				pair<int, double> v, pair<int, double> w);


		int get_Direction() { return sign*dir; }

		int get_pType() { return pType; }
		int get_uType() { return uType; }
		int get_vType() { return vType; }
		int get_wType() { return wType; }

		double get_pBCvalue() { return pBCvalue; }
		double get_uBCvalue() { return uBCvalue; }
		double get_vBCvalue() { return vBCvalue; }
		double get_wBCvalue() { return wBCvalue; }

		int get_pBCIdx() { return pBCIdx; }
		int get_uBCIdx() { return uBCIdx; }
		int get_vBCIdx() { return vBCIdx; }
		int get_wBCIdx() { return wBCIdx; }

		int get_pBCcorIdx() { return pBCcorIdx; }
		int get_uBCcorIdx() { return uBCcorIdx; }
		int get_vBCcorIdx() { return vBCcorIdx; }
		int get_wBCcorIdx() { return wBCcorIdx; }

	private:

		int dir, sign;

		int pType, uType, vType, wType;

		double pBCvalue, uBCvalue, vBCvalue, wBCvalue;

		int pBCIdx, pBCcorIdx;
		int uBCIdx, uBCcorIdx;
		int vBCIdx, vBCcorIdx;
		int wBCIdx, wBCcorIdx;

};


