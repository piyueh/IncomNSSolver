/*
 * The class to store data of each boundary.
 */


# pragma once


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
		Boundary(const array<int ,3>, const unsigned int, const int, 
				const pair<int, double>, const pair<int, double>, 
				const pair<int, double>, const pair<int, double>);

		Boundary & operator=(const Boundary &);

		const int & get_Dir() const { return dir; }
		const int & get_Sign() const { return sign; }

		const int & get_pType() const { return pType; }
		const int & get_uType() const { return uType; }
		const int & get_vType() const { return vType; }
		const int & get_wType() const { return wType; }

		const double & get_pBCvalue() const { return pBCvalue; }
		const double & get_uBCvalue() const { return uBCvalue; }
		const double & get_vBCvalue() const { return vBCvalue; }
		const double & get_wBCvalue() const { return wBCvalue; }

		const int & get_pBCIdx() const { return pBCIdx; }
		const int & get_uBCIdx() const { return uBCIdx; }
		const int & get_vBCIdx() const { return vBCIdx; }
		const int & get_wBCIdx() const { return wBCIdx; }

		const int & get_pBCcorIdx() const { return pBCcorIdx; }
		const int & get_uBCcorIdx() const { return uBCcorIdx; }
		const int & get_vBCcorIdx() const { return vBCcorIdx; }
		const int & get_wBCcorIdx() const { return wBCcorIdx; }

	private:

		int dir, sign;

		array<int, 4> BCtype, BCIdx, BCcorIdx;
		array<double, 4> BCvalues;

		int & pType = BCtype[0], & uType = BCtype[1], 
			& vType = BCtype[2], & wType = BCtype[3];
 
		int & pBCIdx = BCIdx[0], & pBCcorIdx =BCcorIdx[0];
		int & uBCIdx = BCIdx[1], & uBCcorIdx =BCcorIdx[1];
		int & vBCIdx = BCIdx[2], & vBCcorIdx =BCcorIdx[2];
		int & wBCIdx = BCIdx[3], & wBCcorIdx =BCcorIdx[3];

		double & pBCvalue = BCvalues[0], & uBCvalue = BCvalues[1], 
			   & vBCvalue = BCvalues[2], & wBCvalue = BCvalues[3];
};


