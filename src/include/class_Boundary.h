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
		 * {Nx, Ny, Nz}, Dir, sign, 
		 * {p BC, Value}, {u BC, Value}, {v BC, Value}, {w BC, Value}
		 */
		Boundary(CaryI3, CUI, CI, CPairID, CPairID, CPairID, CPairID);

		Boundary & operator=(const Boundary &);

		CI & get_Dir() const { return dir; }
		CI & get_Sign() const { return sign; }

		CI & get_pType() const { return pType; }
		CI & get_uType() const { return uType; }
		CI & get_vType() const { return vType; }
		CI & get_wType() const { return wType; }

		CD & get_pBCvalue() const { return pBCvalue; }
		CD & get_uBCvalue() const { return uBCvalue; }
		CD & get_vBCvalue() const { return vBCvalue; }
		CD & get_wBCvalue() const { return wBCvalue; }

		CI & get_pBCIdx() const { return pBCIdx; }
		CI & get_uBCIdx() const { return uBCIdx; }
		CI & get_vBCIdx() const { return vBCIdx; }
		CI & get_wBCIdx() const { return wBCIdx; }

		CI & get_pBCcorIdx() const { return pBCcorIdx; }
		CI & get_uBCcorIdx() const { return uBCcorIdx; }
		CI & get_vBCcorIdx() const { return vBCcorIdx; }
		CI & get_wBCcorIdx() const { return wBCcorIdx; }


		int updGhost(CI &, CI &, CI &, A3Dd &, CD &);

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

		array<function<void(CI &, CI &, A3Dd &, CD &)>, 4> updOneGh;

		int InitUpdGh();
};


