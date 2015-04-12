/*
 * The class to store data of each boundary.
 */
class Boundary
{
	public:

		Boundary() = default;

		/*
		 * Input:
		 * Nx, Ny, Nz, Direction, P BC, P BC Value, V BC, V BC Value
		 */
		Boundary(int, int, int, string, int, double, int, double);


		int get_Ncells() { return Ncells; }

		string get_Direction() { return direction; }

		int get_pType() { return pType; }
		int get_vType() { return vType; }

		double get_pBCvalue() { return pBCvalue; }
		double get_vBCvalue() { return vBCvalue; }

		int getCell(int idx) { return Cell[idx]; }
		int getOppCell(int idx) { return OppCell[idx]; }

		vector<int>::const_iterator bgCell() { return Cell.cbegin(); }
		vector<int>::const_iterator edCell() { return Cell.cend(); }
		vector<int>::const_iterator bgOppCell() { return OppCell.cbegin(); }
		vector<int>::const_iterator edOppCell() { return OppCell.cend(); }

		void print();

	private:

		int Ncells;

		string direction;

		int pType, vType;

		double pBCvalue, vBCvalue;

		vector<int> Cell, OppCell;

};

