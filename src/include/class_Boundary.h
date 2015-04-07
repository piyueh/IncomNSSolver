# include <string>
# include <vector>

using namespace std;


/*
 * The class to store data of each boundary.
 */
class Boundary
{
	public:

		Boundary() = default;
		Boundary(int, string, int, double, vector<int>);
		Boundary(int, string, int, double, vector<int>, vector<int>);

		vector<int>::const_iterator bgCell();
		vector<int>::const_iterator edCell();
		vector<int>::const_iterator bgOppCell();
		vector<int>::const_iterator edOppCell();

		int getNcells();
		int getType();
		double getBCvalue();
		string getDirection();

		int getCell(int);
		int getOppCell(int);

		void print();

	private:

		int Ncells;
		string direction;
		int Type;
		double value;
		vector<int> Cell;
		vector<int> OppCell;

};
