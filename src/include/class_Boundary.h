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
		Boundary(int, string, vector<int>, int, double);

		vector<int>::const_iterator bgCell();
		vector<int>::const_iterator edCell();

		int getNcells();
		int getType();
		double getBCvalue();
		string getDirection();

		void print();

	private:

		int Ncells;
		string direction;
		vector<int> Cell;
		int Type;
		double value;

};
