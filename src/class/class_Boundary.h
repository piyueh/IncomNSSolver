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
		void print();

	private:
		int Ncells;
		string direction;
		vector<int> Cell;
		int Type;
		double value;

};
