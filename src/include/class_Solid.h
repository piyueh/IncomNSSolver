class Solid
{
	public:

		Solid(const array<double, 2> & Center, CD & Radius, const Mesh & mesh);

		int updVelocity(A3Dd &u, A3Dd &v, A3Dd &w);

		int output(string fileName);

	private:

		array<double, 2> center;
		double R;

		int NPu, NPv, NPw;

		vector<array<int, 3>> uIdxBC, vIdxBC, wIdxBC;

		vector<array<int, 3>> uCorIdxBC, vCorIdxBC, wCorIdxBC;

		vector<array<double, 2>> uDist, vDist, wDist;

		Array3D<int> uFlag, vFlag, wFlag;

		int determineFlags(const Mesh &);
		int determineIdx(const Mesh &);

};
