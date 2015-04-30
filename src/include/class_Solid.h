class Solid
{
	public:
		Solid(const array<double, 2> & Center, CD & Radius, const Mesh & mesh);

	private:

		array<double, 2> center;
		double R;

		vector<array<int, 3>> uIdxBC, vIdxBC, wIdxBC;

		vector<array<int, 3>> uCorIdxBC, vCorIdxBC, wCorIdxBC;

		vector<array<double, 2>> uDist, vDist, wDist;

		Array3D<int> uFlag, vFlag, wFlag;

		int determineFlags(const Mesh &);
		int determineIdx(const Mesh &);

};
