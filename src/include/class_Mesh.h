/*
 * The class to store the variables of computation 
 * domain, and mesh
 */
class Mesh
{
	public:

		Mesh() = default;
		Mesh(int N[3], double L[3]){ InitFlowField(N, L); };

		int InitMesh(int [3], double [3]);

		int addBC(string, int, double, int, double);

	private:

		int Nx, Ny, Nz;
		int NCells;

		int Nxu, Nyu, Nzu;
		int Nxv, Nyv, Nzv;
		int Nxw, Nyw, Nzw;

		double Lx, Ly, Lz;
		double dx, dy, dz;

		map<string, Boundary> BCs;	
		
		Array1D<double> xp, yp, zp;
		Array1D<double> xu, yu, zu;
		Array1D<double> xv, yv, zv;
		Array1D<double> xw, yw, zw;
};
