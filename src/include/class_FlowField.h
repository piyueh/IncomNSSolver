/*
 * The class to store the variables of computation 
 * domain, mesh, and primative variables.
 */
class Discretization
{
	public:

		

	private:

		int Nx, Ny, Nz;

		int NCells;
		int Nyz;

		double Lx, Ly, Lz;
		double dx, dy, dz;

		vector<Boundary> BCs;	
		
		ArrayXd x, y, z;

		VectorXd u, v, w;
		VectorXd p;
}
