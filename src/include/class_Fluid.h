/*
 * Fluid properties
 */
class Fluid
{
	public:
		Fluid(double a, double b): nu(a), rho(b) {};
		// viscosity
		double nu;

		// density
		double rho;
};
