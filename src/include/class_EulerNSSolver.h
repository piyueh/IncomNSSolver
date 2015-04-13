/*
 * The Navier-Stoke solver using the basic 1st-order Euler time marching.
 */
class NSSolverEuler
{
	public:

	private:

		double time;

		Array3D<double> u, v, w, p;
	p.initShape(-1, Nx, -1, Nx, -1, Nz);
	u.initShape(-1, Nxu, -1, Nyu, -1, Nzu);
	v.initShape(-1, Nxv, -1, Nyv, -1, Nzv);
	w.initShape(-1, Nxw, -1, Nyw, -1, Nzw);
	p.setZeros();
	u.setZeros(); v.setZeros(); w.setZeros();
};
