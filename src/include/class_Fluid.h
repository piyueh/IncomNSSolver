/*
 * Fluid properties
 */


# pragma once


class Fluid
{
	public:

		Fluid(double a, double b): nu(a), rho(b) {};

		Fluid(const string & fName)
		{
			ifstream file(fName);
			string Parm;

			while(getline(file, Parm))
			{
				istringstream line(Parm);
				string var;

				line >> var;

				if (var == "nu") line >> nu;
				if (var == "rho") line >> rho;
			}

			file.close();
		}

		// viscosity
		double nu;

		// density
		double rho;
};
