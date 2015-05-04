# include "include/IncomNSSolver.h"


int main(int argc, char *argv[])
{

	map<string, string> Files;

	for(int i=1; i<argc; ++i)
	{
		if ((strcmp(argv[i], "-f") == 0) ||  (strcmp(argv[i] ,"-m") == 0) || 
				(strcmp(argv[i] ,"-d") == 0) || (strcmp(argv[i] ,"-c") == 0)) 
		{
			Files[string(argv[i])] = string(argv[i+1]);
			i += 1;
		} 
		else
		{
			cerr << "Invalid parameter: " << argv[i] << endl;
			return 1;
		}
	}

	Fluid fluid(Files["-f"]);
	Mesh mesh(Files["-m"]);
	Data data(Files["-d"], mesh);

	Solid cylinder({10, 10}, 0.5, mesh);

	cout << mesh << endl;

	NSSolver solver(mesh, fluid, data, cylinder, Files["-c"]);

	cout << solver << endl;

	solver.solve();
	
	data.output("Data.txt");

	return 0;
}


