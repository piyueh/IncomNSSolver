/*
 * All header files.
 */


# pragma once


# include <iostream>
# include <iomanip>
# include <fstream>
# include <stdexcept>
# include <string>
# include <vector>
# include <map>
# include <utility>
# include <functional>

using namespace std;

# include <eigen3/Eigen/Dense>
# include <eigen3/Eigen/Sparse>

using namespace Eigen;

# include "class_Array3D.h"
# include "class_Fluid.h"
# include "class_Boundary.h"
# include "class_Mesh.h"
# include "class_PoissonSolver.h"
# include "class_NSSolverEuler.h"


ostream &operator<<(ostream &os, Boundary &BC);


int tripleLoop(int &, int &, int &, int &, int &, int &,  
		function<void(int &, int &, int &)>);

int dualLoop(int &, int &, int &, int &,  
		function<void(int &, int &)>);
