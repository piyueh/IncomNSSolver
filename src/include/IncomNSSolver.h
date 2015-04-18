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


typedef const int CI;
typedef const unsigned int CUI;
typedef const double CD;
typedef Array3D<double> A3Dd;
typedef const pair<int, double> CPairID;
typedef const array<int, 3> CaryI3;


# include "class_Fluid.h"
# include "class_Boundary.h"
# include "class_Mesh.h"
# include "class_PoissonSolver.h"
# include "class_NSSolverEuler.h"




template<typename T> ostream & operator<<(ostream &os, vector<T> x);
template<typename T> ostream & operator<<(ostream &os, Array3D<T> &A);
ostream &operator<<(ostream &os, Boundary &BC);
ostream &operator<<(ostream &os, Mesh &mesh);
ostream &operator<<(ostream &os, NSSolverEuler &solver);


int tripleLoop(CI &, CI &, CI &, CI &, CI &, CI &, 
		function<void(int &, int &, int &)>);

int dualLoop(CI &, CI &, CI &, CI &, function<void(CI &, CI &)>);
