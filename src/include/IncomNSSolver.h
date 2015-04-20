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
# include <cmath>

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
typedef array<int, 3> aryI3;
typedef array<double, 3> aryD3;
typedef const array<int, 3> CaryI3;
typedef const array<double, 3> CaryD3;
typedef vector<double> VD;


# include "class_Fluid.h"
# include "class_Boundary.h"
# include "class_Mesh.h"
# include "class_Data.h"
# include "class_PoissonSolver.h"
# include "class_NSSolverEuler.h"




template<typename T> ostream & operator<<(ostream &os, vector<T> x);
template<typename T> ostream & operator<<(ostream &os, Array3D<T> &A);
ostream &operator<<(ostream &os, Boundary &BC);
ostream &operator<<(ostream &os, Mesh &mesh);
ostream &operator<<(ostream &os, NSSolver &solver);


int tripleLoop(CI &, CI &, CI &, CI &, CI &, CI &, 
		function<void(CI &, CI &, CI &)>);

int tripleLoop(CI &, CI &, CI &, CI &, CI &, CI &, 
		CD &, function<void(CI &, CI &, CI &, CD &)>);

int dualLoop(CI &, CI &, CI &, CI &, function<void(CI &, CI &)>);

int TGVortex(CI & Nxu, CI & Nyu, CI & Nzu, A3Dd & u, VD & xu, VD & yu,
		CI & Nxv, CI & Nyv, CI & Nzv, A3Dd & v, VD & xv, VD & yv, CD & t);

double evalRelErr(VectorXd & x, VectorXd & xe);
