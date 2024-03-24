// Header file for Gauss-Seidel Solver Class

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "GSSmooth.H"

class GSSolver {
 private:
  int max_it;
  int sweeps;

 public:
  SLS& sls;
  GSSmooth& smoother;
  GSSolver(GSSmooth& smoother_object, int sweeps_, int maxit);
  double solve(double relTol);
};
