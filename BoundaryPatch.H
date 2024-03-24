#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Define BoundaryPatch class
class BoundaryPatch {
 public:
  std::string name;      // name of boundary patch
  Eigen::ArrayXi faces;  // Faces adjoining boundary patch
  double bc;             // Optional boundary condition
  int cell;              // owner cell of boundary patch
};
