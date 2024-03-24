#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXi;

// Define BoundaryPatch class
class BoundaryPatch {
public:
    std::string name;
    ArrayXi faces;
    int cell;
    double bc;
};
