#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Define Cell class header
class Cell {
 public:
  Eigen::ArrayXi faces;      // Indices of faces enclosing cell
  std::vector<int> nbrs;     // Indices of neighbours enclosing cell
  Eigen::Vector3d centroid;  // volumetric centroid of cell
  double vol;                // volume of cell
  Eigen::Vector3d vel;       // velocity at cell centroid
  double p;                  // pressure at cell centroid
  Eigen::Vector3d grad_p;    // gradient of pressure at cell centroid
  Eigen::Vector3d grad_vx;   // Velocity gradients at cell centroid
  Eigen::Vector3d grad_vy;
  Eigen::Vector3d grad_vz;
};
