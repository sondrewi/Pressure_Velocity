// Header file for Sparse Linear System Class

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SparseMat.H"

using Field = Mesh::Field;

class SLS {
 private:
  Eigen::ArrayXd x, b, res, Ax;

 public:
  SparseMat& A;
  Mesh& mesh;
  SparseAddress& spa;
  SLS(SparseMat& SM);
  void set_convection_mat(std::vector<double>& boundary_values,
                          std::vector<int>& bc_type, bool up_diff);
  void set_temporal_mat(double dt, Field field);
  void set_Diffusion_mat(const std::vector<double>& boundary_values,
                         const std::vector<int>& bc_type,
                         const Eigen::ArrayXd& scale_vec, Field field,
                         bool lap);
  void set_res();
  double calc_res_norm(int norm);
  double norm_factor();
  double Scaled_norm(int norm);
  double get_x_entry(int i) const;
  void set_x_entry(int i, double value);
  double get_b_entry(int i) const;
  void set_b_entry(int i, double value);
  void add_to_b_entry(int i, double value);
  void set_Ax();
  Eigen::ArrayXd get_Ax();
  Eigen::ArrayXd get_res();
  double calc_b_norm(int norm);
  void reset_x();
  void reset_b();
  Eigen::ArrayXd underrelax_LHS(double alpha);
  void underrelax_RHS(double alpha, Eigen::ArrayXd& diag, Field field);
};
