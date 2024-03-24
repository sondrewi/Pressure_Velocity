// Headr file for SparseMat class

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "SparseAddress.H"

class SparseMat {
 private:
  std::vector<double> entries;

 public:
  SparseAddress& spa;
  Mesh& mesh;
  SparseMat(SparseAddress& sa);
  void add_to_entry(int i, int j, double entry);
  bool symmetric() const;
  bool diagonal_dominance() const;
  bool spd() const;
  double compute_det(int n, int start_row, int start_col,
                     std::vector<int> skip) const;
  double get_entry2D(int i, int j) const;
  double get_entry1D(int i) const;
  Eigen::VectorXd operator*(const Eigen::ArrayXd& v) const;
  void reset();
  const Eigen::ArrayXd mat_underrelax(double alpha);
};
