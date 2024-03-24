// Header file for Sparse Address class

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Mesh.H"

class SparseAddress {
 private:
  int nr_equations;
  int matrix_band;

  std::vector<int> col_index;
  std::vector<int> row_index;

 public:
  Mesh& mesh;
  SparseAddress(Mesh& the_mesh);
  int get_nr_eqs() const;
  int nr_off_diag() const;
  int get_col_index(int) const;
  int get_row_index(int) const;
  int mat_band() const;
  std::pair<int, int> oneD_to_twoD(int i) const;
  int twoD_to_oneD(int i, int j, bool trunc = false) const;
};
