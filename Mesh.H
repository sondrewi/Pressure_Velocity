// Header file for Mesh.cpp file

#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "BoundaryPatch.H"
#include "Cell.H"
#include "Face.H"
#include "Point.H"

class Mesh {
 private:
  std::vector<Point> points;  // vector of points making up mesh
  std::vector<Face> faces;    // vector of faces in the mesh
  std::vector<Cell> cells;    // vector of cells in the mesh
  std::vector<BoundaryPatch>
      boundaryPatches;  // vector of boundary patches enclosing mesh

  void readPoints(const std::string& fileName);
  void readFaces(const std::string& fileName);
  void readCells(const std::string& fileName);
  void readBoundaryPatches(const std::string& fileName);
  void assign_face_owners();
  void assign_cell_neighbours();
  void assign_boundary_cells();
  Eigen::VectorXd getCellVolumes() const;
  Eigen::Vector3d face_raw_avg(int face_nr) const;
  void Set_Sf(int face_nr);
  std::vector<std::array<int, 2> > triangle_list(int face_nr) const;
  std::vector<Eigen::Vector3d> triangle_midpoints(int face_nr) const;
  Eigen::VectorXd triangle_areas(int face_nr) const;
  void set_face_centroid(int face_nr);
  Eigen::Vector3d cell_raw_avg(int cell_nr) const;
  Eigen::VectorXd pyramid_volumes(int cell_nr) const;
  void set_cell_vol(int cell_nr);
  void set_cell_centroid(int cell_nr);
  void set_face_interpolation_factor(int face_nr);
  void set_delta_coef(int face_nr);
  void set_cell_velocity(int cell_nr, double vx, double vy, double vz);

 public:
  // Mesh will have a given number of fields. They are set to vx, vy, vz and p
  // (pressure) in this case but this list can be expanded and altered depending
  // on application
  enum Field { vx, vy, vz, p };
  // Constructor to read mesh data from files
  Mesh(const std::string& pointsFile, const std::string& facesFile,
       const std::string& cellsFile, const std::string& boundaryFile);
  int getNumPoints() const;
  int getNumFaces() const;
  int getNumCells() const;
  int getNumBoundaryPatches() const;
  double get_cell_vol(int cell_nr) const;
  Eigen::Vector3d get_cell_velocity(int cell_nr) const;
  bool get_face_empty(int face_nr) const;
  int get_face_owner(int face_nr) const;
  int get_face_nbr(int face_nr) const;
  Eigen::Vector3d get_face_normal(int face_nr) const;
  int face_boundary(int face_nr) const;
  Eigen::Vector3d get_cell_centroid(int cell_nr) const;
  Eigen::Vector3d get_face_centroid(int face_nr) const;
  double get_face_interpolation_factor(int face_nr) const;
  double get_delta_coef(int face_nr) const;
  const Eigen::Vector3d get_cell_grad(int cell_nr, Field field) const;
  std::vector<int> get_cell_nbrs(int cell_nr) const;
  double get_cell_netflux(int cell_nr) const;
  double get_face_flux(int face_nr) const;
  double get_cell_val(int cell_nr, Field field) const;
  void calc_gauss_gradients(std::vector<int>& bc_type,
                            std::vector<double>& boundary_values, Field field);
  void calc_face_flux();
  void set_velocity(Eigen::MatrixXd& v);
  void set_pressure(Eigen::ArrayXd& p);
  void set_cell_pressure(int cell_nr, double p);
  void correct_velocity(Eigen::VectorXd scaler);
  const Eigen::ArrayXd get_field(Field field) const;
  void correct_face_flux(Eigen::VectorXd scaler, std::vector<int> b_type,
                         std::vector<double> boundary_values_p);
  void calc_face_flux(std::vector<int> b_type,
                      Eigen::MatrixXd boundary_values_v);
  void set_cell_grad(int cell_nr, Field field, Eigen::Vector3d grad);
};
