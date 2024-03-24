/*Solver file for the revised SIMPLE algorithm. The main
function below follows the revised algorithm as outlined in the associated
report.*/

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "GSSolver.H"

std::string points_file = "points.dat";
std::string faces_file = "faces.dat";
std::string cells_file = "cells.dat";
std::string boundary_file = "boundary.dat";

std::vector<double> calc_step_norm(Eigen::MatrixXd& old_v,
                                   Eigen::MatrixXd& new_v,
                                   Eigen::ArrayXd& old_p,
                                   Eigen::ArrayXd& new_p);

int main(void) {
  // Initiate mesh, sparse addressing object, sparse matrix objects and sparse
  // linear systems for velocity and pressure.
  Mesh mesh(points_file, faces_file, cells_file, boundary_file);
  SparseAddress SA(mesh);
  SparseMat A(SA);
  SparseMat B(SA);
  SLS SLS_v(A);
  SLS SLS_p(B);

  // Prompt user to enter Reynolds number at which to run simulation:
  std::cout << "Enter Reynolds Number: ";
  double Re;
  std::cin >> Re;

  std::cout << "Enter tolrance for absolute change between solutions: ";
  double tol;
  std::cin >> tol;

  // Convert Reynolds number to a n int
  int int_Re = static_cast<int>(Re);
  std::string Re_string = std::to_string(int_Re);

  // Calculate kinematic viscosity for this case of lid-driven cavity:
  double nu_d = 0.1 / Re;

  // Initiate file to contain infinity norm of residuals in inner loops
  // and the absolute change, between interations, in inner and outer solutions
  // File will thus have 9 columns: residuals, norm of change in outer
  // solutions, norm of change in inner solutions, in each variable
  std::ofstream norms("Norms_Re" + Re_string + ".dat");

  // Initiate vectors for the max absolute change, between iterations,
  // in each variable for outer and inner solutions
  std::vector<double> step_max_delta_o = {1, 1, 1};
  std::vector<double> step_max_delta_i;

  // Define boundary conditions. Note that a boundary condition type 0
  // is defined to correspond to Dirichlet type while 1 will correspond to
  // Neumann (see SLS.cpp set_diffusion_mat function for use of this convention)
  Eigen::MatrixXd boundary_values_v = Eigen::MatrixXd::Zero(3, 4);
  boundary_values_v(0, 2) = 1.0;

  // Boundary values for p and velocity. Boundary patches were indexed as
  //  left, right, top, bottom for the lid-driven cavity
  std::vector<double> boundary_values_p = {0, 0, 0, 0};
  std::vector<double> boundary_values_vx = {0, 0, 1, 0};
  std::vector<double> boundary_values_vy = {0, 0, 0, 0};
  std::vector<double> boundary_values_vz = {0, 0, 0, 0};

  // Boundary types for pressure and velocity
  std::vector<int> btype_v = {0, 0, 0, 0};
  std::vector<int> btype_p = {1, 1, 1, 1};

  // Set under-relaxation factor for velocity
  const double alpha_u = 0.9;

  // Initiate the smoother and solver objects
  // Here we carry out only 200 iterations for each linear system
  GSSmooth smoother_v(SLS_v);
  GSSolver solver_v(smoother_v, 200, 1);

  GSSmooth smoother_p(SLS_p);
  GSSolver solver_p(smoother_p, 200, 1);

  // Set up large matrix of solutions.
  Eigen::MatrixXd v_sols = Eigen::MatrixXd::Zero(mesh.getNumCells(), 3);
  Eigen::ArrayXd p_sols = Eigen::ArrayXd::Zero(mesh.getNumCells());

  // Initiate Arrays that will hold diagonal coefficients in velocity system
  // and their inverse
  Eigen::ArrayXd diag, diag_inv;

  // Set initial zero solutions in mesh
  mesh.set_velocity(v_sols);
  mesh.set_pressure(p_sols);
  mesh.calc_face_flux(btype_v, boundary_values_v);

  // Define vector of kinematic viscosities. This vector could contain ´varying
  // viscosities throughout the domain for different applications
  Eigen::VectorXd nu = nu_d * Eigen::VectorXd::Ones(mesh.getNumCells());

  // Define cell volume:
  double vol = mesh.get_cell_vol(1004);

  // Define counter to count number of iterations required for convergence
  int counter = 0;

  // Define old outer solutions that will be used to gauge change in solutions
  // from one iteration to the next
  Eigen::MatrixXd old_v_o = Eigen::MatrixXd::Ones(mesh.getNumCells(), 3);
  Eigen::ArrayXd old_p_o = Eigen::ArrayXd::Ones(mesh.getNumCells());

  // Define old inner solutions that will be used to gauge change in solutions
  // from one iteration to the next
  Eigen::MatrixXd old_v_i = Eigen::MatrixXd::Ones(mesh.getNumCells(), 3);
  Eigen::ArrayXd old_p_i = Eigen::ArrayXd::Ones(mesh.getNumCells());

  // Solve with tolerance of 10^(-5) on the max absolute change in outer
  // solutions between iterations in all variables
  while (*max_element(step_max_delta_o.begin(), step_max_delta_o.end()) > tol) {
    counter += 1;

    // Get the old outer solutions
    old_v_o.col(0) = mesh.get_field(Field::vx).matrix();
    old_v_o.col(1) = mesh.get_field(Field::vy).matrix();
    old_p_o = mesh.get_field(Field::p);

    // We will not solve for the vx predictor
    //  Reset A, x and b vector in sparse system to zero
    SLS_v.reset_x();
    SLS_v.reset_b();
    A.reset();

    // Discretise convection and diffusion operators as described in
    // SLS.cpp and add discretisations to A and b as appropriate.
    SLS_v.set_Diffusion_mat(boundary_values_vx, btype_v, nu, Field::vx, false);
    SLS_v.set_convection_mat(boundary_values_vx, btype_v, true);

    // Implicitly under-relax system
    diag = SLS_v.underrelax_LHS(alpha_u);
    SLS_v.underrelax_RHS(alpha_u, diag, Field::vx);
    diag_inv = 1.0 / diag;

    // Set initial Ax vector and residual
    SLS_v.set_Ax();
    SLS_v.set_res();

    // Solve using 200 sweeps of Gauss-Seidel
    double res_norm_vx = solver_v.solve(1e-15);

    // Assign to first column of v_sols the vx predictor
    for (int j = 0; j < mesh.getNumCells(); j++) {
      double sol = SLS_v.get_x_entry(j);
      v_sols(j, 0) = sol;
    }

    // Reset x and b vector in sparse system to zero.
    // A-matrix will be same for y_velocity. RHS b will
    // be zero for y-velocity due to zero boundary conditions
    SLS_v.reset_x();
    SLS_v.reset_b();

    // Under-relax RHS of vy equation
    SLS_v.underrelax_RHS(alpha_u, diag, Field::vy);

    SLS_v.set_Ax();
    SLS_v.set_res();

    // Solve for y-velocity predictor
    double res_norm_vy = solver_v.solve(1e-15);

    // Set second column in v_sols
    for (int j = 0; j < mesh.getNumCells(); j++) {
      double sol = SLS_v.get_x_entry(j);
      v_sols(j, 1) = sol;
    }

    // Set mesh velocity using predictors (z-velocity always zero)
    mesh.set_velocity(v_sols);

    // Calculate face fluxes using above solution to velocities
    mesh.calc_face_flux(btype_v, boundary_values_v);

    /*Now calculate pressure predictor from pressure equation*/
    B.reset();
    SLS_p.reset_x();
    SLS_p.reset_b();

    // Set the Laplacian discretisation matrix B. The last argument indicates
    // Laplacian
    SLS_p.set_Diffusion_mat(boundary_values_p, btype_p, diag_inv * vol,
                            Field::p, true);

    // Explicitly set RHS of pressure equation using cell fluxes
    for (int j = 0; j < mesh.getNumCells(); j++) {
      double cell_flux = mesh.get_cell_netflux(j);
      SLS_p.add_to_b_entry(j, cell_flux);
    }

    SLS_p.set_Ax();
    SLS_p.set_res();

    // Solve for pressure predictor
    double res_norm_p = solver_p.solve(1e-15);

    // Get p solutions from sparse linear system
    for (int j = 0; j < mesh.getNumCells(); j++) {
      double sol = SLS_p.get_x_entry(j);
      p_sols(j) = sol;
    }

    // Calaculate convergence parameter for inner solutions
    step_max_delta_i = calc_step_norm(old_v_i, v_sols, old_p_i, p_sols);
    old_v_i.col(0) = v_sols.col(0);
    old_v_i.col(1) = v_sols.col(1);
    old_p_i = p_sols;

    // Set pressure in mesh and correct face fluxes
    mesh.set_pressure(p_sols);
    mesh.correct_face_flux(diag_inv * vol, btype_p, boundary_values_p);

    // Explicitly under-relax pressure
    p_sols = old_p_o + (1.1 - alpha_u) * (p_sols - old_p_o);

    // Set mesh pressure again
    mesh.set_pressure(p_sols);

    // Calculate volume integrals of pressure gradients
    mesh.calc_gauss_gradients(btype_p, boundary_values_p, Field::p);

    // Correct velocity
    mesh.correct_velocity(diag_inv);

    // Gather new solutions for velocity
    v_sols.col(0) = mesh.get_field(Field::vx).matrix();
    v_sols.col(1) = mesh.get_field(Field::vy).matrix();

    // In each variable, calculate the maximal absolute change
    // compared to last outer solution
    step_max_delta_o = calc_step_norm(old_v_o, v_sols, old_p_o, p_sols);

    // Write to norm file:
    norms << res_norm_vx << "  " << res_norm_vy << "  " << res_norm_p << "  "
          << step_max_delta_o[0] << "  " << step_max_delta_o[1] << "  "
          << step_max_delta_o[2] << "  " << step_max_delta_i[0] << "  "
          << step_max_delta_i[1] << "  " << step_max_delta_i[2] << std::endl;

    // Output max absolute changes in each variable every 10 iterations
    if (counter % 10 == 0) {
      std::cout << "Iteration " << counter << " finished:" << std::endl;
      std::cout << "Max delta vx:" << step_max_delta_o[0] << std::endl;
      std::cout << "Max delta vy:" << step_max_delta_o[1] << std::endl;
      std::cout << "Max delta p:" << step_max_delta_o[2] << std::endl;
      std::cout << std::endl;
    }
  }

  // Initiate file for final solutions
  std::ofstream Sols("Sols_Re" + Re_string + ".dat");

  for (int i = 0; i < mesh.getNumCells(); i++) {
    Sols << mesh.get_cell_centroid(i)(0) << "  " << mesh.get_cell_centroid(i)(1)
         << "  " << v_sols(i, 0) << "  " << v_sols(i, 1) << "  " << p_sols(i)
         << std::endl;
  }

  std::cout << "Simulation Finished in " << counter << " iterations"
            << std::endl;
}

// Calculate infinity norm in change between new and old solutions
std::vector<double> calc_step_norm(Eigen::MatrixXd& old_v,
                                   Eigen::MatrixXd& new_v,
                                   Eigen::ArrayXd& old_p,
                                   Eigen::ArrayXd& new_p) {
  std::vector<double> step_norms(3);

  Eigen::ArrayXd vx_old = old_v.col(0).array();
  Eigen::ArrayXd vx_new = new_v.col(0).array();

  Eigen::ArrayXd vy_old = old_v.col(1).array();
  Eigen::ArrayXd vy_new = new_v.col(1).array();

  step_norms[0] = (vx_new - vx_old).abs().maxCoeff();
  step_norms[1] = (vy_new - vy_old).abs().maxCoeff();
  step_norms[2] = (new_p - old_p).abs().maxCoeff();

  return step_norms;
}
