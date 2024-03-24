#include "SLS.H"

// Constructor for sparse linear system
SLS::SLS(SparseMat& SM) : A(SM), mesh(SM.mesh), spa(SM.spa) {
  x = Eigen::ArrayXd::Zero(mesh.getNumCells());
  b = Eigen::ArrayXd::Zero(mesh.getNumCells());
  Ax = Eigen::ArrayXd::Zero(mesh.getNumCells());
  res = Eigen::ArrayXd::Zero(mesh.getNumCells());
}

// Function to set the residual
void SLS::set_res() {
  for (int i = 0; i < x.size(); i++) {
    res(i) = b(i) - Ax(i);
  }
}

// Function to calculate the residual norm (1-norm, 2-norm, or inf norm)
double SLS::calc_res_norm(int norm) {
  double s = 0;

  if (norm == 1) {
    for (int i = 0; i < res.size(); i++) {
      s += abs(res(i));
    }
  }

  else if (norm == 2) {
    for (int i = 0; i < res.size(); i++) {
      s += pow(res(i), 2);
    }
    s = sqrt(s);
  }

  else {
    for (int i = 0; i < res.size(); i++) {
      if (abs(res(i)) > s) {
        s = abs(res(i));
      }
    }
  }

  return s;
}

Eigen::ArrayXd SLS::get_Ax() { return Ax; }
Eigen::ArrayXd SLS::get_res() { return res; }

// Function to calculate the RHS norm (1-norm, 2-norm, or inf norm)
double SLS::calc_b_norm(int norm) {
  double s = 0;

  if (norm == 1) {
    for (int i = 0; i < b.size(); i++) {
      s += abs(b(i));
    }
  }

  else if (norm == 2) {
    for (int i = 0; i < b.size(); i++) {
      s += pow(b(i), 2);
    }
    s = sqrt(s);
  }

  else {
    for (int i = 0; i < b.size(); i++) {
      if (abs(b(i)) > s) {
        s = abs(b(i));
      }
    }
  }

  return s;
}

// Function to under-relax LHS of sparse linear system
// alpha is under-relaxation factor
Eigen::ArrayXd SLS::underrelax_LHS(double alpha) {
  Eigen::ArrayXd diag = A.mat_underrelax(alpha);

  return diag;
}

// Under-relax RHS of sparse linear system
// Add weighted values of the old field for which we are solving
void SLS::underrelax_RHS(double alpha, Eigen::ArrayXd& diag, Field field) {
  b += ((1 - alpha) / alpha) * diag * mesh.get_field(field);
}

// Find factor n that quantifies deviation of entries in residual
// relative to a reference state of the residual (an average)
double SLS::norm_factor() {
  double n = 0;
  double avg = 0;

  for (int i = 0; i < res.size(); i++) {
    avg += x[i];
  }

  avg /= x.size();

  Eigen::VectorXd x_ref = avg * Eigen::VectorXd::Ones(x.size());
  Eigen::VectorXd Ax_ref = A * x_ref;

  for (int i = 0; i < res.size(); i++) {
    n += abs(Ax(i) - Ax_ref(i)) - abs(b(i) - Ax_ref(i));
  }

  return n;
}

// When scaling the residual by the factor n computed above, it gives a measure
// of the residual taking into account the scale of the solution
double SLS::Scaled_norm(int norm) {
  double n = norm_factor();

  return calc_res_norm(2) / n;
}

// Reset x vector to zero
void SLS::reset_x() { x = Eigen::VectorXd::Zero(mesh.getNumCells()); }

// Reset b vector to zero
void SLS::reset_b() { b = Eigen::VectorXd::Zero(mesh.getNumCells()); }

/*Calculate first-order coefficients for matrix-vector representation
of volume integral of time-derivative for finite volume method.
At a given cell centroid, we take dø/dt ≈ (ø_new - ø_old)/dt.
When Ø denotes the field, we define a matrix A and a vector b such
that dØ/dt ≈ AØ - b. For the above linear approximation, we get
a diagonal marix A for which A_ii = cell_volume/dt and vector b
for which b_i = ø_old * (cell_volume/dt)*/
void SLS::set_temporal_mat(double dt, Field field) {
  int nr_cells = mesh.getNumCells();
  Eigen::ArrayXd old_sol = mesh.get_field(field);

  for (int i = 0; i < nr_cells; i++) {
    A.add_to_entry(i, i, (mesh.get_cell_vol(i) / dt));
    b(i) += (mesh.get_cell_vol(i) / dt) * old_sol(i);
  }
}

/*Calculate first-order coefficients for matrix-vector representation
of volume integral of diffusion term for finite volume method.
At a given face, we take
Sf dot grad(ø)_f = (Delta_f + k_f) dot (grad(ø)_f)
                 = (|Delta_f|/|d_f|)*(ø_N - ø_P) + (k_f dot ((grad(ø))_f),
where we decomposed Sf into a component Delta_f paralell to the vector between
cell centroids (d_f) and one component paralell to the face (non-orthogonal
correction). In very last expression, we interpolate as
grad(ø)_f = f_x * (grad(ø))_O + (1 - f_x) * (grad(ø))_N
where grad(ø))_O and (grad(ø))_N denote gradients at owner
and neighbour-cells, respectively. Using divergence theorem
on each cell, it is now easy to see what terms belong in A and b for
a discretisation of the form
integral_V (div(gamma * grad(Ø)) dV ≈ AØ - b, where gamma is diffusion
coefficient.

**Note that we populate matrices A and b by looping over the FACES of the mesh.
However, the expression AØ - b is evaluated for the cell centroid values ø.

Dirichlet boundary condition:
For a dirichlet boundary condition at a given boundary patch (of which face f is
part) we take Sf dot ((grad(ø))_f) = |Sf|*(ø_b - ø_O)/|d_f|, where ø_b denotes
boundary condition for the specific patch.

Neumann boundary condition:
For a Neumann boundary condition we can simply take Sf dot ((grad(ø))_f) = g_b,
where g_b denotes the boundary condition for the specified patch

Diffusion is Negative of Laplacian. Scale vec contains viscosity in case of
diffusion and coefficient under divergence for laplacian.
*/
void SLS::set_Diffusion_mat(const std::vector<double>& boundary_values,
                            const std::vector<int>& bc_type,
                            const Eigen::ArrayXd& scale_vec, Field field,
                            bool lap) {
  double lc = lap ? -1.0 : 1.0;  // negative coefficients for case of laplacian

  int nr_faces = mesh.getNumFaces();
  int nr_cells = mesh.getNumCells();
  double alpha, f_x;
  Eigen::Vector3d k, PN, Sf, gradf;
  int cell1, cell2, boundary;  // neighbouring cell index and potential index of
                               // boundary patch

  for (int i = 0; i < nr_faces; i++) {
    // Only evaluate diffusion coefficients by evaluating
    // over non-empty faces
    if (!mesh.get_face_empty(i)) {
      cell1 = mesh.get_face_owner(i);
      Sf = mesh.get_face_normal(i);

      boundary = mesh.face_boundary(i);

      double scale1 = scale_vec(cell1);

      // In case face is not on boundary
      if (boundary < 0) {
        cell2 = mesh.get_face_nbr(i);
        double scale2 = scale_vec(cell2);
        double f_x = mesh.get_face_interpolation_factor(i);

        // Compute vector between cell centroids
        // PN = mesh.get_cell_centroid(cell2) - mesh.get_cell_centroid(cell1);

        // Compute correction k_f
        // alpha = Sf.dot(Sf) / PN.dot(Sf);
        // k = Sf - alpha * PN;

        // Compute coefficient |Delta_f|/|d_f|
        double coef = (Sf.norm()) * mesh.get_delta_coef(i);

        double face_scale = f_x * scale1 + (1 - f_x) * scale2;

        // Cell communication coefficients in A can be set directly
        // We subtract the dependence of cells upon themselves
        // as is implied by above diffusion discretisation
        A.add_to_entry(cell1, cell2, -lc * face_scale * coef);
        A.add_to_entry(cell2, cell1, -lc * face_scale * coef);
        A.add_to_entry(cell1, cell1, lc * face_scale * coef);
        A.add_to_entry(cell2, cell2, lc * face_scale * coef);

        // Due to numerical accuracy, set lower bound for |k_f|,
        // below which we force k_f to zero
        /*if (k.norm() > 1e-10) {
          // Get face interpolation factor
          f_x = mesh.get_face_interpolation_factor(i);

          // Account for non-orthogonal correction
          // Only old gradients are involved in this calculation
          // so we subtract from b-vector (source term)
          gradf = f_x * mesh.get_cell_grad(cell1, field) +
                  (1 - f_x) * mesh.get_cell_grad(cell2, field);
          b(cell1) -= k.dot(gradf);
          b(cell2) -= k.dot(gradf);
        }*/
      }

      // Now consider faces at boundary. Note: the boundary patch belonging to
      // a face i is indexed in faces[i].boundary. Hence
      // boundary_values[faces[i].boundary] indicates value

      else if (bc_type[boundary] == 0) {
        // Dirichlet boundary condition as specified above
        A.add_to_entry(cell1, cell1,
                       lc * scale1 * Sf.norm() * mesh.get_delta_coef(i));
        b(cell1) += lc * scale1 * (Sf.norm()) * (mesh.get_delta_coef(i)) *
                    boundary_values[boundary];
      }

      else {
        // Neumann boundary condition as specified above
        b(cell1) += lc * scale1 * boundary_values[boundary] * Sf.norm();
      }
    }
  }
}

/*Calculate first-order coefficients for matrix-vector representation
of volume integral of convection term for finite volume method.
At a given face, we take
(Sf dot F)ø_f = (Sf dot u_f)ø_f,
where u_f and ø_f denote the velocity and intensive variable
interpolated onto the face. We can interpolate ø_f in D ways:
1st method uses central differencing, i.e. ø_f = f_x*ø_P + (1-f_x)*ø_N
2nd methud uses upwind differencing, i.e. ø_f = pos(F)*ø_P + neg(F)*ø_N
We specify the type of interpolation used through up_diff (default is falses)
Using divergence theorem on each cell, it is now easy to see
what terms belong in A and b for a discretisation of the form
integral_V (div(u*Ø))dV = AØ - b

Dirichlet boundary condition:
For a dirichlet boundary condition at a given boundary patch (of which face f
is part) we take Fø_f = F*ø_b, where ø_b denotes boundary condition for the
specific patch.

Neumann boundary condition:
For a Neumann boundary condition we make an interpolation of the form
 ø_b = ø_P + |d_b|g_b, where g_b denotes given boundary gradient at specified
patch and ø_P denotes centroid value of intensive property at boundary cell
*/
void SLS::set_convection_mat(std::vector<double>& boundary_values,
                             std::vector<int>& bc_type, bool up_diff) {
  int nr_faces = mesh.getNumFaces();
  int nr_cells = mesh.getNumCells();
  double alpha, f_x;
  Eigen::Vector3d k, PN, Sf;
  int cell1, cell2, boundary_id;  // neighbouring cell index and potential
                                  // index of boundary patch

  // Loop over faces
  for (int i = 0; i < nr_faces; i++) {
    // only consider faces that are not empty
    if (!mesh.get_face_empty(i)) {
      cell1 = mesh.get_face_owner(i);
      cell2 = mesh.get_face_nbr(i);
      Sf = mesh.get_face_normal(i);

      boundary_id = mesh.face_boundary(i);

      // Take flux as
      double F = mesh.get_face_flux(i);

      if (up_diff == false) {
        f_x = mesh.get_face_interpolation_factor(i);
      } else {
        f_x = (F > 0);
      }

      if (boundary_id < 0) {
        // In the case of upwind differncing. If flux is out of owner cell (F
        // > 0), we want ø_f = ø_Owner. Hence, when flux is out of owner cell,
        // f_x=1.
        A.add_to_entry(cell1, cell1, f_x * F);

        // In the case of upwind differncing. If flux is out of owner cell (F
        // > 0), we want ø_f = ø_Owner. Hence, when flux is out of owner cell,
        // (1-f_x)=0.
        A.add_to_entry(cell1, cell2, (1 - f_x) * F);

        // If flux is out of owner-cell, we do not want ø_f to depend on
        // neighbour cell (cell2). Flux is positive out of owner for F > 0 ->
        // f_x = 1, so (1-f_x)=0. Conversely, if flux is out of nbr cell, then
        // -F > 0
        A.add_to_entry(cell2, cell2, -(1 - f_x) * F);
        A.add_to_entry(cell2, cell1, (-1) * f_x * F);
      }

      // dirichlet
      else if (bc_type[boundary_id] == 0) {
        b(cell1) -= F * boundary_values[boundary_id];
      }

      // neumann
      else {
        A.add_to_entry(cell1, cell1, F);
        b(cell1) -=
            F * (1 / mesh.get_delta_coef(i)) * boundary_values[boundary_id];
      }
    }
  }
}

// Various getter and set functions for the vectors
double SLS::get_x_entry(int i) const { return x[i]; }
double SLS::get_b_entry(int i) const { return b[i]; }
void SLS::set_x_entry(int i, double value) { x[i] = value; }

void SLS::set_b_entry(int i, double value) { b[i] = value; }

void SLS::add_to_b_entry(int i, double value) { b[i] += value; }

void SLS::set_Ax() { Ax = A * x; }
