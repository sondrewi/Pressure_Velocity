#include "GSSolver.H"

/*Gauss-Seidel Solver class. The solver holds references to a smoother object
that performs the actual sweeps of the Gauss-Seidel algorithm and a sparse
linear system instance sls The solver keeps track of the residual in the
associated linear system and stops calling the smoother object when this is
sufficiently low*/

// Constructor for Gauss-Seidel solver. Specify number of sweeps between each
// time that residual is evaluated and max number of iterations
GSSolver::GSSolver(GSSmooth& smoother_object, int sweeps_, int maxit)
    : sls(smoother_object.sls), smoother(smoother_object) {
  sweeps = sweeps_;
  max_it = maxit;
}

double GSSolver::solve(double relTol) {
  int iter_count = 0;
  double prev_resid_norm = sls.calc_res_norm(2);
  double b_norm = sls.calc_b_norm(2);

  // While relative tolerance time b-norm is smaller than residual norm or max
  // iterations not exceeded, carry out a new sets of sweeps of Gauss-Seidel
  // method
  while (b_norm * relTol < sls.calc_res_norm(2) && iter_count < max_it) {
    iter_count += 1;

    smoother.smooth(sweeps);

    // Set new residual and Ax vector based on current guess x
    sls.set_Ax();
    sls.set_res();
    prev_resid_norm = sls.calc_res_norm(2);
  }
  return prev_resid_norm;
}
