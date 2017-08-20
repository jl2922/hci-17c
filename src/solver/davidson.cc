#include "davidson.h"

#include <Eigen/Dense>

size_t Davidson::diagonalize(
    const std::vector<double>& initial_vector, std::size_t max_iterations) {
  const double TOLERANCE = 1.0e-7;

  if (n == 1) {
    lowest_eigenvalue = diagonal[0];
    lowest_eigenvector = std::vector<double>(1, 1.0);
    diagonalized = true;
    return 0;
  }

  const std::size_t iterations = std::min(n, max_iterations);
  double lowest_eigenvalue = 0.0;
  double lowest_eigenvalue_prev = 0.0;
  double residual_norm = 0.0;

  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(n, iterations);

  if (initial_vector.size() != n) {
    v(0, 0) = 1.0;  // Start from HF.
  } else {
    for (std::size_t i = 0; i < n; i++) v(i, 0) = initial_vector[i];
    v.col(0).normalize();
  }

  Eigen::MatrixXd Hv = Eigen::MatrixXd::Zero(n, iterations);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(n);  // Lowest eigenvector so far.
  Eigen::VectorXd Hw = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd h_krylov = Eigen::MatrixXd::Zero(iterations, iterations);
  Eigen::MatrixXd h_overwrite;
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(iterations);
  std::size_t len_work = 3 * iterations - 1;
  Eigen::VectorXd work(len_work);
  bool converged = false;
  std::vector<double> tmp_v(n);

  // Get diagonal elements.
  Eigen::VectorXd diag_elems(n);
  for (std::size_t i = 0; i < n; i++) diag_elems[i] = diagonal[i];

  // First iteration.
  for (std::size_t i = 0; i < n; i++) tmp_v[i] = v(i, 0);
  const auto& tmp_Hv = apply_hamiltonian(tmp_v);
  for (std::size_t i = 0; i < n; i++) Hv(i, 0) = tmp_Hv[i];
  lowest_eigenvalue = v.col(0).dot(Hv.col(0));
  h_krylov(0, 0) = lowest_eigenvalue;
  w = v.col(0);
  Hw = Hv.col(0);
  if (verbose) printf("Davidson Iteration #1. Eigenvalue: %#.15g\n", lowest_eigenvalue);

  residual_norm = 1.0;  // So at least one iteration is done.
  std::size_t n_iter = std::min(n, iterations);
  size_t n_diagonalize = 1;

  for (std::size_t it = 1; it < n_iter; it++) {
    // Compute residual.
    for (std::size_t j = 0; j < n; j++) {
      v(j, it) = (Hw(j, 0) - lowest_eigenvalue * w(j, 0)) / (lowest_eigenvalue - diag_elems(j));
      if (fabs(lowest_eigenvalue - diag_elems[j]) < 1.0e-8) v(j, it) = -1.0;
    }

    // If residual is small, converge.
    residual_norm = v.col(it).norm();
    if (residual_norm < 1.0e-6) converged = true;

    // Orthogonalize and normalize.
    for (std::size_t i = 0; i < it; i++) {
      double norm = v.col(it).dot(v.col(i));
      v.col(it) -= norm * v.col(i);
    }
    v.col(it).normalize();

    // Apply H once.
    for (std::size_t i = 0; i < n; i++) tmp_v[i] = v(i, it);
    const auto& tmp_Hv2 = apply_hamiltonian(tmp_v);
    for (std::size_t i = 0; i < n; i++) Hv(i, it) = tmp_Hv2[i];

    // Construct Krylow matrix and diagonalize.
    for (std::size_t i = 0; i <= it; i++) {
      h_krylov(i, it) = v.col(i).dot(Hv.col(it));
      h_krylov(it, i) = h_krylov(i, it);
    }

    len_work = 3 * it + 2;
    h_overwrite = h_krylov.leftCols(it + 1).topRows(it + 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(
        h_krylov.leftCols(it + 1).topRows(it + 1));
    const auto& eigenvalues = eigenSolver.eigenvalues();
    const auto& eigenvectors = eigenSolver.eigenvectors();
    lowest_eigenvalue = eigenvalues[0];
    std::size_t lowest_id = 0;
    for (std::size_t i = 1; i <= it; i++) {
      if (eigenvalues[i] < lowest_eigenvalue) {
        lowest_eigenvalue = eigenvalues[i];
        lowest_id = i;
      }
    }
    w = v.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);
    Hw = Hv.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);

    if (it > 1 && fabs(lowest_eigenvalue - lowest_eigenvalue_prev) < TOLERANCE) {
      converged = true;
      break;
    } else {
      lowest_eigenvalue_prev = lowest_eigenvalue;
      n_diagonalize++;
      if (verbose)
        printf("Davidson Iteration #%zu. Eigenvalue: %#.15g\n", n_diagonalize, lowest_eigenvalue);
    }

    if (converged) break;
  }

  this->lowest_eigenvalue = lowest_eigenvalue;
  lowest_eigenvector.resize(n);
  for (std::size_t i = 0; i < n; i++) lowest_eigenvector[i] = w(i);
  diagonalized = true;

  return n_diagonalize;
}