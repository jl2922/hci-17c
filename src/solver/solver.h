#ifndef SOLVER_H_
#define SOLVER_H_

#include <boost/functional/hash.hpp>
#include "../std.h"
#include "../wavefunction/wavefunction.h"
#include "helper_strings.h"

class Solver {
 protected:
  size_t n_up;
  size_t n_dn;
  double max_abs_H;
  double energy_hf;
  double energy_var;

  Wavefunction wf;

  virtual void solve() {}

  virtual double hamiltonian(const Det&, const Det&) const = 0;

  void variation(const double eps);

  virtual std::list<Det> find_connected_dets(const Det&, const double eps) const = 0;

 private:
  bool converged;

  Det generate_hf_det();

  double diagonalize(const bool has_new_dets);

  std::vector<double> apply_hamiltonian(const std::vector<double>&, HelperStrings&);
};

#endif