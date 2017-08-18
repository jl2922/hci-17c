#include "heg_solver.h"

#include <boost/functional/hash.hpp>
#include "../config.h"
#include "../parallel.h"
#include "../time.h"
#include "omp.h"

void HEGSolver::solve() {
  n_up = Config::get<size_t>("n_up");
  n_dn = Config::get<size_t>("n_dn");
  rcut_vars = Config::get_array<double>("rcut_vars");
  eps_vars = Config::get_array<double>("eps_vars");
}