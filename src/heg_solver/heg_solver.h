#ifndef HEG_SOLVER_H_
#define HEG_SOLVER_H_

#include <boost/functional/hash.hpp>
#include "../solver/solver.h"
#include "../std.h"

class HEGSolver : public Solver {
 public:
  static void run() { HEGSolver::get_instance().solve(); }

  void solve() override;

 private:
  std::vector<double> rcut_vars;
  std::vector<double> eps_vars;

  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }
};

#endif