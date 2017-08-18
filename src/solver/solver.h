#ifndef SOLVER_H_
#define SOLVER_H_

#include <boost/functional/hash.hpp>
#include "../std.h"

class Solver {
 protected:
  size_t n_up;
  size_t n_dn;

  virtual void solve() {}
};

#endif