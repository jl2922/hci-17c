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
  double k_unit;
  double H_unit;
  std::vector<double> rcut_vars;
  std::vector<double> eps_vars;
  std::vector<std::array<int8_t, 3>> k_points;
  std::unordered_map<std::array<int8_t, 3>, std::size_t, boost::hash<std::array<int8_t, 3>>> k_lut;
  std::unordered_map<
      std::array<int8_t, 3>,
      std::vector<std::pair<std::array<int8_t, 3>, double>>,
      boost::hash<std::array<int8_t, 3>>>
      same_spin_hci_queue;  // O(k_points^2).
  std::vector<std::pair<std::array<int8_t, 3>, double>> opposite_spin_hci_queue;  // O(k_points).

  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }

  void setup(const double rcut);

  void generate_hci_queue(const double rcut);

  double hamiltonian(const Det&, const Det&) const override;

  int get_gamma_exp(const SpinDet&, const std::vector<uint16_t>& eor) const;

  std::list<Det> find_connected_dets(const Det&, const double eps) const override;

  std::list<OrbitalPair> get_pq_pairs(const Det&, const Orbital dn_offset) const;
};

#endif