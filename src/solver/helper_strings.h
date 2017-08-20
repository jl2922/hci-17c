#ifndef HCI_HELPER_STRINGS_H_
#define HCI_HELPER_STRINGS_H_

#include <boost/functional/hash.hpp>
#include "../std.h"
#include "../wavefunction/wavefunction.h"

class HelperStrings {
 public:
  HelperStrings(std::function<double(Det, Det)>& hamiltonian) : hamiltonian(hamiltonian) {}

  void setup(const std::vector<Det>& dets);

  std::vector<std::pair<size_t, double>> find_connections(const std::size_t i);

 private:
  std::vector<Det> dets;

  std::vector<std::vector<std::pair<size_t, double>>> cached_connections;

  std::vector<bool> cached;

  size_t cache_size;

  std::function<double(Det, Det)>& hamiltonian;

  // alpha and beta strings, O(n_dets).
  std::unordered_map<
      Orbitals,
      std::pair<std::vector<size_t>, std::vector<size_t>>,
      boost::hash<Orbitals>>
      ab;

  // alpha-m1 and beta-m1 strings, O(n_dets * n_elecs).
  std::unordered_map<
      Orbitals,
      std::pair<std::vector<size_t>, std::vector<size_t>>,
      boost::hash<Orbitals>>
      ab_m1;

  // Whether has been included in the potential connections.
  std::vector<std::vector<bool>> connected;

  // Whether the variational dets are one-up excitations of the det passed in.
  std::vector<std::vector<bool>> one_up;

  // Setup alpha and beta.
  void setup_ab();

  // Setup alpha-m1 and beta-m1.
  void setup_ab_m1();
};

#endif