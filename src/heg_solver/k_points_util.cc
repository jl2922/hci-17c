#include "k_points_util.h"

#include <boost/functional/hash.hpp>
#include "../array_math.h"

size_t KPointsUtil::get_n_k_points(const double rcut) {
  size_t count = 0;
  const int8_t n_max = floor(rcut);
  for (int8_t i = -n_max; i <= n_max; i++) {
    for (int8_t j = -n_max; j <= n_max; j++) {
      for (int8_t k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > pow(rcut, 2)) continue;
        count++;
      }
    }
  }
  return count;
}

std::vector<std::array<int8_t, 3>> KPointsUtil::generate_k_points(const double rcut) {
  std::vector<std::array<int8_t, 3>> k_points;
  const int8_t n_max = floor(rcut);
  for (int8_t i = -n_max; i <= n_max; i++) {
    for (int8_t j = -n_max; j <= n_max; j++) {
      for (int8_t k = -n_max; k <= n_max; k++) {
        if (i * i + j * j + k * k > pow(rcut, 2)) continue;
        k_points.push_back(std::array<int8_t, 3>({i, j, k}));
      }
    }
  }
  std::stable_sort(
      k_points.begin(),
      k_points.end(),
      [](const std::array<int8_t, 3>& a, const std::array<int8_t, 3>& b) -> bool {
        return squared_norm(a) < squared_norm(b);
      });
  return k_points;
}

std::vector<std::array<int8_t, 3>> KPointsUtil::get_k_diffs(
    const std::vector<std::array<int8_t, 3>>& k_points) {
  // Generate all possible differences between two different k points.
  std::unordered_set<std::array<int8_t, 3>, boost::hash<std::array<int8_t, 3>>> k_diffs_set;
  std::vector<std::array<int8_t, 3>> k_diffs;
  const size_t n_orbs = k_points.size();
  for (size_t p = 0; p < n_orbs; p++) {
    for (size_t q = 0; q < n_orbs; q++) {
      if (p == q) continue;
      const auto& diff_pq = k_points[q] - k_points[p];
      if (k_diffs_set.count(diff_pq) == 1) continue;
      k_diffs.push_back(diff_pq);
      k_diffs_set.insert(diff_pq);
    }
  }

  // Sort k_diffs into ascending order so that later sorting hci queue will be faster.
  std::stable_sort(
      k_diffs.begin(),
      k_diffs.end(),
      [](const std::array<int8_t, 3>& a, const std::array<int8_t, 3>& b) -> bool {
        return squared_norm(a) < squared_norm(b);
      });

  return k_diffs;
}