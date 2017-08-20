#ifndef K_POINTS_UTIL_H_
#define K_POINTS_UTIL_H_

#include "../std.h"

class KPointsUtil {
 public:
  static size_t get_n_k_points(const double);

  static std::vector<std::array<int8_t, 3>> generate_k_points(const double);

  static std::vector<std::array<int8_t, 3>> get_k_diffs(const std::vector<std::array<int8_t, 3>>&);
};
#endif