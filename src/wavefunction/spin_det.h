#ifndef SPIN_DET_H_
#define SPIN_DET_H_

#include "../std.h"
#include "types.h"

class SpinDet {
 public:
  enum EncodeScheme { FIXED, VARIABLE };

  void set_orb(const Orbital orb_id, const bool occ);

  bool get_orb(const Orbital orb_id) const {
    return std::binary_search(elecs.begin(), elecs.end(), orb_id);
  }

  size_t get_n_elecs() const { return elecs.size(); }

  void from_eor(const SpinDet&, const SpinDet&);

  const Orbitals get_elec_orbs() const { return elecs; }

  const Orbitals encode(const EncodeScheme scheme = VARIABLE) const {
    if (scheme == FIXED) return elecs;
    return encode_variable();
  }

  void decode(const Orbitals& code, const EncodeScheme scheme = VARIABLE) {
    if (scheme == FIXED) {
      elecs = code;
    } else {
      decode_variable(code);
    }
  }

  friend bool operator==(const SpinDet&, const SpinDet&);

  friend bool operator!=(const SpinDet&, const SpinDet&);

  friend std::ostream& operator<<(std::ostream&, const SpinDet&);

 private:
  Orbitals elecs;

  const Orbitals encode_variable() const;

  void decode_variable(const Orbitals& code);
};

bool operator==(const SpinDet&, const SpinDet&);

bool operator!=(const SpinDet&, const SpinDet&);

std::ostream& operator<<(std::ostream&, const SpinDet&);

#endif