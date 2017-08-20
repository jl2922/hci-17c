#include "spin_det.h"
#include "gtest/gtest.h"
#include "types.h"

TEST(SpinDetTest, SetAndGetOrbitals) {
  SpinDet spin_det;
  EXPECT_FALSE(spin_det.get_orb(5));
  spin_det.set_orb(5, true);
  EXPECT_TRUE(spin_det.get_orb(5));
  EXPECT_EQ(spin_det.get_n_elecs(), 1);
  EXPECT_EQ(spin_det.get_elec_orbs().size(), 1);
  EXPECT_EQ(spin_det.get_elec_orbs()[0], 5);
}

TEST(SpinDetTest, EncodeAndDecodeFixed) {
  SpinDet spin_det1;
  spin_det1.set_orb(2, true);
  spin_det1.set_orb(3, true);
  SpinDet spin_det2;
  spin_det2.decode(spin_det1.encode(SpinDet::EncodeScheme::FIXED), SpinDet::EncodeScheme::FIXED);
  EXPECT_TRUE(spin_det1 == spin_det2);
}

TEST(SpinDetTest, EncodeAndDecodeVariable) {
  SpinDet spin_det1;
  spin_det1.set_orb(2, true);
  spin_det1.set_orb(3, true);
  SpinDet spin_det2;
  Orbitals code = spin_det1.encode(SpinDet::EncodeScheme::VARIABLE);
  spin_det2.decode(code, SpinDet::EncodeScheme::VARIABLE);
  Orbitals expected_code({2, 0, 1, 2, 3});
  for (size_t i = 0; i < expected_code.size(); i++) {
    EXPECT_EQ(code[i], expected_code[i]);
  }
  EXPECT_TRUE(spin_det1 == spin_det2);

  spin_det1.set_orb(1, true);
  code = spin_det1.encode(SpinDet::EncodeScheme::VARIABLE);
  spin_det2.decode(code, SpinDet::EncodeScheme::VARIABLE);
  expected_code = Orbitals({3, 0, 3});
  for (size_t i = 0; i < expected_code.size(); i++) {
    EXPECT_EQ(code[i], expected_code[i]);
  }
  EXPECT_TRUE(spin_det1 == spin_det2);

  spin_det1.set_orb(0, true);
  code = spin_det1.encode(SpinDet::EncodeScheme::VARIABLE);
  spin_det2.decode(code, SpinDet::EncodeScheme::VARIABLE);
  expected_code = Orbitals({4});
  for (size_t i = 0; i < expected_code.size(); i++) {
    EXPECT_EQ(code[i], expected_code[i]);
  }
  EXPECT_TRUE(spin_det1 == spin_det2);
}

TEST(SpinDetTest, FromEOR) {
  SpinDet spin_det1, spin_det2, spin_det3;
  spin_det1.set_orb(1, true);
  spin_det1.set_orb(2, true);
  spin_det2.set_orb(2, true);
  spin_det2.set_orb(3, true);
  spin_det3.from_eor(spin_det1, spin_det2);
  EXPECT_EQ(spin_det3.get_n_elecs(), 2);
  EXPECT_TRUE(spin_det3.get_orb(1));
  EXPECT_TRUE(spin_det3.get_orb(3));
}