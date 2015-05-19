/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "../include/energymomentumtensor.h"
#include "../include/fourvector.h"
#include "unittest.h"

using namespace Smash;

TEST(assign) {
  EnergyMomentumTensor Tmn;
  for (size_t i = 0; i < 10; i++) {
    COMPARE(Tmn[i], 0.0);
  }
  // exact assignment with binary-representable numbers
  EnergyMomentumTensor T({0.25, -37.5, 0.0125, 1.0, 2.0, 3.3, 4.5, 6.9, 10.0, 1.e6});
  COMPARE(T[0], 0.25);
  COMPARE(T[1], -37.5);
  COMPARE(T[2], 0.0125);
  COMPARE(T[9], 1000000.0);
  FUZZY_COMPARE(T[5], 3.3);
}

TEST(arithmetic) {
  EnergyMomentumTensor A({0.25, -37.5, 0.0125, 1.0, 2.0, 3.3, 4.5, 6.9, 10.0, 1.e6});
  EnergyMomentumTensor B({1.25, -36.5, 1.0125, 2.0, 3.0, 4.3, 5.5, 7.9, 11.0, 999999.0});
  EnergyMomentumTensor C = A - B;
  FUZZY_COMPARE(C[0], -1.0);
  FUZZY_COMPARE(C[1], -1.0);
  FUZZY_COMPARE(C[9], 1.0);
  C += B;
  FUZZY_COMPARE(C[0], A[0]);
  FUZZY_COMPARE(C[9], A[9]);
  FUZZY_COMPARE(C[5], A[5]);
  C -= A;
  FUZZY_COMPARE(C[4], 0.0);
  FUZZY_COMPARE(C[3], 0.0);
  FUZZY_COMPARE(C[7], 0.0);
  EnergyMomentumTensor D = A;
  D *= 2;
  FUZZY_COMPARE(D[2], 2.0 * A[2]);
  FUZZY_COMPARE(D[8], 2.0 * A[8]);
  FUZZY_COMPARE(D[9], 2.0 * A[9]);
  D /= 2;
  FUZZY_COMPARE(D[3], A[3]);
  FUZZY_COMPARE(D[7], A[7]);
  FUZZY_COMPARE(D[9], A[9]);
}

TEST(indices) {
  using se = Smash::EnergyMomentumTensor;
  VERIFY(se::tmn_index(0,0) == 0);
  VERIFY(se::tmn_index(0,1) == 1);
  VERIFY(se::tmn_index(0,2) == 2);
  VERIFY(se::tmn_index(0,3) == 3);
  VERIFY(se::tmn_index(1,1) == 4);
  VERIFY(se::tmn_index(1,2) == 5);
  VERIFY(se::tmn_index(1,3) == 6);
  VERIFY(se::tmn_index(2,2) == 7);
  VERIFY(se::tmn_index(2,3) == 8);
  VERIFY(se::tmn_index(3,3) == 9);
  for (std::int8_t i = 0; i < 4; i++) {
    for (std::int8_t j = 0; j < i; j++) {
      VERIFY(se::tmn_index(i,j) == se::tmn_index(j,i));
    }
  }
}

TEST(Landau_frame) {
  const FourVector p1 = FourVector(1.0, 0.1, 0.2, 0.3);
  EnergyMomentumTensor T1, T3;
  T1.add_particle(p1);
  T3.add_particle(p1);
  T3.add_particle(FourVector(2.0, 0.3, 0.1, 0.4));
  T3.add_particle(FourVector(3.0, 1.3, 0.3, 0.7));
  ThreeVector v;

  // One particle: is u ~ momentum?
  FourVector u = T1.landau_frame_4velocity();
  for (size_t i = 1; i < 4; i++) {
    COMPARE_RELATIVE_ERROR(u[0]/p1[0], -u[i]/p1[i], 1.e-15);
  }

  // Three particles: check if boost to Landau frame gives T^{0i} = 0
  // By definition of Landau frame T^{0i} = 0 should be fulfilled.
  EnergyMomentumTensor T3L = T3.boosted(T3.landau_frame_4velocity());
  for (size_t i = 1; i < 4; i++) {
    COMPARE_ABSOLUTE_ERROR(T3L[i], 0.0, 1.e-15);
  }
}
