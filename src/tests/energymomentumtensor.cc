/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/energymomentumtensor.h"
#include "../include/smash/fourvector.h"

using namespace smash;

TEST(assign) {
  EnergyMomentumTensor Tmn;
  for (size_t i = 0; i < 10; i++) {
    COMPARE(Tmn[i], 0.0);
  }
  // exact assignment with binary-representable numbers
  EnergyMomentumTensor T(
      {0.25, -37.5, 0.0125, 1.0, 2.0, 3.3, 4.5, 6.9, 10.0, 1.e6});
  COMPARE(T[0], 0.25);
  COMPARE(T[1], -37.5);
  COMPARE(T[2], 0.0125);
  COMPARE(T[9], 1000000.0);
  FUZZY_COMPARE(T[5], 3.3);
}

TEST(arithmetic) {
  EnergyMomentumTensor A(
      {0.25, -37.5, 0.0125, 1.0, 2.0, 3.3, 4.5, 6.9, 10.0, 1.e6});
  EnergyMomentumTensor B(
      {1.25, -36.5, 1.0125, 2.0, 3.0, 4.3, 5.5, 7.9, 11.0, 999999.0});
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
  using se = smash::EnergyMomentumTensor;
  COMPARE(se::tmn_index(0, 0), 0);
  COMPARE(se::tmn_index(0, 1), 1);
  COMPARE(se::tmn_index(0, 2), 2);
  COMPARE(se::tmn_index(0, 3), 3);
  COMPARE(se::tmn_index(1, 1), 4);
  COMPARE(se::tmn_index(1, 2), 5);
  COMPARE(se::tmn_index(1, 3), 6);
  COMPARE(se::tmn_index(2, 2), 7);
  COMPARE(se::tmn_index(2, 3), 8);
  COMPARE(se::tmn_index(3, 3), 9);
  for (std::int8_t i = 0; i < 4; i++) {
    for (std::int8_t j = 0; j < i; j++) {
      COMPARE(se::tmn_index(i, j), se::tmn_index(j, i));
    }
  }
}

TEST_CATCH(invalid_index1, std::invalid_argument) {
  using se = smash::EnergyMomentumTensor;
  se::tmn_index(4, 0);
}

TEST_CATCH(invalid_index2, std::invalid_argument) {
  using se = smash::EnergyMomentumTensor;
  se::tmn_index(0, 4);
}

TEST_CATCH(invalid_index3, std::invalid_argument) {
  using se = smash::EnergyMomentumTensor;
  se::tmn_index(-1, -1);
}

TEST(add_particle) {
  using se = smash::EnergyMomentumTensor;
  EnergyMomentumTensor T;
  const FourVector p = FourVector(1.0, 0.1, 0.2, 0.3);
  T.add_particle(p);
  for (std::int8_t i = 0; i < 4; i++) {
    for (std::int8_t j = 0; j < 4; j++) {
      COMPARE(T[se::tmn_index(i, j)], p[i] * p[j] / p[0]);
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
    COMPARE_RELATIVE_ERROR(u[0] / p1[0], -u[i] / p1[i], 1.e-15);
  }

  // Three particles: check if boost to Landau frame gives T^{0i} = 0
  // By definition of Landau frame T^{0i} = 0 should be fulfilled.
  EnergyMomentumTensor T3L = T3.boosted(T3.landau_frame_4velocity());
  for (size_t i = 1; i < 4; i++) {
    COMPARE_ABSOLUTE_ERROR(T3L[i], 0.0, 2.e-15);
  }
}

TEST(Landau_frame_values) {
  EnergyMomentumTensor T(
      {100., 1.0, 10., 3.3, 30.0, 4.3, 5.5, 29.9, 11.0, 40.0});
  /* Executing this mathematica code:

     T = SetPrecision[{{100., 1.0, 10., 3.3}, {1.0, 30.0, 4.3, 5.5},
          {10.0, 4.3, 29.9, 11.0}, {3.3, 5.5, 11.0, 40.0}}, 20]
     g = SetPrecision[DiagonalMatrix[{1.0, -1.0, -1.0, -1.0}], 20]
     u = Eigenvectors[g.T][[1]]
     u = -u/Sqrt[u.g.u]

     results in u = {1.0030526944248855612, -0.004483908837502199607,
                    -0.07605947885039531450, -0.017594261324820867270}
   */
  const FourVector u = T.landau_frame_4velocity();
  // Allow 8ulp, seen to fail for 3 ulp
  UnitTest::setFuzzyness<double>(8);
  FUZZY_COMPARE(u[0], 1.0030526944248855612);
  FUZZY_COMPARE(u[1], -0.004483908837502199607);
  FUZZY_COMPARE(u[2], -0.07605947885039531450);
  FUZZY_COMPARE(u[3], -0.017594261324820867270);

  /* Executing further mathematica code:

     a = SetPrecision[Flatten[{1.0, u[[2 ;; 4]]/(1.0 + u[[1]])}], 20]
     L = SetPrecision[{u, u[[2]]*a, u[[3]]*a, u[[4]]*a} +
                      DiagonalMatrix[{0.0, 1.0, 1.0, 1.0}], 20]
     L.T.L

     results in the values I compare */
  EnergyMomentumTensor TL = T.boosted(u);
  FUZZY_COMPARE(TL[0], 99.17936538699464299);
  COMPARE_ABSOLUTE_ERROR(TL[1], 0.0, 1.e-14) << TL[1];
  COMPARE_ABSOLUTE_ERROR(TL[2], 0.0, 1.e-14) << TL[2];
  COMPARE_ABSOLUTE_ERROR(TL[3], 0.0, 1.e-14) << TL[3];
  FUZZY_COMPARE(TL[4], 29.995527036975529013);
  FUZZY_COMPARE(TL[5], 4.239712597343882893);
  FUZZY_COMPARE(TL[6], 5.483845238049431987);
  FUZZY_COMPARE(TL[7], 29.141747611031320028);
  FUZZY_COMPARE(TL[8], 10.787129594442447275);
  FUZZY_COMPARE(TL[9], 39.94209073898776673);
}
