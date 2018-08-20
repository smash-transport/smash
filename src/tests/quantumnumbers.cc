/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "../include/smash/quantumnumbers.h"

using namespace smash;

TEST(init_particle_types) {
  ParticleType::create_type_list(
      "σ 0.123 -.0 + 123\n"
      "σ3 0.123 -.0 + -1234568\n"
      "σ2 0.245 2.3 + 2346\n");
}

TEST(size) { COMPARE(sizeof(QuantumNumbers), 56u); }

TEST(assign_empty) {
  QuantumNumbers emptyset;
  COMPARE(emptyset.momentum().x0(), 0);
  COMPARE(emptyset.momentum().x1(), 0);
  COMPARE(emptyset.momentum().x2(), 0);
  COMPARE(emptyset.momentum().x3(), 0);
  COMPARE(emptyset.charge(), 0);
  COMPARE(emptyset.isospin3(), 0);
  COMPARE(emptyset.strangeness(), 0);
  COMPARE(emptyset.charmness(), 0);
  COMPARE(emptyset.bottomness(), 0);
  COMPARE(emptyset.baryon_number(), 0);
}

TEST(assign_full) {
  FourVector P(1, 2, 3, 4);
  QuantumNumbers fullset(P, 5, 6, 7, 8, 9, 0);
  COMPARE(fullset.momentum(), FourVector(1, 2, 3, 4));
  COMPARE(fullset.charge(), 5);
  COMPARE(fullset.isospin3(), 6);
  COMPARE(fullset.strangeness(), 7);
  COMPARE(fullset.charmness(), 8);
  COMPARE(fullset.bottomness(), 9);
  COMPARE(fullset.baryon_number(), 0);
}

TEST(compare) {
  FourVector P(1, 2, 3, 4);
  FourVector Q(2, 2, 3, 4);
  QuantumNumbers A(P, 5, 6, 7, 8, 9, 0);
  QuantumNumbers B = A;
  QuantumNumbers C(Q, 5, 6, 7, 8, 9, 0);
  QuantumNumbers D(Q, 1, 6, 7, 8, 9, 0);
  QuantumNumbers E(P, 5, 1, 7, 8, 9, 0);
  QuantumNumbers F(P, 5, 1, 3, 8, 9, 0);
  QuantumNumbers G(Q, 5, 6, 7, -8, 9, 0);
  QuantumNumbers H(Q, 5, 6, 7, -8, 12358, 0);
  QuantumNumbers I(Q, 5, 6, 7, 8, 9, -78);
  VERIFY(A == B);
  // these differ in FourVector:
  VERIFY(A != C);
  // these differ in Charge:
  VERIFY(C != D);
  // these differ in Isospin:
  VERIFY(A != E);
  // these differ in Strangeness:
  VERIFY(E != F);
  // these differ in Charmness:
  VERIFY(D != G);
  // these differ in Bottomness:
  VERIFY(G != H);
  // these differ in Baryon Number:
  VERIFY(C != I);
  // these all differ in more than one entry (basically, these pairs
  // have been chosen randomly):
  VERIFY(A != D);
  VERIFY(A != E);
  VERIFY(A != F);
  VERIFY(D != I);
}

TEST(difference) {
  FourVector P(1, 2, 3, 4);
  FourVector Q(2, 3, 4, 4);
  QuantumNumbers A(P, 5, 6, 7, 8, 9, 0);
  QuantumNumbers H(Q, 5, 6, 1, -8, 12358, -15);
  QuantumNumbers diff(P - Q, 0, 0, 6, 16, -12349, 15);
  COMPARE(diff, A - H);
}

TEST(report_deviations) {
  FourVector P(1, 2, 3, 4);
  FourVector Q(2, 3, 4, 4);
  // FourVector Q(1,2,4,4);
  QuantumNumbers A(P, 5, 6, 7, 8, 9, 0);
  QuantumNumbers H(Q, 5, 6, 1, -8, 12358, -15);
  COMPARE(A.report_deviations(H),
          "Conservation law violations detected (old vs. new)\n"
          "Deviation in Four-Momentum:\n"
          " P_0: 1.000000e+00 vs. 2.000000e+00; Δ = -1.000000e+00\n"
          " P_1: 2.000000e+00 vs. 3.000000e+00; Δ = -1.000000e+00\n"
          " P_2: 3.000000e+00 vs. 4.000000e+00; Δ = -1.000000e+00\n"
          "Deviation in Strangeness:\n"
          " 7 vs. 1\n"
          "Deviation in Charmness:\n"
          " 8 vs. -8\n"
          "Deviation in Bottomness:\n"
          " 9 vs. 12358\n"
          "Deviation in Baryon Number:\n"
          " 0 vs. -15\n");
  // small deviation in FourVector should satisfy ==:
  FourVector R(2 + 2e-13, 3, 4, 4);
  QuantumNumbers J(R, 5, 6, 1, -8, 12358, -15);
  COMPARE(H.report_deviations(J), "");
}

TEST(count_from_particles) {
  // we will successively fill a particle list and compare the quantum
  // numbers with it. Finally, we will see if report_deviations also
  // works for particle lists.

  // create particle with fake PdgCode (this would be equivalent to a
  // rho^+, but that should be "213").
  ParticleData particle(ParticleType::find(PdgCode("123")));
  FourVector P(1, 2, 3, 4);
  particle.set_4momentum(P);
  // create particle list:
  Particles list;
  list.insert(particle);

  QuantumNumbers onlyone(list);
  QuantumNumbers check1(P, 1, 2, 0, 0, 0, 0);
  // won't print anything if the following VERIFY succeeds
  std::printf("%s", check1.report_deviations(onlyone).c_str());
  COMPARE(onlyone, check1);

  ParticleData particleQ(ParticleType::find(PdgCode("123")));
  FourVector Q(2, 3, 4, 5);
  particleQ.set_4momentum(Q);
  list.insert(particleQ);

  QuantumNumbers two(list);
  QuantumNumbers check2(P + Q, 2, 4, 0, 0, 0, 0);
  std::printf("%s", check2.report_deviations(two).c_str());
  COMPARE(two, check2);

  ParticleData particleR(ParticleType::find(PdgCode("2346")));
  FourVector R(3, 4, 5, 6);
  particleR.set_4momentum(R);
  list.insert(particleR);

  QuantumNumbers three(list);
  QuantumNumbers check3(P + Q + R, 3, 5, -1, 1, 0, 1);
  std::printf("%s", check3.report_deviations(three).c_str());
  COMPARE(three, check3);

  ParticleData particleS(ParticleType::find(PdgCode("-1234568")));
  FourVector S(-6, -9, -12, -15);
  particleS.set_4momentum(S);
  list.insert(particleS);

  QuantumNumbers four(list);
  QuantumNumbers check4(P + Q + R + S, 2, 5, -1, 0, 1, 0);
  std::printf("%s", check4.report_deviations(four).c_str());
  COMPARE(four, check4);

  //
  COMPARE(three.report_deviations(list),
          "Conservation law violations detected (old vs. new)\n"
          "Deviation in Four-Momentum:\n"
          " P_0: 6.000000e+00 vs. 0.000000e+00; Δ = 6.000000e+00\n"
          " P_1: 9.000000e+00 vs. 0.000000e+00; Δ = 9.000000e+00\n"
          " P_2: 1.200000e+01 vs. 0.000000e+00; Δ = 1.200000e+01\n"
          " P_3: 1.500000e+01 vs. 0.000000e+00; Δ = 1.500000e+01\n"
          "Deviation in Charge:\n"
          " 3 vs. 2\n"
          "Deviation in Charmness:\n"
          " 1 vs. 0\n"
          "Deviation in Bottomness:\n"
          " 0 vs. 1\n"
          "Deviation in Baryon Number:\n"
          " 1 vs. 0\n");
}
