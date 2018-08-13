/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/clebschgordan.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

TEST(coefficient) {
  /* spins are two times the actual values,
   * so j = 1 for spin-1/2 particle, 2 for spin-1 particle, etc.
   * Ordering of spins in array:
   *  0: j1, 1: j2, 2: j3, 3: m1, 4: m2, 5: m3
   */
  int spin[7][3];
  int spinz[7][3];
  double correct_coefficient[7];

  spin[0][0] = 1;
  spin[0][1] = 1;
  spin[0][2] = 2;
  spinz[0][0] = 1;
  spinz[0][1] = 1;
  spinz[0][2] = 2;
  correct_coefficient[0] = 1.;

  spin[1][0] = 1;
  spin[1][1] = 1;
  spin[1][2] = 2;
  spinz[1][0] = 1;
  spinz[1][1] = -1;
  spinz[1][2] = 0;
  correct_coefficient[1] = 1 / sqrt(2.);

  spin[2][0] = 2;
  spin[2][1] = 1;
  spin[2][2] = 1;
  spinz[2][0] = 2;
  spinz[2][1] = -1;
  spinz[2][2] = 1;
  correct_coefficient[2] = sqrt(2. / 3.);

  spin[3][0] = 2;
  spin[3][1] = 1;
  spin[3][2] = 3;
  spinz[3][0] = -2;
  spinz[3][1] = 1;
  spinz[3][2] = -1;
  correct_coefficient[3] = sqrt(1. / 3.);

  spin[4][0] = 2;
  spin[4][1] = 2;
  spin[4][2] = 2;
  spinz[4][0] = 0;
  spinz[4][1] = 2;
  spinz[4][2] = 2;
  correct_coefficient[4] = -1 / sqrt(2.);

  spin[5][0] = 2;
  spin[5][1] = 2;
  spin[5][2] = 2;
  spinz[5][0] = 0;
  spinz[5][1] = 0;
  spinz[5][2] = 0;
  correct_coefficient[5] = 0.;

  spin[6][0] = 2;
  spin[6][1] = 2;
  spin[6][2] = 4;
  spinz[6][0] = 2;
  spinz[6][1] = -2;
  spinz[6][2] = 0;
  correct_coefficient[6] = 1 / sqrt(6.);
  for (int i = 0; i < 7; i++) {
    double cg = clebsch_gordan(spin[i][0], spin[i][1], spin[i][2], spinz[i][0],
                               spinz[i][1], spinz[i][2]);
    FUZZY_COMPARE(cg, correct_coefficient[i])
        << '\n'  // Using double quotes here produces an error(?!)
        << "J1: " << spin[i][0] << " Jz1: " << spinz[i][0] << "\n"
        << "J2: " << spin[i][1] << " Jz2: " << spinz[i][1] << "\n"
        << "J3: " << spin[i][2] << " Jz3: " << spinz[i][2] << "\n"
        << "CG: " << cg << " Correct: " << correct_coefficient[i];
  }
}

const double tolerance = 1.0e-7;

TEST(iso_clebsch_2to1) {
  const ParticleType &pip = ParticleType::find(0x211);
  const ParticleType &piz = ParticleType::find(0x111);
  const ParticleType &pim = ParticleType::find(-0x211);
  const ParticleType &rho_p = ParticleType::find(0x213);
  const ParticleType &rho_z = ParticleType::find(0x113);
  const ParticleType &rho_m = ParticleType::find(-0x213);
  const ParticleType &sigma = ParticleType::find(0x9000221);
  const ParticleType &proton = ParticleType::find(0x2212);
  const ParticleType &neutron = ParticleType::find(0x2112);
  const ParticleType &Delta_pp = ParticleType::find(0x2224);
  const ParticleType &Delta_p = ParticleType::find(0x2214);
  const ParticleType &Delta_z = ParticleType::find(0x2114);
  const ParticleType &Delta_m = ParticleType::find(0x1114);

  // π π -> X
  double iso_cg = isospin_clebsch_gordan_sqr_2to1(pip, piz, rho_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 2., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, pip, rho_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 2., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to1(pip, pim, rho_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 2., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pim, pip, rho_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 2., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, piz, rho_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pip, pim, sigma);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pim, pip, sigma);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, piz, sigma);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, pim, rho_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 2., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pim, piz, rho_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 2., tolerance);

  // π N -> X
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pip, proton, Delta_pp);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to1(pip, neutron, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pip, neutron, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 2. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, proton, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 2. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, proton, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, neutron, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 2. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(piz, neutron, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pim, proton, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 3., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to1(pim, proton, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 2. / 3., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to1(pim, neutron, Delta_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1., tolerance);
}

TEST(iso_clebsch_3to1) {
  const ParticleType &pip = ParticleType::find(0x211);
  const ParticleType &piz = ParticleType::find(0x111);
  const ParticleType &pim = ParticleType::find(-0x211);
  const ParticleType &omega = ParticleType::find(0x223);

  // ω -> π⁺π⁻π⁰ : all permutations are equally likely
  double iso_cg = isospin_clebsch_gordan_sqr_3to1(pip, piz, pim, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 6., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_3to1(pip, pim, piz, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 6., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_3to1(piz, pip, pim, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 6., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_3to1(piz, pim, pip, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 6., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_3to1(pim, pip, piz, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 6., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_3to1(pim, piz, pip, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 6., tolerance);

  // ω -> 3π⁰ : forbidden
  iso_cg = isospin_clebsch_gordan_sqr_3to1(piz, piz, piz, omega);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);
}

TEST(iso_clebsch_2to2) {
  const ParticleType &proton = ParticleType::find(0x2212);
  const ParticleType &neutron = ParticleType::find(0x2112);
  const ParticleType &Delta_pp = ParticleType::find(0x2224);
  const ParticleType &Delta_p = ParticleType::find(0x2214);
  const ParticleType &Delta_z = ParticleType::find(0x2114);
  const ParticleType &Delta_m = ParticleType::find(0x1114);

  // N N -> N N
  double iso_cg =
      isospin_clebsch_gordan_sqr_2to2(proton, proton, proton, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, proton, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.5, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, neutron, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.5, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(neutron, neutron, neutron, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1., tolerance);

  // N N -> N N with given total isospsin
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, proton, proton, 2);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, proton, proton, 0);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, proton, neutron, 2);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.25, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, proton, neutron, 0);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.25, tolerance);
  iso_cg =
      isospin_clebsch_gordan_sqr_2to2(neutron, neutron, neutron, neutron, 2);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1., tolerance);
  iso_cg =
      isospin_clebsch_gordan_sqr_2to2(neutron, neutron, neutron, neutron, 0);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);

  // N N -> N Delta
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_pp, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3. / 4., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_p, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 4., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_p, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 4., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_z, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1. / 4., tolerance);

  // N N -> Delta Delta
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_pp, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3. / 10., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_z, Delta_pp);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3. / 10., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_p, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 4. / 10., tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_pp, Delta_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 7. / 20., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_m, Delta_pp);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 7. / 20., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_p, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3. / 20., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_z, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3. / 20., tolerance);

  // invalid cases
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, proton, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_z, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_m, Delta_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0., tolerance);
}

TEST(range) {
  const auto &proton = ParticleType::find(0x2212);
  const auto &neutron = ParticleType::find(0x2112);
  const auto &Delta_pp = ParticleType::find(0x2224);
  const auto &Delta_z = ParticleType::find(0x2114);
  const auto &Lambda = ParticleType::find(0x3122);
  const auto &Lambda1520 = ParticleType::find(0x3124);
  const auto &Sigma1385_z = ParticleType::find(0x3214);

  {
    auto range = I_tot_range(Delta_z, Delta_pp);
    std::vector<int> isospins_from_range;
    std::copy(range.begin(), range.end(),
              std::back_inserter(isospins_from_range));
    std::vector<int> isospins_expected = {6, 4, 2};
    COMPARE(isospins_from_range.size(), isospins_expected.size());
    for (size_t i = 0; i < isospins_expected.size(); i++) {
      COMPARE(isospins_from_range[i], isospins_expected[i]);
    }
  }

  {
    auto range = I_tot_range(proton, Lambda);
    std::vector<int> isospins_from_range;
    std::copy(range.begin(), range.end(),
              std::back_inserter(isospins_from_range));
    std::vector<int> isospins_expected = {1};
    COMPARE(isospins_from_range.size(), isospins_expected.size());
    for (size_t i = 0; i < isospins_expected.size(); i++) {
      COMPARE(isospins_from_range[i], isospins_expected[i]);
    }
  }

  {
    auto range = I_tot_range(proton, proton, Lambda1520, Delta_pp);
    std::vector<int> isospins_from_range;
    std::copy(range.begin(), range.end(),
              std::back_inserter(isospins_from_range));
    std::vector<int> isospins_expected = {};
    COMPARE(isospins_from_range.size(), isospins_expected.size());
    for (size_t i = 0; i < isospins_expected.size(); i++) {
      COMPARE(isospins_from_range[i], isospins_expected[i]);
    }
  }

  {
    auto range = I_tot_range(proton, neutron, Sigma1385_z, proton);
    std::vector<int> isospins_from_range;
    std::copy(range.begin(), range.end(),
              std::back_inserter(isospins_from_range));
    std::vector<int> isospins_expected = {};
    COMPARE(isospins_from_range.size(), isospins_expected.size());
    for (size_t i = 0; i < isospins_expected.size(); i++) {
      COMPARE(isospins_from_range[i], isospins_expected[i]);
    }
  }
}
