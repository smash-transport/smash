/*
 *
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include "../include/clebschgordan.h"

using namespace Smash;

TEST(init_particle_types) {
  Test::create_actual_particletypes();
}


TEST(coefficient) {
  /* spins are two times the actual values,
   * so j = 1 for spin-1/2 particle, 2 for spin-1 particle, etc.
   * Ordering of spins in array:
   *  0: j1, 1: j2, 2: j3, 3: m1, 4: m2, 5: m3
   */
  int spin[7][3];
  int spinz[7][3];
  float correct_coefficient[7];

  spin[0][0] = 1;
  spin[0][1] = 1;
  spin[0][2] = 2;
  spinz[0][0] = 1;
  spinz[0][1] = 1;
  spinz[0][2] = 2;
  correct_coefficient[0] = 1.0;

  spin[1][0] = 1;
  spin[1][1] = 1;
  spin[1][2] = 2;
  spinz[1][0] = 1;
  spinz[1][1] = -1;
  spinz[1][2] = 0;
  correct_coefficient[1] = 1 / sqrt(2.0);

  spin[2][0] = 2;
  spin[2][1] = 1;
  spin[2][2] = 1;
  spinz[2][0] = 2;
  spinz[2][1] = -1;
  spinz[2][2] = 1;
  correct_coefficient[2] = sqrt(2.0 / 3.0);

  spin[3][0] = 2;
  spin[3][1] = 1;
  spin[3][2] = 3;
  spinz[3][0] = -2;
  spinz[3][1] = 1;
  spinz[3][2] = -1;
  correct_coefficient[3] = sqrt(1.0 / 3.0);

  spin[4][0] = 2;
  spin[4][1] = 2;
  spin[4][2] = 2;
  spinz[4][0] = 0;
  spinz[4][1] = 2;
  spinz[4][2] = 2;
  correct_coefficient[4] = -1 / sqrt(2.0);

  spin[5][0] = 2;
  spin[5][1] = 2;
  spin[5][2] = 2;
  spinz[5][0] = 0;
  spinz[5][1] = 0;
  spinz[5][2] = 0;
  correct_coefficient[5] = 0.0;

  spin[6][0] = 2;
  spin[6][1] = 2;
  spin[6][2] = 4;
  spinz[6][0] = 2;
  spinz[6][1] = -2;
  spinz[6][2] = 0;
  correct_coefficient[6] = 1 / sqrt(6.0);
  for (int i = 0; i < 7; i++) {
    float cg = clebsch_gordan(spin[i][0], spin[i][1], spin[i][2],
                              spinz[i][0], spinz[i][1], spinz[i][2]);
    COMPARE(cg, correct_coefficient[i])
      << '\n' // Using double quotes here produces an error(?!)
      << "J1: " << spin[i][0] << " Jz1: " << spinz[i][0] << "\n"
      << "J2: " << spin[i][1] << " Jz2: " << spinz[i][1] << "\n"
      << "J3: " << spin[i][2] << " Jz3: " << spinz[i][2] << "\n"
      << "CG: " << cg
      << " Correct: " << correct_coefficient[i];
  }
}


const float tolerance = 1.0e-7;

TEST (iso_clebsch_2to1) {
  const ParticleType &pip = ParticleType::find(0x211);
  const ParticleType &piz = ParticleType::find(0x111);
  const ParticleType &pim = ParticleType::find(-0x211);
  const ParticleType &rho_p = ParticleType::find(0x213);
  const ParticleType &rho_z = ParticleType::find(0x113);
  const ParticleType &rho_m = ParticleType::find(-0x213);
  const ParticleType &sigma = ParticleType::find(0x9000221);
  const ParticleType &proton  = ParticleType::find(0x2212);
  const ParticleType &neutron = ParticleType::find(0x2112);
  const ParticleType &Delta_pp = ParticleType::find(0x2224);
  const ParticleType &Delta_p  = ParticleType::find(0x2214);
  const ParticleType &Delta_z  = ParticleType::find(0x2114);
  const ParticleType &Delta_m  = ParticleType::find(0x1114);

  // π π -> X
  float iso_cg = isospin_clebsch_gordan_2to1(pip, piz, rho_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/2.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(piz, pip, rho_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/2.f, tolerance);

  iso_cg = isospin_clebsch_gordan_2to1(pip, pim, rho_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/2.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pim, pip, rho_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/2.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(piz, piz, rho_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 0.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pip, pim, sigma);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pim, pip, sigma);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(piz, piz, sigma);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);

  iso_cg = isospin_clebsch_gordan_2to1(piz, pim, rho_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/2.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pim, piz, rho_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/2.f, tolerance);


  // π N -> X
  iso_cg = isospin_clebsch_gordan_2to1(pip, proton, Delta_pp);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f, tolerance);

  iso_cg = isospin_clebsch_gordan_2to1(pip, neutron, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pip, neutron, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 2.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(piz, proton, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 2.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(piz, proton, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);

  iso_cg = isospin_clebsch_gordan_2to1(piz, neutron, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 2.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(piz, neutron, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pim, proton, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f/3.f, tolerance);
  iso_cg = isospin_clebsch_gordan_2to1(pim, proton, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 2.f/3.f, tolerance);

  iso_cg = isospin_clebsch_gordan_2to1(pim, neutron, Delta_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg*iso_cg, 1.f, tolerance);
}


TEST (iso_clebsch_2to2) {
  const ParticleType &proton  = ParticleType::find(0x2212);
  const ParticleType &neutron = ParticleType::find(0x2112);
  const ParticleType &Delta_pp = ParticleType::find(0x2224);
  const ParticleType &Delta_p  = ParticleType::find(0x2214);
  const ParticleType &Delta_z  = ParticleType::find(0x2114);
  const ParticleType &Delta_m  = ParticleType::find(0x1114);

  // N N -> N N
  float iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, proton, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, proton, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.5f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, neutron, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.5f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(neutron, neutron, neutron, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1.f, tolerance);

  // N N -> N Delta
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_pp, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3.f/4.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_p, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1.f/4.f, tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_p, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1.f/4.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_z, proton);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 1.f/4.f, tolerance);

  // N N -> Delta Delta
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_pp, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3.f/10.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_z, Delta_pp);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3.f/10.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_p, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 4.f/10.f, tolerance);

  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_pp, Delta_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 7.f/20.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_m, Delta_pp);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 7.f/20.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_p, Delta_z);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3.f/20.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, neutron, Delta_z, Delta_p);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 3.f/20.f, tolerance);

  // invalid cases
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, proton, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_z, neutron);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.f, tolerance);
  iso_cg = isospin_clebsch_gordan_sqr_2to2(proton, proton, Delta_m, Delta_m);
  COMPARE_ABSOLUTE_ERROR(iso_cg, 0.f, tolerance);
}
