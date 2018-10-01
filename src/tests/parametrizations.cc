/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/clebschgordan.h"
#include "../include/smash/parametrizations.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

constexpr double tolerance = 1.0e-7;

TEST(clebsch_strangeness_exchange) {
  const auto& proton = ParticleType::find(0x2212);
  const auto& neutron = ParticleType::find(0x2112);
  const auto& K_m = ParticleType::find(-0x321);
  const auto& Lambda = ParticleType::find(0x3122);
  const auto& pi_p = ParticleType::find(0x211);
  const auto& pi_z = ParticleType::find(0x111);
  const auto& pi_m = ParticleType::find(-0x211);
  const auto& Sigma_z = ParticleType::find(0x3212);
  const auto& Sigma_m = ParticleType::find(0x3112);

  const double sqrts = 2.0;

  {
    const auto ratio = kminusn_piminussigma0(sqrts) / kminusp_pi0sigma0(sqrts);
    const auto expected_ratio =
        isospin_clebsch_gordan_sqr_2to2(K_m, proton, pi_z, Sigma_z) /
        isospin_clebsch_gordan_sqr_2to2(K_m, neutron, pi_m, Sigma_z);
    COMPARE_ABSOLUTE_ERROR(ratio, expected_ratio, tolerance);
  }
  {
    const auto ratio =
        kminusn_pi0sigmaminus(sqrts) / kminusp_piplussigmaminus(sqrts);
    const auto expected_ratio =
        isospin_clebsch_gordan_sqr_2to2(K_m, proton, pi_p, Sigma_m) /
        isospin_clebsch_gordan_sqr_2to2(K_m, neutron, pi_z, Sigma_m);
    COMPARE_ABSOLUTE_ERROR(ratio, expected_ratio, tolerance);
  }
  {
    const auto ratio = kminusn_piminuslambda(sqrts) / kminusp_pi0lambda(sqrts);
    const auto expected_ratio =
        isospin_clebsch_gordan_sqr_2to2(K_m, proton, pi_z, Lambda) /
        isospin_clebsch_gordan_sqr_2to2(K_m, neutron, pi_m, Lambda);
    COMPARE_ABSOLUTE_ERROR(ratio, expected_ratio, tolerance);
  }
}

TEST(clebsch_kaon_charge_exchange) {
  const auto& proton = ParticleType::find(0x2212);
  const auto& neutron = ParticleType::find(0x2112);
  const auto& K_p = ParticleType::find(0x321);
  const auto& K_z = ParticleType::find(0x311);

  const auto cg1 = isospin_clebsch_gordan_sqr_2to2(neutron, K_p, K_z, proton);
  const auto cg2 = isospin_clebsch_gordan_sqr_2to2(proton, K_z, K_p, neutron);
  // We assume they are same in crosssections.cc.
  COMPARE_ABSOLUTE_ERROR(cg1, cg2, tolerance);
}
