/*
 *
 *    Copyright (c) 2016-2018,2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include "setup.h"

#include "../include/smash/clebschgordan.h"
#include "../include/smash/parametrizations.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

constexpr double tolerance = 1.0e-7;
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
