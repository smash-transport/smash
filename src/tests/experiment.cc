/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

TEST(create_box) { VERIFY(!!Test::experiment("General: {Modus: Box}")); }

TEST(create_collider) {
  VERIFY(!!Test::experiment("General: {Modus: Collider}"));
}

TEST(create_sphere) { VERIFY(!!Test::experiment("General: {Modus: Sphere}")); }

TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  Test::experiment("General: {Modus: Invalid}");
}
