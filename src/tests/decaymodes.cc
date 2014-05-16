/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"

#include "../include/decaymodes.h"

using namespace Smash;

TEST_CATCH(add_no_particles, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode({}, 1.f);
}

TEST_CATCH(add_one_particle, DecayModes::InvalidDecay) {
  DecayModes m;
  m.add_mode({0}, 1.f);
}

TEST(add_two_particles) {
  DecayModes m;
  VERIFY(m.is_empty());
  m.add_mode({0, 1}, 1.f);
  VERIFY(!m.is_empty());
}
