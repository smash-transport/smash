/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/scatteraction.h"
#include "../include/smash/scatteractionsfinder.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

TEST(init_particle_types) { Test::create_smashon_particletypes(); }

// test collision_time:
// parallel momenta => impossible collision
TEST(impossible_collision) {
  const auto a =
      Test::smashon(Position{1., 1., 1., 1.}, Momentum{0.1, 0.3, -0.1, 0.2});
  const auto b =
      Test::smashon(Position{2., 2., 2., 2.}, Momentum{0.1, 0.3, -0.1, 0.2});

  VERIFY(ScatterActionsFinder::collision_time(a, b) < 0.0);
}

// test particle_distance:
// particles with null momenta
TEST(particle_distance) {
  const auto a =
      Test::smashon(Position{1., 1., 1., 1.}, Momentum{0.1, 0., 0., 0.});
  const auto b =
      Test::smashon(Position{2., 2., 2., 2.}, Momentum{0.1, 0., 0., 0.});

  ScatterAction act(a, b, 0.);
  const auto distance_squared = act.transverse_distance_sqr();
  VERIFY(distance_squared >= 0.);
  VERIFY(distance_squared <= 100.);
}

// particles with finite momenta
TEST(finite_momenta) {
  const auto a =
      Test::smashon(Position{1., 1., 1., 1.}, Momentum{0.1, 10.0, 9.0, 8.0});
  const auto b = Test::smashon(Position{2., 2., 2., 2.},
                               Momentum{0.1, -10.0, -90.0, -80.0});
  ScatterAction act(a, b, 0.);
  VERIFY(act.transverse_distance_sqr() >= 0.);
}
