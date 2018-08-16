/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

#include "../include/smash/boxmodus.h"
#include "../include/smash/collidermodus.h"
#include "../include/smash/modusdefault.h"
#include "../include/smash/potentials.h"
#include "../include/smash/propagation.h"
#include "../include/smash/spheremodus.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

TEST(init_particle_types) { Test::create_smashon_particletypes(); }

// create a particle list with various interesting particles. We will
// assume a box of 5 fm length and a time step (for propagation) of 1 fm.
static Test::ParticlesPtr create_box_particles() {
  return Test::create_particles(
      {// particle that doesn't move:
       Test::smashon(Position{0.0, 0.6, 0.7, 0.8},
                     Momentum{4.0, 0.0, 0.0, 0.0}),
       // particle that moves with speed of light
       Test::smashon(Position{0.0, 0.7, 0.8, 0.9},
                     Momentum{sqrt(0.03), 0.1, -.1, 0.0}),
       // particle that moves slowly:
       Test::smashon(Position{0.0, 0.1, 0.2, 0.3},
                     Momentum{sqrt(1.14), 0.1, 0.2, -.3}),
       // particle that will cross a box boundary at high x:
       Test::smashon(Position{0.0, 4.5, 0.0, 0.0},
                     Momentum{0.11, 0.1, 0.0, 0.0}),
       // particle that will cross a box boundary at low y:
       Test::smashon(Position{0.0, 0.0, 0.2, 0.0},
                     Momentum{0.11, 0.0, -.1, 0.0}),
       // particle that will cross a box boundary at low x and high z:
       Test::smashon(Position{0.0, 0.2, 0.0, 4.8},
                     Momentum{0.51, -.3, 0.0, 0.4})});
}

TEST(propagate_default_no_potentials) {
  auto Pdef = create_box_particles();
  propagate_straight_line(Pdef.get(), 1.0, {});
  // after propagation: Momenta should be unchanged.
  // but positions should be updated.
  auto it = Pdef->begin();
  COMPARE(it->momentum(), FourVector(4.0, 0.0, 0.0, 0.0));
  COMPARE(it->position(), FourVector(1.0, 0.6, 0.7, 0.8));
  ++it;
  COMPARE(it->momentum(), FourVector(sqrt(0.03), 0.1, -.1, 0.0));
  COMPARE(it->position(), FourVector(1.0, 0.7 + 0.1 / std::sqrt(0.03),
                                     0.8 - 0.1 / std::sqrt(0.03), 0.9));
  ++it;
  COMPARE(it->momentum(), FourVector(sqrt(1.14), 0.1, 0.2, -.3));
  COMPARE(it->position(),
          FourVector(1., 0.1 + 0.1 / std::sqrt(1.14),
                     0.2 + 0.2 / std::sqrt(1.14), 0.3 - 0.3 / std::sqrt(1.14)));
  ++it;
  COMPARE(it->momentum(), FourVector(0.11, 0.1, 0.0, 0.0));
  COMPARE(it->position(), FourVector(1.0, 4.5 + 0.1 / 0.11, 0.0, 0.0));
  ++it;
  COMPARE(it->momentum(), FourVector(0.11, 0.0, -.1, 0.0));
  COMPARE(it->position(), FourVector(1.0, 0.0, 0.2 - 0.1 / 0.11, 0.0));
  ++it;
  COMPARE(it->momentum(), FourVector(0.51, -.3, 0.0, 0.4));
  COMPARE(it->position(),
          FourVector(1.0, 0.2 - 0.3 / 0.51, 0.0, 4.8 + 0.4 / 0.51));
}
