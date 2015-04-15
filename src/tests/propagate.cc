/*
 *
 *    Copyright (c) 2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include "../include/boxmodus.h"
#include "../include/collidermodus.h"
#include "../include/modusdefault.h"
#include "../include/potentials.h"
#include "../include/propagation.h"
#include "../include/spheremodus.h"

using namespace Smash;
using Smash::Test::Position;
using Smash::Test::Momentum;

TEST(init_particle_types) {
  Test::create_smashon_particletypes();
}

// create a particle list with various interesting particles. We will
// assume a box of 5 fm length and a time step (for propagation) of 1 fm.
static Test::ParticlesPtr create_box_particles() {
  return Test::create_particles(
      {// particle that doesn't move:
       Test::smashon(Position{0.0, 0.6, 0.7, 0.8},
                     Momentum{4.0, 0.0, 0.0, 0.0}),
       // particle that moves with speed of light
       Test::smashon(Position{0.5, 0.7, 0.8, 0.9},
                     Momentum{sqrt(0.03), 0.1, -.1, 0.0}),
       // particle that moves slowly:
       Test::smashon(Position{0.7, 0.1, 0.2, 0.3},
                     Momentum{sqrt(1.14), 0.1, 0.2, -.3}),
       // particle that will cross a box boundary at high x:
       Test::smashon(Position{1.2, 4.5, 0.0, 0.0},
                     Momentum{0.11, 0.1, 0.0, 0.0}),
       // particle that will cross a box boundary at low y:
       Test::smashon(Position{1.8, 0.0, 0.2, 0.0},
                     Momentum{0.11, 0.0, -.1, 0.0}),
       // particle that will cross a box boundary at low x and high z:
       Test::smashon(Position{2.2, 0.2, 0.0, 4.8},
                     Momentum{0.51, -.3, 0.0, 0.4})});
}

static Potentials create_potential() {
  Configuration conf(TEST_CONFIG_PATH);
  conf["Potentials"]["Skyrme"]["Skyrme_A"] = -209.2;
  conf["Potentials"]["Skyrme"]["Skyrme_B"] = 156.4;
  conf["Potentials"]["Skyrme"]["Skyrme_Tau"] = 1.35;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  return Potentials(conf["Potentials"], param);
}

TEST(propagate_default_no_potentials) {
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  auto Pdef = create_box_particles();
  propagate_straight_line(Pdef.get(), param);
  // after propagation: Momenta should be unchanged.
  COMPARE(Pdef->data(0).momentum(), FourVector(4.0, 0.0, 0.0, 0.0));
  COMPARE(Pdef->data(1).momentum(), FourVector(sqrt(0.03), 0.1, -.1, 0.0));
  COMPARE(Pdef->data(2).momentum(), FourVector(sqrt(1.14), 0.1, 0.2, -.3));
  COMPARE(Pdef->data(3).momentum(), FourVector(0.11, 0.1, 0.0, 0.0));
  COMPARE(Pdef->data(4).momentum(), FourVector(0.11, 0.0, -.1, 0.0));
  COMPARE(Pdef->data(5).momentum(), FourVector(0.51, -.3, 0.0, 0.4));
  // positions should be updated:
  COMPARE(Pdef->data(0).position(), FourVector(0.0, 0.6, 0.7, 0.8));
  COMPARE(Pdef->data(1).position(), FourVector(0.0, 0.7 + 0.1/std::sqrt(0.03),
                                                   0.8 - 0.1/std::sqrt(0.03), 0.9));
  COMPARE(Pdef->data(2).position(), FourVector(0., 0.1 + 0.1 / std::sqrt(1.14),
                                                  0.2 + 0.2 / std::sqrt(1.14),
                                                  0.3 - 0.3 / std::sqrt(1.14)));
  COMPARE(Pdef->data(3).position(), FourVector(0.0, 4.5 + 0.1/0.11, 0.0, 0.0));
  COMPARE(Pdef->data(4).position(), FourVector(0.0, 0.0, 0.2 - 0.1/0.11, 0.0));
  COMPARE(Pdef->data(5).position(), FourVector(0.0, 0.2 - 0.3/0.51, 0.0, 4.8 + 0.4/0.51));
}

TEST(propagate_collider) {
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  SphereModus s(
       "Sphere: { Radius: 10, Start_Time: 0.0, Sphere_Temperature: 0.2, "
       "Init_Multiplicities: {661: 500} }",
       param);
  ColliderModus c(
      "Collider: { Projectile: { Particles: {661: 1} }, "
      "Target: { Particles: {661: 1} }, Sqrtsnn: 1.0 }",
      param);

  auto Pdef = create_box_particles();
  auto Pcol = create_box_particles();
  Potentials pot = create_potential();
  propagate(Pdef.get(), param, pot, s);
  propagate(Pcol.get(), param, pot, c);
  // Collider and Default modus should do the same everywhere:
  for (size_t i = 0; i < 6; i++) {
    COMPARE(Pdef->data(i).momentum(), Pcol->data(i).momentum());
    COMPARE(Pdef->data(i).position(), Pcol->data(i).position());
  }
}

TEST(propagate_box) {
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  SphereModus s(
       "Sphere: { Radius: 10, Start_Time: 0.0, Sphere_Temperature: 0.2, "
       "Init_Multiplicities: {661: 500} }",
       param);
  BoxModus b(
      "Box: { Initial_Condition: 1, Length: 5.0, Temperature: 0.13, "
      "Start_Time: 0.2, Init_Multiplicities: {661: 10} }",
      param);
  auto Pdef = create_box_particles();
  auto Pbox = create_box_particles();
  Potentials pot = create_potential();
  propagate(Pdef.get(), param, pot, s);
  propagate(Pbox.get(), param, pot, b);
  // Now wrapping, i.e. imposing initial conditions is separated from
  // propagation, so box should produce the same results that sphere
  for (size_t i = 0; i < 6; i++) {
    COMPARE(Pdef->data(i).momentum(), Pbox->data(i).momentum());
    COMPARE(Pdef->data(i).position(), Pbox->data(i).position());
  }
}
