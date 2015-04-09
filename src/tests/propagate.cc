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
#include "../include/spheremodus.h"

using namespace Smash;
using Smash::Test::Position;
using Smash::Test::Momentum;

TEST(init_particle_types) {
  Test::create_smashon_particletypes();
}

// create a particle list with various interesting particles. We will
// assume a box of 5 fm length and a time step (for propagation) of 1
// fm.
static Test::ParticlesPtr create_box_particles() {
  return Test::create_particles(
      {// particle that doesn't move:
       Test::smashon(Position{0.0, 0.6, 0.7, 0.8},
                     Momentum{4.0, 0.0, 0.0, 0.0}),
       // particle that moves with speed of light
       Test::smashon(Position{0.5, 0.7, 0.8, 0.9},
                     Momentum{sqrt(0.02), 0.1, -.1, 0.0}),
       // particle that moves slowly:
       Test::smashon(Position{0.7, 0.1, 0.2, 0.3},
                     Momentum{sqrt(1.13), 0.1, 0.2, -.3}),
       // particle that will cross a box boundary at high x:
       Test::smashon(Position{1.2, 4.5, 0.0, 0.0},
                     Momentum{0.1, 0.1, 0.0, 0.0}),
       // particle that will cross a box boundary at low y:
       Test::smashon(Position{1.8, 0.0, 0.2, 0.0},
                     Momentum{0.1, 0.0, -.1, 0.0}),
       // particle that will cross a box boundary at low x and high z:
       Test::smashon(Position{2.2, 0.2, 0.0, 4.8},
                     Momentum{0.5, -.3, 0.0, 0.4})});
}

TEST(propagate_default) {
  ModusDefault m;
  auto Pdef = create_box_particles();
  OutputsList out;
  // clock, output interval, cross-section, testparticles
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  Potentials* pot = NULL;
  m.propagate(Pdef.get(), param, out, pot);
  // after propagation: Momenta should be unchanged.
  COMPARE(Pdef->data(0).momentum(), FourVector(4.0, 0.0, 0.0, 0.0));
  COMPARE(Pdef->data(1).momentum(), FourVector(sqrt(0.02), 0.1, -.1, 0.0));
  COMPARE(Pdef->data(2).momentum(), FourVector(sqrt(1.13), 0.1, 0.2, -.3));
  COMPARE(Pdef->data(3).momentum(), FourVector(0.1, 0.1, 0.0, 0.0));
  COMPARE(Pdef->data(4).momentum(), FourVector(0.1, 0.0, -.1, 0.0));
  COMPARE(Pdef->data(5).momentum(), FourVector(0.5, -.3, 0.0, 0.4));
  // positions should be updated:
  COMPARE(Pdef->data(0).position(), FourVector(0.0, 0.6, 0.7, 0.8));
  COMPARE(Pdef->data(1).position(),
          FourVector(0.0, 0.7 + std::sqrt(0.5), 0.8 - std::sqrt(0.5), 0.9));
  COMPARE(Pdef->data(2).position(),
          FourVector(0., 0.1 + 0.1 / std::sqrt(1.13),
                     0.2 + 0.2 / std::sqrt(1.13), 0.3 - 0.3 / std::sqrt(1.13)));
  COMPARE(Pdef->data(3).position(), FourVector(0.0, 4.5 + 1.0, 0.0, 0.0));
  COMPARE(Pdef->data(4).position(), FourVector(0.0, 0.0, 0.2 - 1.0, 0.0));
  COMPARE(Pdef->data(5).position(), FourVector(0.0, 0.2 - 0.6, 0.0, 4.8 + 0.8));
}

TEST(propagate_box) {
  ModusDefault m;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  BoxModus b(
      "Box: { Initial_Condition: 1, Length: 5.0, Temperature: 0.13, "
      "Start_Time: 0.2, Init_Multiplicities: {661: 10} }",
      param);
  auto Pdef = create_box_particles();
  auto Pbox = create_box_particles();
  OutputsList out;
  // clock, output interval, cross-section, testparticles
  Potentials* pot = NULL;
  m.propagate(Pdef.get(), param, out, pot);
  b.propagate(Pbox.get(), param, out, pot);
  // Box and Default modus should do the same to momentum:
  COMPARE(Pdef->data(0).momentum(), Pbox->data(0).momentum());
  COMPARE(Pdef->data(1).momentum(), Pbox->data(1).momentum());
  COMPARE(Pdef->data(2).momentum(), Pbox->data(2).momentum());
  COMPARE(Pdef->data(3).momentum(), Pbox->data(3).momentum());
  COMPARE(Pdef->data(4).momentum(), Pbox->data(4).momentum());
  COMPARE(Pdef->data(5).momentum(), Pbox->data(5).momentum());
  // stop, fast and slow are propagated equally:
  COMPARE(Pdef->data(0).position(), Pbox->data(0).position());
  COMPARE(Pdef->data(1).position(), Pbox->data(1).position());
  COMPARE(Pdef->data(2).position(), Pbox->data(2).position());
  // those particles are expected to be wraped around in Box Modus
  COMPARE(Pdef->data(3).position(),
          Pbox->data(3).position() + FourVector(0.0, 5.0, 0.0, 0.0));
  COMPARE(Pdef->data(4).position(),
          Pbox->data(4).position() + FourVector(0.0, 0.0, -5., 0.0));
  COMPARE(Pdef->data(5).position(),
          Pbox->data(5).position() + FourVector(0.0, -5., 0.0, 5.0));
}

TEST(propagate_collider) {
  ModusDefault m;
  ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
  ColliderModus c(
      "Collider: { Projectile: { Particles: {661: 1} }, "
      "Target: { Particles: {661: 1} }, Sqrtsnn: 1.0 }",
      param);
  auto Pdef = create_box_particles();
  auto Pcol = create_box_particles();
  OutputsList out;
  Potentials* pot = NULL;
  m.propagate(Pdef.get(), param, out, pot);
  c.propagate(Pcol.get(), param, out, pot);
  // Collider and Default modus should do the same everywhere:
  COMPARE(Pdef->data(0).momentum(), Pcol->data(0).momentum());
  COMPARE(Pdef->data(1).momentum(), Pcol->data(1).momentum());
  COMPARE(Pdef->data(2).momentum(), Pcol->data(2).momentum());
  COMPARE(Pdef->data(3).momentum(), Pcol->data(3).momentum());
  COMPARE(Pdef->data(4).momentum(), Pcol->data(4).momentum());
  COMPARE(Pdef->data(5).momentum(), Pcol->data(5).momentum());
  COMPARE(Pdef->data(0).position(), Pcol->data(0).position());
  COMPARE(Pdef->data(1).position(), Pcol->data(1).position());
  COMPARE(Pdef->data(2).position(), Pcol->data(2).position());
  COMPARE(Pdef->data(3).position(), Pcol->data(3).position());
  COMPARE(Pdef->data(4).position(), Pcol->data(4).position());
  COMPARE(Pdef->data(5).position(), Pcol->data(5).position());
}

TEST(propagate_sphere) {
   ModusDefault m;
   ExperimentParameters param{{0.f, 1.f}, 1.f, 1, 1.0};
   SphereModus s(
       "Sphere: { Radius: 10, Start_Time: 0.0, Sphere_Temperature: 0.2, "
       "Init_Multiplicities: {661: 500} }",
       param);
   auto Pdef = create_box_particles();
   auto Psph = create_box_particles();
   OutputsList out;
   Potentials* pot = NULL;
   m.propagate(Pdef.get(), param, out, pot);
   s.propagate(Psph.get(), param, out, pot);
   // Sphere and Default modus should do the same everywhere:
   COMPARE(Pdef->data(0).momentum(), Psph->data(0).momentum());
   COMPARE(Pdef->data(1).momentum(), Psph->data(1).momentum());
   COMPARE(Pdef->data(2).momentum(), Psph->data(2).momentum());
   COMPARE(Pdef->data(3).momentum(), Psph->data(3).momentum());
   COMPARE(Pdef->data(4).momentum(), Psph->data(4).momentum());
   COMPARE(Pdef->data(5).momentum(), Psph->data(5).momentum());
   COMPARE(Pdef->data(0).position(), Psph->data(0).position());
   COMPARE(Pdef->data(1).position(), Psph->data(1).position());
   COMPARE(Pdef->data(2).position(), Psph->data(2).position());
   COMPARE(Pdef->data(3).position(), Psph->data(3).position());
   COMPARE(Pdef->data(4).position(), Psph->data(4).position());
   COMPARE(Pdef->data(5).position(), Psph->data(5).position());
}
