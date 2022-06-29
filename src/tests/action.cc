/*
 *
 *    Copyright (c) 2015,2017-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <vir/test.h>  // This include has to be first

#include "../include/smash/scatteraction.h"
#include "../include/smash/scatteractionsfinder.h"
#include "setup.h"

using namespace smash;
using smash::Test::Momentum;
using smash::Test::Position;

TEST(init_particle_types) { Test::create_stable_smashon_particletypes(); }

// test collision_time:
// parallel momenta => impossible collision
TEST(impossible_collision) {
  const auto a =
      Test::smashon(Position{1., 1., 1., 1.}, Momentum{0.1, 0.3, -0.1, 0.2}, 1);
  const auto b =
      Test::smashon(Position{2., 2., 2., 2.}, Momentum{0.1, 0.3, -0.1, 0.2}, 2);

  Configuration config = Test::configuration("");
  ExperimentParameters exp_par = Test::default_parameters();
  ScatterActionsFinder finder(config, exp_par);

  VERIFY(finder.collision_time(a, b, 0.0, {}) < 0.0);
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

TEST(phasespace_five_body) {
  // Sample 5-body phase space repeatly and check if the average angles of one
  // of the outgoing particles matches an isotropic distribution
  const auto a =
      Test::smashon(Position{1., 1., 1., 1.}, Momentum{0.5, 0.0, 0., 0.});
  const auto b =
      Test::smashon(Position{2., 2., 2., 2.}, Momentum{0.5, 0.0, 0., 0.});
  const ParticleType &type_so = ParticleType::find(0x661);

  ScatterAction act(a, b, 0.);

  // Add CollisionBranch with 5 outgoing particles
  CollisionBranchPtr coll_b =
      std::make_unique<CollisionBranch>(type_so, type_so, type_so, type_so,
                                        type_so, 10.0, ProcessType::TwoToFive);
  act.add_collision(std::move(coll_b));

  const int N_samples = 10000;
  double sum_theta1 = 0;
  double sum_phi1 = 0;
  double theta1, phi1;
  for (int i = 0; i < N_samples; i++) {
    // Sample the 5-body phase space of the outgoing particles
    act.generate_final_state();

    theta1 = act.outgoing_particles()[4].momentum().threevec().get_theta();
    sum_theta1 += theta1;

    phi1 = act.outgoing_particles()[4].momentum().threevec().get_phi();
    sum_phi1 += phi1;

    // Comment out, if you want a print-out e.g. for a angle histogram
    // std::cout << "phi1 = " <<
    // act.outgoing_particles()[0].momentum().threevec().get_phi() << '\n';
    // std::cout << "theta1 = " << theta1 << '\n';
  }

  // Check if outgoing particle's phi and theta average match an uniform
  // distribution
  COMPARE_ABSOLUTE_ERROR(sum_theta1 / N_samples, M_PI / 2.0, 0.1);
  COMPARE_ABSOLUTE_ERROR(sum_phi1 / N_samples, 0.0, 0.1);
}

TEST(phasespace_manybody) {
  double sqrts = 3.0;
  std::vector<double> m{0.4, 0.3, 0.55, 0.12, 0.17};
  std::vector<FourVector> momenta(m.size());
  for (size_t i = 0; i < 1000; i++) {
    Action::sample_manybody_phasespace_impl(sqrts, m, momenta);
    for (size_t j = 0; j < m.size(); j++) {
      // Uncomment for printout
      // for (size_t k = 0; k < 4; k++) std::cout << momenta[j][k] << " ";
    }
    std::cout << std::endl;
  }
}
