/*
 *
 *    Copyright (c) 2014-2015,2017-2020,2022-2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include <filesystem>

#include "setup.h"
#include "smash/collidermodus.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

static Configuration get_common_configuration() {
  return Configuration{R"(
    General:
      Modus: ToBeSet
      End_Time: 20.1
      Nevents: 1
      Randomseed: 1
    )"};
}

static Configuration get_collider_configuration() {
  auto config = get_common_configuration();
  config.set_value(InputKeys::gen_modus, "Collider");
  /* The 'Collisions_Within_Nucleus' key is here used with its default value to
   * make sure Experiment parses it correctly (it is in the "Collider" section,
   * but it is taken by ScatterActionsFinderParameters).*/
  config.merge_yaml(R"(
    Modi:
      Collider:
        Projectile:
          Particles: {2212: 79, 2112: 118}
        Target:
          Particles: {2212: 79, 2112: 118}
        E_Kin: 1.23
      Collisions_Within_Nucleus: false
  )");
  return config;
}

TEST(create_box) {
  auto config = get_common_configuration();
  config.set_value(InputKeys::gen_modus, "Box");
  config.merge_yaml(R"(
    Modi:
      Box:
        Initial_Condition: "peaked momenta"
        Length: 10.0
        Temperature: 0.2
        Start_Time: 0.0
        Init_Multiplicities:
          661: 724
  )");
  VERIFY(!!Test::experiment(std::move(config)));
}

TEST(create_collider) {
  VERIFY(!!Test::experiment(get_collider_configuration()));
}

TEST(create_sphere) {
  auto config = get_common_configuration();
  config.set_value(InputKeys::gen_modus, "Sphere");
  config.merge_yaml(R"(
    Modi:
      Sphere:
        Initial_Condition: "thermal momenta"
        Radius: 5.0
        Temperature: 0.2
        Start_Time: 0.0
        Init_Multiplicities:
          661: 724
  )");
  VERIFY(!!Test::experiment(std::move(config)));
}

TEST(create_experiment_with_default_lattice) {
  auto config = get_collider_configuration();
  config.merge_yaml(R"(
    Lattice:
      Automatic: True
  )");
  VERIFY(!!Test::experiment(std::move(config)));
}

TEST_CATCH(create_experiment_with_invalid_modus,
           ExperimentBase::InvalidModusRequest) {
  Test::experiment(Configuration{"General: {Modus: Invalid}"});
}

TEST_CATCH(create_experiment_with_invalid_delta_time, std::invalid_argument) {
  auto config = get_collider_configuration();
  config.set_value(InputKeys::gen_deltaTime, 0.0);
  Test::experiment(std::move(config));
}

TEST_CATCH(create_experiment_with_invalid_end_time, std::invalid_argument) {
  auto config = get_collider_configuration();
  config.set_value(InputKeys::gen_endTime, 0.0);
  Test::experiment(std::move(config));
}

TEST_CATCH(create_experiment_with_invalid_output_interval,
           std::invalid_argument) {
  auto config = get_collider_configuration();
  config.merge_yaml(R"(
    Output:
      Output_Interval: 0.0
      Particles:
        Format: ["Oscar2013"]
  )");
  Test::experiment(std::move(config));
}

TEST_CATCH(run_experiment_beyond_end_time, std::logic_error) {
  auto config = get_collider_configuration();
  auto exp = std::make_unique<Experiment<ColliderModus>>(config, ".");
  exp->run_time_evolution(1000);
}

TEST(access_particles) {
  auto config = get_collider_configuration();
  auto exp = std::make_unique<Experiment<ColliderModus>>(config, ".");
  Particles* part = exp->first_ensemble();
  part->create(0x211);
  ParticleList part_list = part->copy_to_vector();
  VERIFY(part_list.size() == 1);
}

TEST(add_and_remove_particles) {
  /*
   * NOTE: As in the first test of this file Test::create_actual_particletypes
   *       is called, only particles with valid PDG codes can be easily used.
   *       Although ParticleData{ParticleType{"Inv", 0, 0, Parity::Neg, 0x0}}
   *       seems valid it would fail in Debug mode because of an assert in the
   *       overload of the address operator of ParticleType. Therefore, adding
   *       or removing invalid particles is not tested here.
   */

  // Set up collider experiment without setting up initial state (no Au-Au)
  auto config = get_collider_configuration();
  auto exp = std::make_unique<Experiment<ColliderModus>>(config, ".");

  // Neither add nor remove particles -> expect 0 particles
  exp->run_time_evolution(1., ParticleList{}, ParticleList{});
  VERIFY(exp->first_ensemble()->size() == 0);

  // Add 1 pion off shell -> expect pion added on shell
  ParticleData pion_plus{ParticleType::find(pdg::pi_p)};
  pion_plus.set_4momentum(FourVector(1.0, 0.95, 0.0, 0.0));
  pion_plus.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  exp->run_time_evolution(1., ParticleList{pion_plus}, ParticleList{});
  VERIFY(exp->first_ensemble()->size() == 1);
  COMPARE_ABSOLUTE_ERROR(
      exp->first_ensemble()->begin()->momentum().x0(),
      std::sqrt(0.95 * 0.95 + pion_plus.pole_mass() * pion_plus.pole_mass()),
      very_small_double);

  // Remove existing pion -> expect 0 particles
  exp->run_time_evolution(1., ParticleList{},
                          exp->first_ensemble()->copy_to_vector());
  VERIFY(exp->first_ensemble()->size() == 0);

  // Add pion and omega off shell -> expect 2 particles (omega on its shell)
  ParticleData omega{ParticleType::find(pdg::omega)};
  omega.set_4momentum(FourVector(0.783, 0.5, 0.0, 0.0));
  omega.set_4position(FourVector(0.0, 1.0, 0.0, 0.0));
  exp->run_time_evolution(1., ParticleList{pion_plus, omega}, ParticleList{});
  VERIFY(exp->first_ensemble()->size() == 2);
  COMPARE_ABSOLUTE_ERROR(exp->first_ensemble()->back().momentum().abs(),
                         omega.effective_mass(), very_small_double);

  // Remove non existing eta -> still 2 particles
  ParticleData eta{ParticleType::find(pdg::eta)};
  eta.set_4momentum(FourVector(0.548, 0.0, 0.0, 0.0));
  eta.set_4position(FourVector(0.0, 2.0, 0.0, 0.0));
  exp->run_time_evolution(1., ParticleList{}, ParticleList{eta});
  VERIFY(exp->first_ensemble()->size() == 2);

  /* Removing a particle, which is not a direct copy of the
   * particle in ensembles_[0], but has only PDG code, position and momentum
   * set.*/
  exp->run_time_evolution(1., ParticleList{}, ParticleList{pion_plus});
  VERIFY(exp->first_ensemble()->size() == 1);
}

TEST_CATCH(remove_particle_twice, std::logic_error) {
  // Set up collider experiment without setting up initial state (no Au-Au)
  auto config = get_collider_configuration();
  auto exp = std::make_unique<Experiment<ColliderModus>>(config, ".");

  // Add an eta
  ParticleData eta{ParticleType::find(pdg::eta)};
  eta.set_4momentum(FourVector(0.548, 0.0, 0.0, 0.0));
  eta.set_4position(FourVector(0.0, 2.0, 0.0, 0.0));
  exp->run_time_evolution(1., ParticleList{eta}, ParticleList{});

  // Try to remove the eta twice
  exp->run_time_evolution(1., ParticleList{}, ParticleList{eta, eta});
}
