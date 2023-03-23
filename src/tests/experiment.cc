/*
 *
 *    Copyright (c) 2014-2015,2017-2020,2022
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
      End_Time: 100.0
      Nevents: 1
      Randomseed: 1
    )"};
}

static Configuration get_collider_configuration() {
  auto config = get_common_configuration();
  config.set_value({"General", "Modus"}, "Collider");
  config.merge_yaml(R"(
    Modi:
      Collider:
        Projectile:
          Particles: {2212: 79, 2112: 118}
        Target:
          Particles: {2212: 79, 2112: 118}
        E_Kin: 1.23
  )");
  return config;
}

static ParticleData create_pion() {
  return ParticleData{ParticleType::find(0x211)};
}

static ParticleData create_eta() {
  return ParticleData{ParticleType::find(0x221)};
}

static ParticleData create_omega() {
  return ParticleData{ParticleType::find(0x223)};
}

TEST(create_box) {
  auto config = get_common_configuration();
  config.set_value({"General", "Modus"}, "Box");
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
  config.set_value({"General", "Modus"}, "Sphere");
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

TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  Test::experiment(Configuration{"General: {Modus: Invalid}"});
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
  auto config = get_collider_configuration();
  auto exp = std::make_unique<Experiment<ColliderModus>>(config, ".");
  exp->run_time_evolution(1., ParticleList{}, ParticleList{});
  VERIFY(exp->first_ensemble()->size() == 0);

  ParticleData part_1 = create_pion();
  part_1.set_4momentum(FourVector(1.0, 0.95, 0.0, 0.0));
  part_1.set_4position(FourVector(0.0, 0.0, 0.0, 0.0));
  ParticleList P1{part_1};
  exp->run_time_evolution(1., P1, ParticleList{});
  VERIFY(exp->first_ensemble()->size() == 1);

  Particles* part = exp->first_ensemble();
  ParticleList part_list = part->copy_to_vector();
  exp->run_time_evolution(1., ParticleList{}, part_list);
  VERIFY(exp->first_ensemble()->size() == 0);

  ParticleData part_2 = create_omega();
  part_2.set_4momentum(FourVector(0.783, 0.5, 0.0, 0.0));
  part_2.set_4position(FourVector(0.0, 1.0, 0.0, 0.0));
  ParticleList P2{part_1, part_2};
  exp->run_time_evolution(1., P2, ParticleList{});
  VERIFY(exp->first_ensemble()->size() == 2);

  ParticleData part_3 = create_eta();
  part_3.set_4momentum(FourVector(0.548, 0.0, 0.0, 0.0));
  part_3.set_4position(FourVector(0.0, 2.0, 0.0, 0.0));
  ParticleList P3{part_3};
  exp->run_time_evolution(1., ParticleList{}, P3);
  VERIFY(exp->first_ensemble()->size() == 2);
}
