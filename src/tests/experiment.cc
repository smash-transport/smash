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

Configuration get_common_configuration() {
  return Configuration{R"(
    General:
      Modus: ToBeSet
      End_Time: 100.0
      Nevents: 1
      Randomseed: 1
    )"};
}

Configuration get_collider_configuration() {
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
  Test::experiment("General: {Modus: Invalid}");
}

TEST(access_particles) {
  auto config = get_collider_configuration();
  auto exp = std::make_unique<Experiment<ColliderModus>>(config, ".");
  Particles* part = exp->first_ensemble();
  part->create(0x211);
  ParticleList part_list = part->copy_to_vector();
  VERIFY(part_list.size() == 1);
}
