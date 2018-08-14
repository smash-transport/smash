/*
 *
 *    Copyright (c) 2014-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include <boost/filesystem.hpp>

#include "../include/smash/collidermodus.h"
#include "setup.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

TEST(create_box) {
  VERIFY(!!Test::experiment(
      Configuration("General:\n"
                    "  Modus: Box\n"
                    "  End_Time: 100.0\n"
                    "  Nevents: 1\n"
                    "  Randomseed: 1\n"
                    "Collision_Term:\n"
                    "  Strings: False\n"
                    "Modi: \n"
                    "  Box:\n"
                    "    Initial_Condition: \"peaked momenta\"\n"
                    "    Length: 10.0\n"
                    "    Temperature: 0.2\n"
                    "    Start_Time: 0.0\n"
                    "    Init_Multiplicities:\n"
                    "      661: 724\n")));
}

TEST(create_collider) {
  VERIFY(!!Test::experiment(
      Configuration("General:\n"
                    "  Modus: Collider\n"
                    "  End_Time: 100.0\n"
                    "  Nevents: 1\n"
                    "  Randomseed: 1\n"
                    "Collision_Term:\n"
                    "  Strings: False\n"
                    "Modi: \n"
                    "  Collider:\n"
                    "    Projectile: \n"
                    "      Particles: {2212: 79, 2112: 118}\n"
                    "    Target: \n"
                    "      Particles: {2212: 79, 2112: 118}\n"
                    "    E_Kin: 1.23\n")));
}

TEST(create_sphere) {
  VERIFY(!!Test::experiment(
      Configuration("General:\n"
                    "  Modus: Sphere\n"
                    "  End_Time: 100.0\n"
                    "  Nevents: 1\n"
                    "  Randomseed: 1\n"
                    "Collision_Term:\n"
                    "  Strings: False\n"
                    "Modi: \n"
                    "  Sphere:\n"
                    "    Initial_Condition: \"thermal momenta\"\n"
                    "    Radius: 5.0\n"
                    "    Sphere_Temperature: 0.2\n"
                    "    Start_Time: 0.0\n"
                    "    Init_Multiplicities:\n"
                    "      661: 724\n")));
}

TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  Test::experiment("General: {Modus: Invalid}");
}

TEST(access_particles) {
  Configuration config = Test::configuration();
  boost::filesystem::path output_path(".");
  auto exp = make_unique<Experiment<ColliderModus>>(config, output_path);
  Particles* part = exp->particles();
  part->create(0x211);
  ParticleList part_list = part->copy_to_vector();
  VERIFY(part_list.size() == 1);
}
