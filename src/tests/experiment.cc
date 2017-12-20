/*
 *
 *    Copyright (c) 2014-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "setup.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

TEST(create_box) { VERIFY(!!Test::experiment(
                        Configuration("General:\n"
                         "  Modus: Box\n"
                         "  End_Time: 100.0\n"
                         "  Nevents: 1\n"
                         "Modi: \n"
                         "  Box:\n"
                         "    Initial_Condition: \"peaked momenta\"\n"
                         "    Length: 10.0\n"
                         "    Temperature: 0.2\n"
                         "    Start_Time: 0.0\n"
                         "    Init_Multiplicities:\n"
                         "      661: 724\n"))); }




TEST(create_collider) {
  VERIFY(!!Test::experiment("General: {Modus: Collider}"));
}

TEST(create_sphere) { VERIFY(!!Test::experiment(
                        Configuration("General:\n"
                         "  Modus: Sphere\n"
                         "  End_Time: 100.0\n"
                         "  Nevents: 1\n"
                         "Modi: \n"
                         "  Sphere:\n"
                         "    Initial_Condition: \"thermal momenta\"\n"
                         "    Radius: 5.0\n"
                         "    Sphere_Temperature: 0.2\n"
                         "    Start_Time: 0.0\n"
                         "    Init_Multiplicities:\n"
                         "      661: 724\n"))); }

TEST_CATCH(create_invalid, ExperimentBase::InvalidModusRequest) {
  Test::experiment("General: {Modus: Invalid}");
}
