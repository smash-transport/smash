/*
 *
 *    Copyright (c) 2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "vir/test.h"  // This include has to be first

#include "smash/spheremodus.h"

#include "setup.h"
#include "smash/experimentparameters.h"

using namespace smash;

TEST(init_particle_types) { Test::create_actual_particletypes(); }

static Configuration get_sphere_configuration() {
  return Configuration{R"(
    Modus:
      Sphere:
        Radius: 10
        Start_Time: 0
        Temperature: 50
        Init_Multiplicities:
          111: 1
    )"};
}

TEST_CATCH(no_negative_temperature, std::invalid_argument) {
  auto config = get_sphere_configuration();
  config.set_value(InputKeys::modi_sphere_temperature, -42);
  SphereModus sphere(std::move(config), ExperimentParameters{});
}

TEST_CATCH(radial_flow_causal, std::invalid_argument) {
  auto config = get_sphere_configuration();
  config.set_value(InputKeys::modi_sphere_addRadialVelocity, 42);
  SphereModus sphere(std::move(config), ExperimentParameters{});
}

TEST_CATCH(radial_flow_exponent_positive, std::invalid_argument) {
  auto config = get_sphere_configuration();
  config.set_value(InputKeys::modi_sphere_addRadialVelocityExponent, -42);
  SphereModus sphere(std::move(config), ExperimentParameters{});
}
