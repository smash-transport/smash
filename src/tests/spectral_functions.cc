/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"
#include "setup.h"

#include "../include/formfactors.h"
#include "../include/integrate.h"
#include "../include/stringfunctions.h"

using namespace Smash;
using namespace UnitTest;

TEST(spectral_functions) {
  Smash::Test::create_actual_particletypes();
  Smash::Test::create_actual_decaymodes();

  Integrator integrate;
  /* Upper limit for integration in GeV. Hadron masses are usually a few GeV,
   * so this should really be on the safe side. */
  const float max_mass = 100.;
  // error tolerance (max. deviation from one)
  const float error_tolerance = 0.25;
  const float warning_level = 0.21;
  // error tolerance (max. deviation from one) for constant-width SF
  const float error_tolerance_const = 0.003;

  /* Loop over all resonances. */
  for (const ParticleType &type : ParticleType::list_all()) {
    if (type.is_stable()) {
      continue;
    }
    /* Integrate spectral function. */
    const auto result = integrate(type.minimum_mass(), max_mass,
                          [&](double m) {
                            return type.spectral_function(m);
                          });
    const auto result_const = integrate(0., max_mass,
                          [&](double m) {
                            return type.spectral_function_const_width(m);
                          });
    if (result.value() > 1 + warning_level) {
      std::cout << AnsiColor::blue;
    } else if (result.value() < 1 - warning_level) {
      std::cout << AnsiColor::yellow;
    }
    std::cout << fill_right(type.name(), 11) << ": "
              << format(result.value(), nullptr, -1, 4) << " ± " << result.error() << ", "
              << result_const.value() << " ± " << result_const.error()
              << AnsiColor::normal << "\n";
    // check if integral is approximately equal to one
    COMPARE_ABSOLUTE_ERROR(result.value(), 1., error_tolerance);
    COMPARE_ABSOLUTE_ERROR(result_const.value(), 1., error_tolerance_const);
  }
}
