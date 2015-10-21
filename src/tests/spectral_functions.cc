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
  const float error_tolerance_no_norm = 0.25;
  const float warning_level = 0.21;
  // error tolerance (max. deviation from one) for constant-width SF
  const float error_tolerance_const = 0.003;

  /* Loop over all resonances. */
  for (const ParticleType &type : ParticleType::list_all()) {
    if (type.is_stable()) {
      continue;
    }
    /* Integrate spectral function. */
    const auto result_no_norm = integrate(type.minimum_mass(), max_mass,
                          [&](double m) {
                            return type.spectral_function_no_norm(m);
                          });
    const auto result_const = integrate(0., max_mass,
                          [&](double m) {
                            return type.spectral_function_const_width(m);
                          });
    const auto result = integrate(type.minimum_mass(), max_mass,
                          [&](double m) {
                            return type.spectral_function(m);
                          });
    if (result_no_norm.value() > 1 + warning_level) {
      std::cout << AnsiColor::blue;
    } else if (result_no_norm.value() < 1 - warning_level) {
      std::cout << AnsiColor::yellow;
    }
    std::cout << fill_right(type.name(), 11) << ": "
              << format(result_no_norm.value(), nullptr, -1, 4) << " ± " << result_no_norm.error() << ", "
              << result_const.value() << " ± " << result_const.error() << ", "
              << result.value() << " ± " << result.error()
              << AnsiColor::normal << "\n";
    // check if integral is approximately equal to one
    COMPARE_ABSOLUTE_ERROR(result_no_norm.value(), 1., error_tolerance_no_norm);
    COMPARE_ABSOLUTE_ERROR(result_const.value(), 1., error_tolerance_const);
    COMPARE_ABSOLUTE_ERROR(result.value(), 1., result.error());
  }
}
