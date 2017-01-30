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
  // error tolerance (max. deviation from one)
  const double error_tolerance_no_norm = 0.55;
  const double warning_level = 0.21;
  // error tolerance (max. deviation from one) for constant-width SF
  const double error_tolerance_const = 0.0032;

  /* Loop over all resonances. */
  for (const ParticleType &type : ParticleType::list_all()) {
    if (type.is_stable()) {
      continue;
    }
    /* Integrate spectral function.
     * We transform the integrals using m = m_min + (1 - t)/t to make them
     * definite and to avoid numerical problems. */
    const auto result_no_norm = integrate(0., 1.,
                          [&](double t) {
                            return type.spectral_function_no_norm(type.minimum_mass() + (1 - t)/t) / (t*t);
                          });
    const auto result_const = integrate(0., 1.,
                          [&](double t) {
                            return type.spectral_function_const_width((1 - t)/t) / (t*t);
                          });
    const auto result = integrate(0., 1.,
                          [&](double t) {
                            return type.spectral_function(type.minimum_mass() + (1 - t)/t) / (t*t);
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
    COMPARE_ABSOLUTE_ERROR(result.value(), 1., 5*result.error());
    //^ We use a bit higher tolerance, because the numerical algorithm might underestimate
    //  the error.
  }
}
