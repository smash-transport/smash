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
#include "../include/integrate.h"
#include "../include/resonances.h"

using namespace Smash;

TEST(spectral_functions) {
  Test::create_actual_particletypes();
  Test::create_actual_decaymodes();

  Integrator integrate;
  /* Upper limit for integration in GeV. Hadron masses are usually a few GeV,
   * so this should really be on the safe side. */
  const float max_mass = 100.;
  // error tolerance (max. deviation from one)
  const float error_tolerance = 0.07;

  /* Loop over all resonances. */
  for (const ParticleType &type : ParticleType::list_all()) {
    if (type.is_stable()) {
      continue;
    }
    /* Integrate spectral function. The factor 2m comes from dm^2 = 2m dm. */
    const auto result = integrate(type.minimum_mass(), max_mass,
                          [&](double m) {
                            const float m0 = type.mass();
                            const double width = type.total_width(m0);
                            return spectral_function(m, m0, width) * 2.*m;
                          });
    std::cout << type.name() << ": " << result.value() << " ± "
              << result.error() << "\n";
    // check if integral is approximately equal to one
    COMPARE_ABSOLUTE_ERROR(result.value(), 1., error_tolerance);
  }
}
