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
  const float max_mass = 100.;         // upper limit for integration
  const float error_tolerance = 0.07;  // error tolerance (max. deviation from one)

  /* Loop over all resonances. */
  for (const ParticleType &type : ParticleType::list_all()) {
    if (type.is_stable()) {
      continue;
    }
    /* Integrate spectral function. */
    const float min_mass = type.minimum_mass();
    const auto result = integrate(min_mass, max_mass,
                          [&](double m) {
                            const float m0 = type.mass();
                            const double width = type.total_width(m0);
                            return spectral_function(m, m0, width) * 2.*m;
                          });
    std::cout << type.pdgcode() << ": " << result.value() << " ± "
              << result.error() << "\n";
    // check if integral is approximately equal to one
    COMPARE_ABSOLUTE_ERROR(result.value(), 1., error_tolerance);
  }
}
