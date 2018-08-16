/*
 *
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "unittest.h"  // This include has to be first

#include "histogram.h"
#include "setup.h"

#include "../include/smash/formfactors.h"
#include "../include/smash/integrate.h"
#include "../include/smash/kinematics.h"
#include "../include/smash/stringfunctions.h"

using namespace smash;
using namespace UnitTest;

TEST(spectral_functions) {
  smash::Test::create_actual_particletypes();
  smash::Test::create_actual_decaymodes();

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
    const auto result_no_norm = integrate(0., 1., [&](double t) {
      return type.spectral_function_no_norm(type.min_mass_kinematic() +
                                            (1 - t) / t) /
             (t * t);
    });
    const auto result_const = integrate(0., 1., [&](double t) {
      return type.spectral_function_const_width((1 - t) / t) / (t * t);
    });
    const auto result = integrate(0., 1., [&](double t) {
      return type.spectral_function(type.min_mass_kinematic() + (1 - t) / t) /
             (t * t);
    });
    if (result_no_norm.value() > 1 + warning_level) {
      std::cout << AnsiColor::blue;
    } else if (result_no_norm.value() < 1 - warning_level) {
      std::cout << AnsiColor::yellow;
    }
    std::cout << utf8::fill_right(type.name(), 11) << ": "
              << format(result_no_norm.value(), nullptr, -1, 4) << " ± "
              << result_no_norm.error() << ", " << result_const.value() << " ± "
              << result_const.error() << ", " << result.value() << " ± "
              << result.error() << AnsiColor::normal << "\n";
    // check if integral is approximately equal to one
    COMPARE_ABSOLUTE_ERROR(result_no_norm.value(), 1., error_tolerance_no_norm);
    COMPARE_ABSOLUTE_ERROR(result_const.value(), 1., error_tolerance_const);
    COMPARE_ABSOLUTE_ERROR(result.value(), 1., 5 * result.error());
    //^ We use a bit higher tolerance, because the numerical algorithm might
    // underestimate
    //  the error.
  }
}

TEST(mass_sampling) {
  const ParticleType &res = ParticleType::find(0x12212);
  // Dummy reaction NN -> NN(1440) at sqrt(s) = 6 GeV
  const double sqrts = 6.0;
  const double mass_stable = 0.938;
  const int L = 0;
  const double dm_hist = 0.01;
  Histogram1d hist(dm_hist);
  // sample distribution and populate histogram
  const int N_sample = 1000000;
  hist.populate(N_sample, [&]() {
    return res.sample_resonance_mass(mass_stable, sqrts, L);
  });
  // hist.print_to_file("N1440_sampling.dat");
  hist.test([&](double m) {
    const double pcm = pCM(sqrts, mass_stable, m);
    const double bw = blatt_weisskopf_sqr(pcm, L);
    return res.spectral_function(m) * pcm * bw;
  });
}
