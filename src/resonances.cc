/*
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/resonances.h"

#include "include/formfactors.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/random.h"

namespace Smash {


/* Integrand for spectral-function integration with one resonance. */
float spec_func_integrand_1res(float resonance_mass, float sqrts,
                               float stable_mass, const ParticleType &type) {
  if (sqrts <= stable_mass + resonance_mass) {
    return 0.;
  }

  /* Integrand is the spectral function weighted by the CM momentum of the
   * final state. */
  return type.spectral_function(resonance_mass)
       * pCM(sqrts, stable_mass, resonance_mass);
}


/* Integrand for spectral-function integration with two resonances. */
float spec_func_integrand_2res(float sqrts, float res_mass_1, float res_mass_2,
                               const ParticleType &t1, const ParticleType &t2) {
  if (sqrts <= res_mass_1 + res_mass_2) {
    return 0.;
  }

  /* Integrand is the product of the spectral function weighted by the
   * CM momentum of the final state. */
  return t1.spectral_function(res_mass_1)
       * t2.spectral_function(res_mass_2)
       * pCM(sqrts, res_mass_1, res_mass_2);
}


/* Resonance mass sampling for 2-particle final state */
float sample_resonance_mass(const ParticleType &type_res,
                            const float mass_stable,
                            const float cms_energy, int L) {
  /* largest possible mass: Use 'nextafter' to make sure it is not above the
   * physical limit by numerical error. */
  const float max_mass = std::nextafter(cms_energy - mass_stable, 0.f);
  // largest possible cm momentum (from smallest mass)
  const float pcm_max = pCM(cms_energy, mass_stable, type_res.minimum_mass());
  const float blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);

  float mass_res, val;
  // outer loop: repeat if maximum is too small
  do {
    /* The maximum of the spectral-function ratio 'usually' happens at the
    * largest mass. However, this is not always the case, therefore we need
    * an additional fudge factor (purely empirical). */
    const float q_max = type_res.spectral_function(max_mass)
                      / type_res.spectral_function_simple(max_mass)
                      * type_res.max_factor1();
    const float max = blw_max * q_max;  // maximum value for rejection sampling
    // inner loop: rejection sampling
    do {
      // sample mass from a simple Breit-Wigner (aka Cauchy) distribution
      mass_res = Random::cauchy(type_res.mass(), type_res.width_at_pole()/2.f,
                                type_res.minimum_mass(), max_mass);
      // determine cm momentum for this case
      const float pcm = pCM(cms_energy, mass_stable, mass_res);
      const float blw = pcm * blatt_weisskopf_sqr(pcm, L);
      // determine ratio of full to simple spectral function
      const float q = type_res.spectral_function(mass_res)
                    / type_res.spectral_function_simple(mass_res);
      val = q * blw;
    } while (val < Random::uniform(0.f, max));

    // check that we are using the proper maximum value
    if (val > max) {
      const auto &log = logger<LogArea::Resonances>();
      log.debug("maximum is being increased in sample_resonance_mass: ",
                type_res.max_factor1(), " ", val/max, " ", type_res.pdgcode(),
                " ", mass_stable, " ", cms_energy, " ", mass_res);
      type_res.increase_max_factor1(val/max);
    } else {
      break;  // maximum ok, exit loop
    }
  } while (true);

  return mass_res;
}


/* Resonance mass sampling for 2-particle final state with two resonances. */
std::pair<float, float> sample_resonance_masses(const ParticleType &t1,
                                                const ParticleType &t2,
                                                const float cms_energy, int L) {
  /* Sample resonance mass from the distribution
   * used for calculating the cross section. */
  const float max_mass_1 = std::nextafter(cms_energy - t2.minimum_mass(), 0.f);
  const float max_mass_2 = std::nextafter(cms_energy - t1.minimum_mass(), 0.f);
  // largest possible cm momentum (from smallest mass)
  const float pcm_max = pCM(cms_energy, t1.minimum_mass(), t2.minimum_mass());
  const float blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);

  float mass_1, mass_2, val;
  // outer loop: repeat if maximum is too small
  do {
    // maximum value for rejection sampling (determined empirically)
    const float max = blw_max * t1.max_factor2();
    // inner loop: rejection sampling
    do {
      // sample mass from a simple Breit-Wigner (aka Cauchy) distribution
      mass_1 = Random::cauchy(t1.mass(), t1.width_at_pole()/2.f,
                              t1.minimum_mass(), max_mass_1);
      mass_2 = Random::cauchy(t2.mass(), t2.width_at_pole()/2.f,
                              t2.minimum_mass(), max_mass_2);
      // determine cm momentum for this case
      const float pcm = pCM(cms_energy, mass_1, mass_2);
      const float blw = pcm * blatt_weisskopf_sqr(pcm, L);
      // determine ratios of full to simple spectral function
      const float q1 = t1.spectral_function(mass_1)
                    / t1.spectral_function_simple(mass_1);
      const float q2 = t2.spectral_function(mass_2)
                    / t2.spectral_function_simple(mass_2);
      val = q1 * q2 * blw;
    } while (val < Random::uniform(0.f, max));

    if (val > max) {
      const auto &log = logger<LogArea::Resonances>();
      log.debug("maximum is being increased in sample_resonance_masses: ",
                t1.max_factor2(), " ", val/max, " ", t1.pdgcode(), " ",
                t2.pdgcode(), " ", cms_energy, " ", mass_1, " ", mass_2);
      t1.increase_max_factor2(val/max);
    } else {
      break;  // maximum ok, exit loop
    }
  } while (true);

  return {mass_1, mass_2};
}


}  // namespace Smash
