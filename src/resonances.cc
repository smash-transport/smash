/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */
#include "include/resonances.h"

#include <gsl/gsl_sf_coupling.h>

#include "include/formfactors.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/random.h"

namespace Smash {


double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c) {
  const auto &log = logger<LogArea::Resonances>();
  double wigner_3j =  gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  double result = 0.;
  if (std::abs(wigner_3j) > really_small)
    result = std::pow(-1, (j_a-j_b+m_c)/2.) * std::sqrt(j_c + 1) * wigner_3j;

  log.debug("CG: ", result, " I1: ", j_a, " I2: ", j_b, " IR: ", j_c,
            " iz1: ", m_a, " iz2: ", m_b, " izR: ", m_c);

  return result;
}


/* Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a total isospin (I_tot, I_z).
 *
 * \fpPrecision Why \c double?
 */
static double isospin_clebsch_gordan(const ParticleType &p_a,
                                     const ParticleType &p_b,
                                     const int I_tot, const int I_z) {
  return clebsch_gordan(p_a.isospin(), p_b.isospin(), I_tot,
                        p_a.isospin3(), p_b.isospin3(), I_z);
}


double isospin_clebsch_gordan(const ParticleType &t_a, const ParticleType &t_b,
                            const ParticleType &t_c, const ParticleType &t_d) {
  const int I_z = t_a.isospin3() + t_b.isospin3();

  /* Compute total isospin range with given initial and final particles. */
  const int I_max = std::min(t_a.isospin() + t_b.isospin(),
                             t_c.isospin() + t_d.isospin());
  int I_min = std::max(std::abs(t_a.isospin() - t_b.isospin()),
                       std::abs(t_c.isospin() - t_d.isospin()));
  I_min = std::max(I_min, std::abs(I_z));

  /* Loop over total isospin in allowed range.
  * Use decrement of 2, since isospin is multiplied by 2. */
  double isospin_factor = 0.;
  for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
    isospin_factor = isospin_factor +
        isospin_clebsch_gordan(t_c, t_d, I_tot, I_z)
      * isospin_clebsch_gordan(t_a, t_b, I_tot, I_z);
  }
  return isospin_factor;
}


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
                            const float mass_stable, const float cms_energy,
                            int L) {
  /* largest possible mass: Use 'nextafter' to make sure it is not above the
   * physical limit by numerical error. */
  const float max_mass = std::nextafter(cms_energy - mass_stable, 0.f);
  // largest possible cm momentum (from smallest mass)
  const float pcm_max = pCM(cms_energy, mass_stable, type_res.minimum_mass());
  const float blw_max = pcm_max * blatt_weisskopf_sqr(pcm_max, L);
  /* The maximum of the spectral-function ratio 'usually' happens at the
   * largest mass. However, this is not always the case, therefore we need
   * an additional fudge factor (empirically 3.6 happens to be sufficient). */
  const float q_max = type_res.spectral_function(max_mass)
                    / type_res.spectral_function_simple(max_mass) * 3.6;
  const float max = blw_max * q_max;  // maximum value for rejection sampling
  float mass_res, val;
  // Loop: rejection sampling
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
    log.fatal("maximum not correct in sample_resonance_mass: ",
              val, " ", max, " ", type_res.pdgcode(), " ",
              mass_stable, " ", cms_energy, " ", mass_res);
    throw std::runtime_error("Maximum not correct in sample_resonance_mass!");
  }

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
  constexpr float q_max = 14.;        // this value is determined empirically
  const float max = blw_max * q_max;  // maximum value for rejection sampling
  float mass_1, mass_2, val;
  // Loop: rejection sampling
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
    log.error("maximum not correct in sample_resonance_masses: ",
              val, " ", max, " ", t1.pdgcode(), " ", t2.pdgcode(), " ",
              cms_energy, " ", mass_1, " ", mass_2);
  }

  return {mass_1, mass_2};
}


}  // namespace Smash
