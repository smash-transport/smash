/*
 *    Copyright (c) 2013-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

/**
 * \file resonances.h
 * Functions related to resonance production.
 */

#ifndef SRC_INCLUDE_RESONANCES_H_
#define SRC_INCLUDE_RESONANCES_H_

#include "forwarddeclarations.h"
#include "particletype.h"

#include <cstddef>
#include <cstdint>

namespace Smash {

/** Parameters for spectral-function integration via GSL. */
struct IntegrandParameters {
  /// Type of the resonance
  const ParticleType *type;
  /// Mass of second particle
  double m2;
  /// Mandelstam s
  double s;
};

/**
 * Calculate Clebsch-Gordan coefficient
 *
 * \f$(-1)^{j_1 - j_2 + m_3} \sqrt(2 j_3 + 1) \cdot [Wigner 3J symbol] \f$
 * Note that the calculation assumes that the spin/isospin values (j/m)
 * have been multiplied by two (in order to be integer). */
double clebsch_gordan(const int j1, const int j2, const int j3,
                      const int m1, const int m2, const int m3);


/* Calculate isospin Clebsch-Gordan coefficient for two particles p1 and p2
 * coupling to a total isospin (I_tot, I_z). */
inline double isospin_clebsch_gordan(const ParticleType p1,
                                     const ParticleType p2,
                                     const int I_tot, const int I_z) {
  return clebsch_gordan (p1.isospin(), p2.isospin(), I_tot,
                         p1.isospin3(), p2.isospin3(), I_z);
}

/* Calculate isospin Clebsch-Gordan coefficient for two particles p1 and p2
 * coupling to a resonance Res. */
inline double isospin_clebsch_gordan(const ParticleType p1,
                                     const ParticleType p2,
                                     const ParticleType Res) {
  return clebsch_gordan (p1.isospin(), p2.isospin(), Res.isospin(),
                         p1.isospin3(), p2.isospin3(), Res.isospin3());
}


/**
 * Given the types of the two initial particles and a resonance,
 * return the 2-to-1 resonance production cross section.
 *
 * Checks are processed in the following order:
 * 1. Charge conservation
 * 2. Baryon number conservation
 * 3. Clebsch-Gordan
 * 4. Enough energy for all decay channels to be available for the resonance
 * 5. Detailed balance (reverse process exists)
 *
 * \param[in] type_particle1 Type information for the first initial particle.
 * \param[in] type_particle2 Type information for the second initial particle.
 * \param[in] type_resonance Type information for the resonance to be produced.
 * \param[in] mandelstam_s Mandelstam-s of the collision
 * of the two initial particles.
 * \param[in] cm_momentum_squared Square of the center-of-mass momentum of the
 * two initial particles.
 *
 * \return The cross section for the process
 * [initial particle 1] + [initial particle 2] -> resonance.
 */
double two_to_one_formation(const ParticleType &type_particle1,
                            const ParticleType &type_particle2,
                            const ParticleType &type_resonance,
                            double mandelstam_s, double cm_momentum_squared);


/**
 * Spectral function
 * \f$A(m)=\frac{1}{\pi}\frac{m\Gamma(m)}{(m^2-m_0^2)^2+(m\Gamma(m))^2}\f$
 * of the resonance.
 */
double spectral_function(double resonance_mass, double resonance_pole,
                         double resonance_width);

/**
 * Spectral function integrand for GSL integration.
 *
 * The integrand is \f$2m A(m) p_{cm}^f\f$, where \f$m\f$ is the
 * resonance mass, \f$A(m)\f$ is the spectral function
 *  and \f$p_{cm}^f\f$ is the center-of-mass momentum of the final state.
 *
 * \param[in] resonance_mass Actual mass of the resonance.
 * \param[in] parameters Container for the parameters needed
 * by the spectral function: Width of the resonance,
 * pole mass of the resonance, mass of the stable particle in the final state
 * and mandelstam-s of the process.
 */
double spectral_function_integrand(double resonance_mass, void * parameters);

/**
 * Resonance mass sampling for 2-particle final state
 * with *one resonance* and one *stable* particle.
 *
 * \param[in] type_resonance Type of the resonance particle.
 * \param[in] type_stable Type of the stable particle.
 * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
 *
 * \return The mass of the resonance particle.
 */
double sample_resonance_mass(const ParticleType &type_resonance,
                             const ParticleType &type_stable,
                             const float cms_energy);

/**
 * Function for 1-dimensional GSL integration.
 *
 * \param[in] integrand_function Function of 1 variable to be integrated over.
 * \param[in] parameters Container for possible parameters
 * needed by the integrand.
 * \param[in] lower_limit Lower limit of the integral.
 * \param[in] upper_limit Upper limit of the integral.
 * \param[out] integral_value Result of integration.
 * \param[out] integral_error Uncertainty of the result.
 */
void quadrature_1d(double (*integrand_function)(double, void *),
                          IntegrandParameters *parameters, double lower_limit,
                          double upper_limit, double *integral_value,
                          double *integral_error);

}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
