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

#include <cstddef>
#include <cstdint>

#include "forwarddeclarations.h"
#include "particletype.h"

namespace Smash {

/**
 * Calculate Clebsch-Gordan coefficient
 *
 * \f$(-1)^{j_a - j_b + m_c} \sqrt(2 j_c + 1) \cdot [Wigner 3J symbol] \f$
 * Note that the calculation assumes that the spin/isospin values (j/m)
 * have been multiplied by two (in order to be integer).
 *
 * \fpPrecision Why \c double?
 */
double clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c);


/* Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a total isospin (I_tot, I_z).
 *
 * \fpPrecision Why \c double?
 */
inline double isospin_clebsch_gordan(const ParticleType &p_a,
                                     const ParticleType &p_b,
                                     const int I_tot, const int I_z) {
  return clebsch_gordan (p_a.isospin(), p_b.isospin(), I_tot,
                         p_a.isospin3(), p_b.isospin3(), I_z);
}

/* Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a resonance Res.
 *
 * \fpPrecision Why \c double?
 */
inline double isospin_clebsch_gordan(const ParticleType &p_a,
                                     const ParticleType &p_b,
                                     const ParticleType &Res) {
  return clebsch_gordan (p_a.isospin(), p_b.isospin(), Res.isospin(),
                         p_a.isospin3(), p_b.isospin3(), Res.isospin3());
}


/**
 * Spectral function
 * \f$A(m)=\frac{1}{\pi}\frac{m\Gamma(m)}{(m^2-m_0^2)^2+(m\Gamma(m))^2}\f$
 * of the resonance.
 *
 * \fpPrecision Why \c double?
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
 * \param[in] srts sqrt(s) of the process
 * \param[in] stable_mass mass of the stable particle in the final state
 * \param[in] type type of the resonance
 *
 * \fpPrecision Why \c double?
 */
double spectral_function_integrand(double resonance_mass, double srts,
                                   double stable_mass,
                                   const ParticleType &type);

/**
 * Resonance mass sampling for 2-particle final state
 * with *one resonance* and one *stable* particle.
 *
 * \param[in] type_resonance Type of the resonance particle.
 * \param[in] mass_stable Mass of the stable particle.
 * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
 *
 * \return The mass of the resonance particle.
 *
 * \fpPrecision Why \c double?
 */
float sample_resonance_mass(const ParticleType &type_resonance,
                            const float mass_stable, const double cms_energy);

}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
