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


/**
 * Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a resonance Res.
 *
 * \fpPrecision Why \c double?
 */
inline double isospin_clebsch_gordan(const ParticleType &p_a,
                                     const ParticleType &p_b,
                                     const ParticleType &Res) {
  return clebsch_gordan(p_a.isospin(), p_b.isospin(), Res.isospin(),
                        p_a.isospin3(), p_b.isospin3(), Res.isospin3());
}


/**
 * Calculate isospin Clebsch-Gordan coefficient for three particles p_a, p_b and
 * p_c coupling to a total isospin (I_tot, I_z) Res.
 *
 * Note that the coefficients also depend on the coupled isospin I_ab.
 *
 * \fpPrecision Why \c double?
 */
inline double isospin_clebsch_gordan(const ParticleType &p_a,
                                     const ParticleType &p_b,
                                     const ParticleType &p_c,
                                     int I_tot, int I_z, int I_ab) {
    const int I_abz = p_a.isospin3() + p_b.isospin3();
    return clebsch_gordan(I_ab, p_c.isospin(), I_tot,
                          I_abz, p_c.isospin3(), I_z)
         * clebsch_gordan(p_a.isospin(), p_b.isospin(), I_ab,
                          p_a.isospin3(), p_b.isospin3(), I_abz);
}

/**
 * Calculate isospin Clebsch-Gordan coefficient for a 2-to-2 reaction
 * A + B -> C + D.
 *
 * \fpPrecision Why \c double?
 */
double isospin_clebsch_gordan(const ParticleType &t_a, const ParticleType &t_b,
                              const ParticleType &t_c, const ParticleType &t_d);


/**
 * Spectral function integrand for GSL integration, with one resonance in the
 * final state (the second particle is stable).
 *
 * The integrand is \f$ A(m) p_{cm}^f \f$, where \f$ m \f$ is the
 * resonance mass, \f$ A(m) \f$ is the spectral function
 *  and \f$ p_{cm}^f \f$ is the center-of-mass momentum of the final state.
 *
 * \param[in] resonance_mass Actual mass of the resonance.
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] stable_mass mass of the stable particle in the final state
 * \param[in] type type of the resonance
 */
float spec_func_integrand_1res(float resonance_mass, float sqrts,
                               float stable_mass, const ParticleType &type);


/**
 * Spectral function integrand for GSL integration, with two resonances in the
 * final state.
 *
 * The integrand is \f$ A_1(m_1) A_2(m_2) p_{cm}^f \f$, where \f$ m_1 \f$ and
 * \f$ m_2 \f$ are the resonance masses, \f$ A_1 \f$ and \f$ A_2 \f$ are the
 * spectral functions and \f$ p_{cm}^f \f$ is the center-of-mass momentum of
 * the final state.
 *
 * \param[in] sqrts Center-of-mass energy, i.e. sqrt of Mandelstam s.
 * \param[in] res_mass_1 Actual mass of the first resonance.
 * \param[in] res_mass_2 Actual mass of the second resonance.
 * \param[in] t1 Type of the first resonance.
 * \param[in] t2 Type of the second resonance.
 */
float spec_func_integrand_2res(float sqrts, float res_mass_1, float res_mass_2,
                               const ParticleType &t1, const ParticleType &t2);


/**
 * Resonance mass sampling for 2-particle final state
 * with *one resonance* and one *stable* particle.
 *
 * \param[in] type_res Type of the resonance particle.
 * \param[in] mass_stable Mass of the stable particle.
 * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
 * \param[in] L relative angular momentum of the final-state particles
 *
 * \return The mass of the resonance particle.
 */
float sample_resonance_mass(const ParticleType &type_res,
                            const float mass_stable, const float cms_energy,
                            int L = 0);

/**
 * Resonance mass sampling for 2-particle final state with two resonances.
 *
 * \param[in] type_res_1 Type of the first resonance.
 * \param[in] type_res_2 Type of the second resonance.
 * \param[in] cms_energy center-of-mass energy of the 2-particle final state.
 * \param[in] L relative angular momentum of the final-state particles
 *
 * \return The masses of the resonance particles.
 */
std::pair<float, float> sample_resonance_masses(const ParticleType &type_res_1,
                                                const ParticleType &type_res_2,
                                                const float cms_energy,
                                                int L = 0);


}  // namespace Smash

#endif  // SRC_INCLUDE_RESONANCES_H_
