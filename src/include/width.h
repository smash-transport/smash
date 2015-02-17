/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_WIDTH_H_
#define SRC_INCLUDE_WIDTH_H_

#include <cmath>

#include "particletype.h"

namespace Smash {

/**
 * Return the center-of-mass momentum of two particles,
 * given sqrt(s) and their masses.
 *
 * \param srts sqrt(s) of the process [GeV].
 * \param mass_a Mass of first particle [GeV].
 * \param mass_b Mass of second particle [GeV].
 */
template <typename T>
T pCM(const T srts, const T mass_a, const T mass_b) noexcept {
  const auto s = srts * srts;
  const auto mass_a_sqr = mass_a * mass_a;
  const auto x = s + mass_a_sqr - mass_b * mass_b;
  return std::sqrt(x * x * (T(0.25) / s) - mass_a_sqr);
}


/**
 * Get the mass-dependent width of a two-body decay into stable particles
 * according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
 *
 * \param mass Actual mass of the decaying particle [GeV].
 * \param poleMass Pole mass of the decaying particle [GeV].
 * \param mass_a Mass of the first daughter particle [GeV].
 * \param mass_b Mass of the second daughter particle [GeV].
 * \param L Angular momentum of the decay.
 * \param partial_width_at_pole Partial width at the pole mass [GeV].
 */
float width_Manley_stable(const float mass, const float poleMass,
                          const float mass_a, const float mass_b,
                          const int L, const float partial_width_at_pole);

/**
 * Get the mass-dependent width of a two-body decay into one stable and one
 * unstable particle according to Manley/Saleski, Phys. Rev. D 45 (1992) 4002.
 *
 * \param mass Actual mass of the decaying particle [GeV].
 * \param poleMass Pole mass of the decaying particle [GeV].
 * \param mass_stable Mass of the stable daughter particle [GeV].
 * \param type_unstable Type of the unstable daughter particle [GeV].
 * \param L Angular momentum of the decay.
 * \param partial_width_at_pole Partial width at the pole mass [GeV].
 */
float width_Manley_semistable(const float mass, const float poleMass,
                              const float mass_stable,
                              const ParticleType &type_unstable,
                              const int L, const float partial_width_at_pole);

/**
 * Get the mass-dependent in-width for a resonance formation process from one
 * stable and one unstable particle according to Manley/Saleski (PRD45),
 * see also PhD thesis Effenberger, eq. (2.77).
 *
 * \param mass Actual mass of the produced resonance [GeV].
 * \param poleMass Pole mass of the produced resonance [GeV].
 * \param mass_stable Mass of the stable incoming particle [GeV].
 * \param mass_unstable Mass of the unstable incoming particle [GeV].
 * \param type_unstable Type of the unstable incoming particle [GeV].
 * \param L Angular momentum of the corresponding decay.
 * \param partial_width_at_pole Partial width at the pole mass [GeV].
 */
float in_width_Manley_semistable(const float mass, const float poleMass,
                                 const float mass_stable,
                                 const float mass_unstable,
                                 const ParticleType &type_unstable,
                                 const int L,
                                 const float partial_width_at_pole);

}  // namespace Smash

#endif  // SRC_INCLUDE_WIDTH_H_
