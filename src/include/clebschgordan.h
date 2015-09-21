/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_CLEBSCHGORDAN_H_
#define SRC_INCLUDE_CLEBSCHGORDAN_H_

#include "particletype.h"

namespace Smash {

/**
 * Calculate Clebsch-Gordan coefficient
 *
 * \f$(-1)^{j_a - j_b + m_c} \sqrt(2 j_c + 1) \cdot [Wigner 3J symbol] \f$
 * Note that the calculation assumes that the spin/isospin values (j/m)
 * have been multiplied by two (in order to be integer).
 */
float clebsch_gordan(const int j_a, const int j_b, const int j_c,
                     const int m_a, const int m_b, const int m_c);


/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for two particles
 * p_a and p_b coupling to a resonance Res.
 */
inline float isospin_clebsch_gordan_sqr_2to1(const ParticleType &p_a,
                                             const ParticleType &p_b,
                                             const ParticleType &Res) {
  const float cg = clebsch_gordan(p_a.isospin(), p_b.isospin(), Res.isospin(),
                                p_a.isospin3(), p_b.isospin3(), Res.isospin3());
  return cg*cg;
}


/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for three particles
 * p_a, p_b and p_c coupling to a total isospin (I_tot, I_z).
 *
 * Note that the coefficients also depend on the coupled isospin I_ab.
 */
inline float isospin_clebsch_gordan_sqr_3to1(const ParticleType &p_a,
                                             const ParticleType &p_b,
                                             const ParticleType &p_c,
                                             int I_tot, int I_z, int I_ab) {
    const int I_abz = p_a.isospin3() + p_b.isospin3();
    const float cg = clebsch_gordan(I_ab, p_c.isospin(), I_tot,
                                    I_abz, p_c.isospin3(), I_z)
                   * clebsch_gordan(p_a.isospin(), p_b.isospin(), I_ab,
                                    p_a.isospin3(), p_b.isospin3(), I_abz);
    return cg*cg;
}

/**
 * Calculate the squared isospin Clebsch-Gordan coefficient for a
 * 2-to-2 reaction A + B -> C + D.
 */
float isospin_clebsch_gordan_sqr_2to2(const ParticleType &t_a,
                                      const ParticleType &t_b,
                                      const ParticleType &t_c,
                                      const ParticleType &t_d);


}  // namespace Smash

#endif  // SRC_INCLUDE_CLEBSCHGORDAN_H_
