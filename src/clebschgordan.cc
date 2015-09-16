/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/clebschgordan.h"

#include <gsl/gsl_sf_coupling.h>

#include "include/constants.h"
#include "include/logging.h"

namespace Smash {


float clebsch_gordan(const int j_a, const int j_b, const int j_c,
                      const int m_a, const int m_b, const int m_c) {
  const double wigner_3j =  gsl_sf_coupling_3j(j_a, j_b, j_c, m_a, m_b, -m_c);
  if (std::abs(wigner_3j) < really_small)
    return 0.;

  const float result = std::pow(-1, (j_a-j_b+m_c)/2.)
                      * std::sqrt(j_c + 1) * wigner_3j;

#ifndef NDEBUG
  const auto &log = logger<LogArea::Resonances>();
  log.debug("CG: ", result, " I1: ", j_a, " I2: ", j_b, " IR: ", j_c,
            " iz1: ", m_a, " iz2: ", m_b, " izR: ", m_c);
#endif

  return result;
}


/* Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a total isospin (I_tot, I_z).
 */
static float isospin_clebsch_gordan_2to1(const ParticleType &p_a,
                                         const ParticleType &p_b,
                                         const int I_tot, const int I_z) {
  return clebsch_gordan(p_a.isospin(), p_b.isospin(), I_tot,
                        p_a.isospin3(), p_b.isospin3(), I_z);
}


float isospin_clebsch_gordan_2to2(const ParticleType &t_a,
                                  const ParticleType &t_b,
                                  const ParticleType &t_c,
                                  const ParticleType &t_d) {
  const int I_z = t_a.isospin3() + t_b.isospin3();

  /* Compute total isospin range with given initial and final particles. */
  const int I_max = std::min(t_a.isospin() + t_b.isospin(),
                             t_c.isospin() + t_d.isospin());
  int I_min = std::max(std::abs(t_a.isospin() - t_b.isospin()),
                       std::abs(t_c.isospin() - t_d.isospin()));
  I_min = std::max(I_min, std::abs(I_z));

  /* Loop over total isospin in allowed range.
   * Use decrement of 2, since isospin is multiplied by 2. */
  float isospin_factor = 0.;
  for (int I_tot = I_max; I_tot >= I_min; I_tot -= 2) {
    isospin_factor = isospin_factor +
        isospin_clebsch_gordan_2to1(t_c, t_d, I_tot, I_z)
      * isospin_clebsch_gordan_2to1(t_a, t_b, I_tot, I_z);
  }
  return isospin_factor;
}


}  // namespace Smash
