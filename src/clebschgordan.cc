/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/clebschgordan.h"

#include <gsl/gsl_sf_coupling.h>
#include <numeric>
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


float isospin_clebsch_gordan_sqr_3to1(const ParticleType &p_a,
                                      const ParticleType &p_b,
                                      const ParticleType &p_c,
                                      const ParticleType &Res) {
  // calculate allowed I_ab
  const auto min_I_ab = std::abs(p_a.isospin() - p_b.isospin());
  const auto max_I_ab = p_a.isospin() + p_b.isospin();
  std::vector<int> possible_I_ab(max_I_ab - min_I_ab + 1);
  std::iota(possible_I_ab.begin(), possible_I_ab.end(), min_I_ab);
  std::vector<int> allowed_I_ab;
  allowed_I_ab.reserve(possible_I_ab.size());
  for (const auto Iab : possible_I_ab) {
    const auto min_I = std::abs(Iab - p_c.isospin());
    const auto max_I = Iab + p_c.isospin();
    if (min_I <= Res.isospin() && Res.isospin() <= max_I) {
      allowed_I_ab.push_back(Iab);
    }
  }
  if (allowed_I_ab.size() != 1) {
    throw std::runtime_error(
        "The coupled 3-body isospin state is not uniquely defined for " +
        Res.name() + " -> " + p_a.name() + " " + p_b.name() + " " + p_c.name());
  }
  const auto I_ab = allowed_I_ab[0];

  const int I_abz = p_a.isospin3() + p_b.isospin3();
  const float cg = clebsch_gordan(I_ab, p_c.isospin(), Res.isospin(),
                                  I_abz, p_c.isospin3(), Res.isospin3())
                  * clebsch_gordan(p_a.isospin(), p_b.isospin(), I_ab,
                                  p_a.isospin3(), p_b.isospin3(), I_abz);
  return cg*cg;
}


float isospin_clebsch_gordan_sqr_2to2(const ParticleType &t_a,
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
    const float cg_in = isospin_clebsch_gordan_2to1(t_a, t_b, I_tot, I_z);
    const float cg_out = isospin_clebsch_gordan_2to1(t_c, t_d, I_tot, I_z);
    isospin_factor = isospin_factor + cg_in*cg_in * cg_out*cg_out;
  }
  return isospin_factor;
}


}  // namespace Smash
