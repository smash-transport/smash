/*
 *    Copyright (c) 2013-2020,2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/clebschgordan.h"

#include <numeric>

#include "gsl/gsl_sf_coupling.h"

#include "smash/clebschgordan_lookup.h"
#include "smash/constants.h"
#include "smash/logging.h"

namespace smash {

/**
 * Calculate isospin Clebsch-Gordan coefficient for two particles p_a and p_b
 * coupling to a total isospin \see ClebschGordan::coefficient for details
 * (I_tot, I_z).
 *
 * \param[in] p_a Information of particle type for first particle
 * \param[in] p_b Information of particle type for second particle
 * \param[out] I_tot Total isospin of the reaction
 * \param[out] I_z Total isospin 3 component of the reaction
 */
static double isospin_clebsch_gordan_2to1(const ParticleType &p_a,
                                          const ParticleType &p_b,
                                          const int I_tot, const int I_z) {
  return ClebschGordan::coefficient(p_a.isospin(), p_b.isospin(), I_tot,
                                    p_a.isospin3(), p_b.isospin3(), I_z);
}

double isospin_clebsch_gordan_sqr_3to1(const ParticleType &p_a,
                                       const ParticleType &p_b,
                                       const ParticleType &p_c,
                                       const ParticleType &Res) {
  // Calculate allowed isospin range for 3->1 reaction I_ab
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
  const double cg =
      ClebschGordan::coefficient(I_ab, p_c.isospin(), Res.isospin(), I_abz,
                                 p_c.isospin3(), Res.isospin3()) *
      ClebschGordan::coefficient(p_a.isospin(), p_b.isospin(), I_ab,
                                 p_a.isospin3(), p_b.isospin3(), I_abz);
  return cg * cg;
}

double isospin_clebsch_gordan_sqr_2to2(const ParticleType &p_a,
                                       const ParticleType &p_b,
                                       const ParticleType &p_c,
                                       const ParticleType &p_d, const int I) {
  const int I_z = p_a.isospin3() + p_b.isospin3();

  /* Loop over total isospin in allowed range. */
  double isospin_factor = 0.;
  for (const int I_tot : I_tot_range(p_a, p_b, p_c, p_d)) {
    if (I < 0 || I_tot == I) {
      const double cg_in = isospin_clebsch_gordan_2to1(p_a, p_b, I_tot, I_z);
      const double cg_out = isospin_clebsch_gordan_2to1(p_c, p_d, I_tot, I_z);
      isospin_factor = isospin_factor + cg_in * cg_in * cg_out * cg_out;
    }
  }
  return isospin_factor;
}

}  // namespace smash
