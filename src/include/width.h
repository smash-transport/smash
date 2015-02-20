/*
 *    Copyright (c) 2012-2014
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_WIDTH_H_
#define SRC_INCLUDE_WIDTH_H_

#include <cmath>

#include "constants.h"
#include "particletype.h"

namespace Smash {

const float interaction_radius = 1. / hbarc;

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
 * Returns the squared Blatt-Weisskopf functions,
 * which govern the mass dependence of the width of a resonance
 * decaying into two particles A and B. See e.g. Effenberger's thesis, page 28.
 *
 * \param x = p_ab*R  with
 *        p_ab = relative momentum of outgoing particles AB and
 *        R = interaction radius
 * \param L Angular momentum of outgoing particles AB.
 */
float BlattWeisskopf(const float x, const int L)
#ifdef NDEBUG
    noexcept
#endif
;

double Post_FF_sqr (double m, double M0, double s0, double L);

}  // namespace Smash

#endif  // SRC_INCLUDE_WIDTH_H_
