/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_KINEMATICS_H_
#define SRC_INCLUDE_KINEMATICS_H_

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
 * Return the squared center-of-mass momentum of two particles,
 * given sqrt(s) and their masses.
 *
 * \param srts sqrt(s) of the process [GeV].
 * \param mass_a Mass of first particle [GeV].
 * \param mass_b Mass of second particle [GeV].
 */
template <typename T>
T pCM_sqr(const T srts, const T mass_a, const T mass_b) noexcept {
  const auto s = srts * srts;
  const auto mass_a_sqr = mass_a * mass_a;
  const auto x = s + mass_a_sqr - mass_b * mass_b;
  return x * x * (T(0.25) / s) - mass_a_sqr;
}


/// nucleon mass
const double mN = 0.938;


/** Convert mandelstam-s to p_lab in a nucleon-nucleon collision. */
inline double plab_from_s_NN(double mandelstam_s) {
  const double mNsqr = mN*mN;
  return std::sqrt((mandelstam_s - 2*mNsqr) * (mandelstam_s - 2*mNsqr)
                   - 4 * mNsqr * mNsqr) / (2 * mN);
}


}  // namespace Smash

#endif  // SRC_INCLUDE_KINEMATICS_H_
