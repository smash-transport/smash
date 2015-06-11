/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_KINEMATICS_H_
#define SRC_INCLUDE_KINEMATICS_H_

#include "constants.h"
#include "experiment.h"

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
  const auto psqr = x * x * (T(0.25) / s) - mass_a_sqr;
  return psqr > 0. ? std::sqrt(psqr) : 0.;
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


/** Convert mandelstam-s to p_lab in a nucleon-nucleon collision.
 *
 * \fpPrecision Why \c double?
 */
inline double plab_from_s_NN(double mandelstam_s) {
  const double mNsqr = mN*mN;
  return std::sqrt((mandelstam_s - 2*mNsqr) * (mandelstam_s - 2*mNsqr)
                   - 4 * mNsqr * mNsqr) / (2 * mN);
}


/**
 * Convert E_kin to mandelstam-s for a fixed-target setup,
 * with a projectile of mass m_P and a kinetic energy e_kin
 * and a target of mass m_T at rest.
 *
 * \fpPrecision Why \c double?
 */
inline double s_from_Ekin(double e_kin, double m_P, double m_T) {
  return m_P*m_P + m_T*m_T + 2 * m_T * (m_P + e_kin);
}


/**
 * Convert p_lab to mandelstam-s for a fixed-target setup,
 * with a projectile of mass m_P and momentum plab
 * and a target of mass m_T at rest.
 *
 * \fpPrecision Why \c double?
 */
inline double s_from_plab(double plab, double m_P, double m_T) {
  return m_P*m_P + m_T*m_T + 2 * m_T * std::sqrt(m_P*m_P + plab*plab);
}

}  // namespace Smash

#endif  // SRC_INCLUDE_KINEMATICS_H_
