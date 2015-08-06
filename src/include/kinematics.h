/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_KINEMATICS_H_
#define SRC_INCLUDE_KINEMATICS_H_

#include <array>

#include "constants.h"

namespace Smash {

/**
 * Return velocity in the center of velocities frame of two particles given
 * their mandelstam s and masses
 *
 * \param mandelstamm s of the collision [GeV^2]
 * \param ma Mass of the first particle [GeV]
 * \param ma Mass of the second particle [GeV]
 */
inline double center_of_velocity_v(float s, float ma, float mb) {
  const float m_sum = ma + mb;
  const float m_dif = ma - mb;
  return std::sqrt((s - m_sum*m_sum) / (s - m_dif*m_dif));
}

/**
 * Return velocity of projectile in the fixed target frame given
 * their mandelstam s of projectile and target and their masses
 *
 * \param mandelstamm s of the collision [GeV^2]
 * \param ma Mass of the projectile [GeV]
 * \param ma Mass of the target [GeV]
 */
inline double fixed_target_projectile_v(float s, float ma, float mb) {
  const float inv_gamma = 2 * ma * mb / (s - ma * ma - mb * mb);
  return std::sqrt(1.0 - inv_gamma*inv_gamma);
}

/**
 * Return the squared center-of-mass momentum of two particles,
 * given s and their masses.
 *
 * \param s mandelstamm s of the process [GeV^2].
 * \param mass_a Mass of first particle [GeV].
 * \param mass_b Mass of second particle [GeV].
 */
template <typename T>
T pCM_sqr_from_s(const T s, const T mass_a, const T mass_b) noexcept {
  const auto mass_a_sqr = mass_a * mass_a;
  const auto x = s + mass_a_sqr - mass_b * mass_b;
  return x * x * (T(0.25) / s) - mass_a_sqr;
}

/**
 * Return the center-of-mass momentum of two particles,
 * given s and their masses.
 *
 * \param s mandelstamm s of the process [GeV^2].
 * \param mass_a Mass of first particle [GeV].
 * \param mass_b Mass of second particle [GeV].
 */
template <typename T>
T pCM_from_s(const T s, const T mass_a, const T mass_b) noexcept {
  const auto psqr = pCM_sqr_from_s(s, mass_a, mass_b);
  return psqr > T(0.) ? std::sqrt(psqr) : T(0.);
}

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
  const auto psqr = pCM_sqr_from_s(srts*srts, mass_a, mass_b);
  return psqr > T(0.) ? std::sqrt(psqr) : T(0.);
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
  return pCM_sqr_from_s(srts * srts, mass_a, mass_b);
}

/**
 * Get the range of mandelstam-t values allowed in a particular 2->2 process,
 * see PDG 2014 booklet, equ. (46.34).
 * \param srts sqrt(s) of the process [GeV].
 * \param m1 Mass of first  incoming particle [GeV].
 * \param m2 Mass of second incoming particle [GeV].
 * \param m3 Mass of first  outgoing particle [GeV].
 * \param m4 Mass of second outgoing particle [GeV].
 * \return array consisiting of {t_min, t_max}
 * Note that both t_min and t_max are negative,
 * with |t_min| < |t_max|, i.e. t_min > t_max.
 */
template <typename T>
std::array<T, 2> get_t_range(const T srts, const T m1, const T m2,
                                           const T m3, const T m4) {
  const T p_i = pCM(srts, m1, m2);  // initial-state CM momentum
  const T p_f = pCM(srts, m3, m4);  // final-state CM momentum
  const T sqrt_t0 = (m1*m1 - m2*m2 - m3*m3 + m4*m4) / (2.*srts);
  const T t0 = sqrt_t0 * sqrt_t0;
  const T t_min = t0 - (p_i-p_f)*(p_i-p_f),
          t_max = t0 - (p_i+p_f)*(p_i+p_f);
  return {t_min, t_max};
}


/** Convert mandelstam-s to p_lab in a nucleon-nucleon collision.
 *
 * \fpPrecision Why \c double?
 */
inline double plab_from_s_NN(double s_NN) {
  const double tmp = s_NN * (s_NN - 4 * nucleon_mass * nucleon_mass);
  return std::sqrt(tmp) / (2 * nucleon_mass);
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
