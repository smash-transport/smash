/*
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_KINEMATICS_H_
#define SRC_INCLUDE_KINEMATICS_H_

#include <array>
#include <sstream>

#include "constants.h"

namespace smash {

/**
 * Return velocity in the center of velocities frame of two particles given
 * their mandelstam s and masses
 *
 * \param s mandelstamm s of the collision [GeV^2]
 * \param ma Mass of the first particle [GeV]
 * \param mb Mass of the second particle [GeV]
 *
 * needs to be double to allow for calculations at LHC energies
 */
inline double center_of_velocity_v(double s, double ma, double mb) {
  const double m_sum = ma + mb;
  const double m_dif = ma - mb;
  return std::sqrt((s - m_sum * m_sum) / (s - m_dif * m_dif));
}

/**
 * Return velocity of projectile in the fixed target frame given
 * their mandelstam s of projectile and target and their masses
 *
 * \param s mandelstamm s of the collision [GeV^2]
 * \param ma Mass of the projectile [GeV]
 * \param mb Mass of the target [GeV]
 */
inline double fixed_target_projectile_v(double s, double ma, double mb) {
  const double inv_gamma = 2 * ma * mb / (s - ma * ma - mb * mb);
  return std::sqrt(1.0 - inv_gamma * inv_gamma);
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
 * \param sqrts sqrt(s) of the process [GeV].
 * \param mass_a Mass of first particle [GeV].
 * \param mass_b Mass of second particle [GeV].
 */
template <typename T>
T pCM(const T sqrts, const T mass_a, const T mass_b) noexcept {
  return pCM_from_s(sqrts * sqrts, mass_a, mass_b);
}

/**
 * Return the squared center-of-mass momentum of two particles,
 * given sqrt(s) and their masses.
 *
 * \param sqrts sqrt(s) of the process [GeV].
 * \param mass_a Mass of first particle [GeV].
 * \param mass_b Mass of second particle [GeV].
 */
template <typename T>
T pCM_sqr(const T sqrts, const T mass_a, const T mass_b) noexcept {
  return pCM_sqr_from_s(sqrts * sqrts, mass_a, mass_b);
}

/**
 * Get the range of mandelstam-t values allowed in a particular 2->2 process,
 * see PDG 2014 booklet, equ. (46.34).
 * \param sqrts sqrt(s) of the process [GeV].
 * \param m1 Mass of first  incoming particle [GeV].
 * \param m2 Mass of second incoming particle [GeV].
 * \param m3 Mass of first  outgoing particle [GeV].
 * \param m4 Mass of second outgoing particle [GeV].
 * \return array consisiting of {t_min, t_max}
 * Note that both t_min and t_max are negative,
 * with |t_min| < |t_max|, i.e. t_min > t_max.
 */
template <typename T>
std::array<T, 2> get_t_range(const T sqrts, const T m1, const T m2, const T m3,
                             const T m4) {
  const T p_i = pCM(sqrts, m1, m2);  // initial-state CM momentum
  const T p_f = pCM(sqrts, m3, m4);  // final-state CM momentum
  const T sqrt_t0 = (m1 * m1 - m2 * m2 - m3 * m3 + m4 * m4) / (2. * sqrts);
  const T t0 = sqrt_t0 * sqrt_t0;
  T t_min = t0 - (p_i - p_f) * (p_i - p_f);
  if (t_min > T(0.0)) {
    assert(t_min < T(1e-15));
    t_min = T(0.0);
  }
  const T t_max = t0 - (p_i + p_f) * (p_i + p_f);
  return {t_min, t_max};
}

/// Helper function for plab_from_s.
static inline void check_energy(double mandelstam_s, double m_sum) {
  if (mandelstam_s < m_sum * m_sum) {
    std::stringstream err;
    err << "plab_from_s: s too small: " << mandelstam_s << " < "
        << m_sum * m_sum;
    throw std::runtime_error(err.str());
  }
}

/// Helper function for plab_from_s.
static inline void check_radicand(double mandelstam_s, double radicand) {
  if (radicand < 0) {
    std::stringstream err;
    err << "plab_from_s: negative radicand: " << mandelstam_s;
    throw std::runtime_error(err.str());
  }
}

/** Convert mandelstam-s to p_lab in a fixed-target collision.
 * This assumes both particles have the given mass.
 */
inline double plab_from_s(double mandelstam_s, double mass) {
  const double radicand = mandelstam_s * (mandelstam_s - 4 * mass * mass);
#ifndef NDEBUG
  const double m_sum = 2 * mass;
  check_energy(mandelstam_s, m_sum);
  check_radicand(mandelstam_s, radicand);
#endif
  return std::sqrt(radicand) / (2 * mass);
}
/** Convert mandelstam-s to p_lab in a fixed-target collision.
 * This assumes both particles have the mass of a nucleon.
 */
inline double plab_from_s(double mandelstam_s) {
  return plab_from_s(mandelstam_s, nucleon_mass);
}
/** Convert mandelstam-s to p_lab in a fixed-target collision.
 * The mass of the projectile and the mass of the target have to be given.
 */
inline double plab_from_s(double mandelstam_s, double m_projectile,
                          double m_target) {
  const double m_sum = m_projectile + m_target;
  const double m_diff = m_projectile - m_target;
  const double radicand =
      (mandelstam_s - m_sum * m_sum) * (mandelstam_s - m_diff * m_diff);
/* This is equivalent to:
const double radicand
    = (mandelstam_s - m_a_sq - m_b_sq) * (mandelstam_s - m_a_sq - m_b_sq)
      - 4 * m_a_sq * m_b_sq;
*/
#ifndef NDEBUG
  check_energy(mandelstam_s, m_sum);
  check_radicand(mandelstam_s, radicand);
#endif
  return std::sqrt(radicand) / (2 * m_target);
}

/**
 * Convert E_kin to mandelstam-s for a fixed-target setup,
 * with a projectile of mass m_P and a kinetic energy e_kin
 * and a target of mass m_T at rest.
 */
inline double s_from_Ekin(double e_kin, double m_P, double m_T) {
  return m_P * m_P + m_T * m_T + 2 * m_T * (m_P + e_kin);
}

/**
 * Convert p_lab to mandelstam-s for a fixed-target setup,
 * with a projectile of mass m_P and momentum plab
 * and a target of mass m_T at rest.
 */
inline double s_from_plab(double plab, double m_P, double m_T) {
  return m_P * m_P + m_T * m_T + 2 * m_T * std::sqrt(m_P * m_P + plab * plab);
}

}  // namespace smash

#endif  // SRC_INCLUDE_KINEMATICS_H_
