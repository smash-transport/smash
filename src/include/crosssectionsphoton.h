/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_CROSSSECTIONSPHOTON_H_
#define SRC_INCLUDE_CROSSSECTIONSPHOTON_H_

#include "cxx14compat.h"
#include "kinematics.h"

namespace smash {
// calculation method for the cross sections
enum class ComputationMethod { Analytic };

// usage would be
// PhotonCrossSection xs<Analytic>; xs.xs_pi_pi_rho(var1, var2...)
template <ComputationMethod method>
class PhotonCrossSection {};

template <>
class PhotonCrossSection<ComputationMethod::Analytic> {
 public:
  static double xs_pi_pi_rho0(const double s, const double m_rho);
  static double xs_pi_pi0_rho(const double s, const double m_rho);
  static double xs_pi0_rho0_pi0(const double s, const double m_rho);
  static double xs_pi_rho0_pi(const double s, const double m_rho);

  static double xs_pi_rho_pi0(const double s, const double m_rho);
  static double xs_pi_rho_pi0_rho_mediated(const double s, const double m_rho);
  static double xs_pi_rho_pi0_omega_mediated(const double s,
                                             const double m_rho);

  static double xs_pi0_rho_pi(const double s, const double m_rho);
  static double xs_pi0_rho_pi_rho_mediated(const double s, const double m_rho);
  static double xs_pi0_rho_pi_omega_mediated(const double s,
                                             const double m_rho);

  static double xs_diff_pi_pi_rho0(const double s, const double t,
                                   const double m_rho);
  static double xs_diff_pi_pi0_rho(const double s, const double t,
                                   const double m_rho);
  static double xs_diff_pi0_rho0_pi0(const double s, const double t,
                                     const double m_rho);
  static double xs_diff_pi_rho0_pi(const double s, const double t,
                                   const double m_rho);

  static double xs_diff_pi_rho_pi0(const double s, const double t,
                                   const double m_rho);
  static double xs_diff_pi_rho_pi0_rho_mediated(const double s, const double t,
                                                const double m_rho);
  static double xs_diff_pi_rho_pi0_omega_mediated(const double s,
                                                  const double t,
                                                  const double m_rho);

  static double xs_diff_pi0_rho_pi(const double s, const double t,
                                   const double m_rho);
  static double xs_diff_pi0_rho_pi_rho_mediated(const double s, const double t,
                                                const double m_rho);
  static double xs_diff_pi0_rho_pi_omega_mediated(const double s,
                                                  const double t,
                                                  const double m_rho);

  static double s_min, s_max, t_min, t_max;

 private:
  /**
   * The choice of these parameters, necessary to determine the photon cross
   * sections, follows from (\iref{Turbide:2006}). Here, different combinations
   * of parameters were proposed and investigated. We decided to use the
   * parameters of set II in their categorization. They are named identically,
   * except for "Const", which corresponds to "C" in \iref{Turbide:2006}.
   */

  constexpr static double Const = 0.059;
  constexpr static double g_POR = 22.6;
  constexpr static double ghat = 6.4483;
  constexpr static double eta1 = 2.3920;
  constexpr static double eta2 = 1.9430;
  constexpr static double delta = -0.6426;
  constexpr static double C4 = -0.14095;
  constexpr static double Gammaa1 = 0.4;
  constexpr static double Pi = M_PI;
};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONSPHOTON_H_
