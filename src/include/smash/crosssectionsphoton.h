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
/**
 * Calculation method for the cross sections. It has only one member at the
 * moment. In the future there will be more options.
 */
enum class ComputationMethod { Analytic };

template <ComputationMethod method>
class CrosssectionsPhoton {};

/**
 * Class to calculate the cross-section of a meson-meson to meson-photon
 * process. This template specialization uses the analytic formulas for the
 * cross-sections.
 */
template <>
class CrosssectionsPhoton<ComputationMethod::Analytic> {
 public:
  /** @name Total cross-section
   * The functions in this group calculate the analytical value of the total
   * cross-section for a photon process.
   *
   * \param[in] s Mandelstam-s [GeV^2]
   * \param[in] m_rho Mass of participating rho-meson [GeV]
   * \returns photon cross-section [mb]
   */
  ///@{
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
  ///@}

  /** @name Differential cross-section
   * The functions in this group calculate the analytical value of the
   * differential cross-section for a photon process.
   *
   * \param[in] s Mandelstam-s [GeV^2]
   * \param[in] t Mandelstam-t [GeV^2]
   * \param[in] m_rho Mass of participating rho-meson [GeV]
   * \returns photon cross-section [mb]
   */
  ///@{
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
  ///@}

  // Todo{JR}: Needed for tabulation.
  /// @cond
  static double s_min, s_max, t_min, t_max;
  /// @endcond

 private:
  /** @name Constants
   * The choice of these parameters, necessary to determine the photon cross
   * sections, follows from (\iref{Turbide:2006}). Here, different combinations
   * of parameters were proposed and investigated. We decided to use the
   * parameters of set II in their categorization. They are named identically,
   * except for "Const", which corresponds to "C" in \iref{Turbide:2006}.
   */
  ///@{
  constexpr static double Const = 0.059;
  constexpr static double g_POR = 22.6;
  constexpr static double ghat = 6.4483;
  constexpr static double eta1 = 2.3920;
  constexpr static double eta2 = 1.9430;
  constexpr static double delta = -0.6426;
  constexpr static double C4 = -0.14095;
  constexpr static double Gammaa1 = 0.4;
  ///@}
  /// Value of \f$ \pi \f$
  constexpr static double Pi = M_PI;
};

}  // namespace smash

#endif  // SRC_INCLUDE_CROSSSECTIONSPHOTON_H_
