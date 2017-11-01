/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include <cmath>
#include "kinematics.h"
#include "particletype.h"
#include "pdgcode.h"
#include "tabulationnd.h"

// calculation method for the cross sections
enum class ComputationMethod { Analytic, Lookup, Parametrized };

// usage would be
// PhotonCrossSection xs<Analytic>; xs.xs_pi_pi_rho(var1, var2...)
template <ComputationMethod method>
class PhotonCrossSection;

template <>
class PhotonCrossSection<ComputationMethod::Analytic> {
 public:
  static double xs_pi_pi_rho0(const double s);
  static double xs_pi_pi0_rho(const double s);
  static double xs_pi0_rho0_pi0(const double s);
  static double xs_pi_rho0_pi(const double s);
  static double xs_pi_rho_pi0(const double s);
  static double xs_pi0_rho_pi(const double s);

  static double xs_diff_pi0_rho0_pi0(const double s, const double t);
  static double xs_diff_pi_pi_rho0(const double s, const double t);
  static double xs_diff_pi_pi0_rho(const double s, const double t);
  static double xs_diff_pi_rho_pi0(const double s, const double t);
  static double xs_diff_pi0_rho_pi(const double s, const double t);
  static double xs_diff_pi_rho0_pi(const double s, const double t);

 private:

  // masses are hardcoded. since we use static tabulation objects we
  // can not make use of particleFinder. 
  constexpr static double to_mb = 0.3894;
  constexpr static double Const = 0.059;
  constexpr static double g_POR = 11.93;
  constexpr static double ma1 = 1.26;
  constexpr static double ghat = 6.4483;
  constexpr static double eta1 = 2.3920;
  constexpr static double eta2 = 1.9430;
  constexpr static double delta = -0.6426;
  constexpr static double C4 = -0.14095;
  constexpr static double Gammaa1 = 0.4;
  constexpr static double Pi = M_PI;
  constexpr static double m_omega = 0.783;

  constexpr static double m_pion_ = 0.139;
  constexpr static double m_rho_ = 0.775;
};

template <>
class PhotonCrossSection<ComputationMethod::Lookup> {
 private:
  static TabulationND<1> tab_pi_pi_rho0_;
  static TabulationND<1> tab_pi_pi0_rho_;
  static TabulationND<1> tab_pi0_rho0_pi0_;
  static TabulationND<1> tab_pi_rho0_pi_;
  static TabulationND<1> tab_pi_rho_pi0_;
  static TabulationND<1> tab_pi0_rho_pi_;

  static TabulationND<2> tab_pi0_rho0_pi0_diff_;
  static TabulationND<2> tab_pi_pi_rho0_diff_;
  static TabulationND<2> tab_pi_pi0_rho_diff_;
  static TabulationND<2> tab_pi_rho_pi0_diff_;
  static TabulationND<2> tab_pi0_rho_pi_diff_;
  static TabulationND<2> tab_pi_rho0_pi_diff_;

 public:
  double xs_pi_pi_rho0(const double s);
  double xs_pi_pi0_rho(const double s);
  double xs_pi0_rho0_pi0(const double s);
  double xs_pi_rho0_pi(const double s);
  double xs_pi_rho_pi0(const double s);
  double xs_pi0_rho_pi(const double s);

  double xs_diff_pi_pi_rho0(const double s, const double t);
  double xs_diff_pi_pi0_rho(const double s, const double t);
  double xs_diff_pi0_rho0_pi0(const double s, const double t);
  double xs_diff_pi_rho0_pi(const double s, const double t);
  double xs_diff_pi_rho_pi0(const double s, const double t);
  double xs_diff_pi0_rho_pi(const double s, const double t);
};

template <>
class PhotonCrossSection<ComputationMethod::Parametrized> {};
