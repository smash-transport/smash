/*
 *
 *    Copyright (c) 2016
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */


#ifndef SRC_INCLUDE_PHOTONCROSSSECTION_H
#define SRC_INCLUDE_PHOTONCROSSSECTION_H


#include <cmath>
#include <memory>

#include <iostream>

#include "logging.h"

#include "kinematics.h"
#include "particletype.h"
#include "pdgcode.h"
#include "tabulationnd.h"
#include "cxx14compat.h"

namespace smash {
// calculation method for the cross sections
enum class ComputationMethod { Analytic, Lookup, Parametrized };

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
  static double xs_pi_rho_pi0_omega_mediated(const double s, const double m_rho);
  
  static double xs_pi0_rho_pi(const double s, const double m_rho);
  static double xs_pi0_rho_pi_rho_mediated(const double s, const double m_rho);
  static double xs_pi0_rho_pi_omega_mediated(const double s, const double m_rho);


  static double xs_diff_pi_pi_rho0(const double s, const double t, const double m_rho);
  static double xs_diff_pi_pi0_rho(const double s, const double t, const double m_rho);
  static double xs_diff_pi0_rho0_pi0(const double s, const double t, const double m_rho);
  static double xs_diff_pi_rho0_pi(const double s, const double t, const double m_rho);
  
  static double xs_diff_pi_rho_pi0(const double s, const double t, const double m_rho);
  static double xs_diff_pi_rho_pi0_rho_mediated(const double s, const double t, const double m_rho);
  static double xs_diff_pi_rho_pi0_omega_mediated(const double s, const double t, const double m_rho);
  
  static double xs_diff_pi0_rho_pi(const double s, const double t, const double m_rho);
  static double xs_diff_pi0_rho_pi_rho_mediated(const double s, const double t, const double m_rho);
  static double xs_diff_pi0_rho_pi_omega_mediated(const double s, const double t, const double m_rho);

  static double s_min, s_max, t_min, t_max;

 private:

  // masses are hardcoded. since we use static tabulation objects we
  // can not make use of particleFinder.
  constexpr static double to_mb = 0.3894;
  constexpr static double Const = 0.059;
  //constexpr static double g_POR = 11.93;
  constexpr static double g_POR = 22.6;
  constexpr static double ma1 = 1.26;
  constexpr static double ghat = 6.4483;
  constexpr static double eta1 = 2.3920;
  constexpr static double eta2 = 1.9430;
  constexpr static double delta = -0.6426;
  constexpr static double C4 = -0.14095;
  constexpr static double Gammaa1 = 0.4;
  constexpr static double Pi = M_PI;
  constexpr static double m_omega_ = 0.783;
  constexpr static double m_pion_ = 0.139;

  // Heaviside step function. This function is not part of the cmath library,
  // we need to define it ourselves. It is necessary for the implementation
  // of a non-stable rho meson in the pi0 + rho0 -> pi0 + gamma channel. For
  // this specific scattering process, there are s and t channels with different
  // theresholds. For s-channels, this is the omega meson mass while it is the
  // sum of pi and rho mass for the t-channel. This is problematic for low rho
  // masses, for example m_rho = 500 MeV. In this case, the s-channel is
  // kinematically not accessible and need thus be excluded from the cross
  // sections. Thus, we need the Heaviside function.

  static double HeavisideTheta(double x) {

  	if (x >= 0.0){
  		return 1.0;
  	}
  	else{
  		return 0.0;
  	}
  }
};

template <>
class PhotonCrossSection<ComputationMethod::Lookup> {
 private:
  static std::unique_ptr<TabulationND<2>> tab_pi_pi_rho0_;
  static std::unique_ptr<TabulationND<2>> tab_pi_pi0_rho_;
  static std::unique_ptr<TabulationND<2>> tab_pi0_rho0_pi0_;
  static std::unique_ptr<TabulationND<2>> tab_pi_rho0_pi_;

  static std::unique_ptr<TabulationND<2>> tab_pi_rho_pi0_;
  static std::unique_ptr<TabulationND<2>> tab_pi_rho_pi0_rho_;
  static std::unique_ptr<TabulationND<2>> tab_pi_rho_pi0_omega_;
  
  
  static std::unique_ptr<TabulationND<2>> tab_pi0_rho_pi_;
  static std::unique_ptr<TabulationND<2>> tab_pi0_rho_pi_rho_;
  static std::unique_ptr<TabulationND<2>> tab_pi0_rho_pi_omega_;

  static std::unique_ptr<TabulationND<3>> tab_pi0_rho0_pi0_diff_;
  static std::unique_ptr<TabulationND<3>> tab_pi_pi_rho0_diff_;
  static std::unique_ptr<TabulationND<3>> tab_pi_pi0_rho_diff_;

  static std::unique_ptr<TabulationND<3>> tab_pi_rho_pi0_diff_;
  static std::unique_ptr<TabulationND<3>> tab_pi_rho_pi0_rho_diff_;
  static std::unique_ptr<TabulationND<3>> tab_pi_rho_pi0_omega_diff_;
  

  static std::unique_ptr<TabulationND<3>> tab_pi0_rho_pi_diff_;
  static std::unique_ptr<TabulationND<3>> tab_pi0_rho_pi_rho_diff_;
  static std::unique_ptr<TabulationND<3>> tab_pi0_rho_pi_omega_diff_;
  
  static std::unique_ptr<TabulationND<3>> tab_pi_rho0_pi_diff_;



 public:

  const double s0_diff = 0.01, s1_diff = 20.0, t0_diff = -18.0;
  const double t1_diff = 5.0, ds_diff = 0.01, dt_diff = 0.01;
  const double s0_tot = 0.01, s1_tot = 20.0, ds_tot = 0.01;
  const double m_rho_0 = 0.775, m_rho_1 = 0.777, dm = 0.01;

  double xs_pi_pi_rho0(const double s, const double m_rho);
  double xs_pi_pi0_rho(const double s, const double m_rho);
  double xs_pi0_rho0_pi0(const double s, const double m_rho);
  double xs_pi_rho0_pi(const double s, const double m_rho);
  
  double xs_pi_rho_pi0(const double s, const double m_rho);
  double xs_pi_rho_pi0_rho_mediated(const double s, const double m_rho);
  double xs_pi_rho_pi0_omega_mediated(const double s, const double m_rho);
  
  double xs_pi0_rho_pi(const double s, const double m_rho);
  double xs_pi0_rho_pi_rho_mediated(const double s, const double m_rho);
  double xs_pi0_rho_pi_omega_mediated(const double s, const double m_rho);


  double xs_diff_pi_pi_rho0(const double s, const double t, const double m_rho);
  double xs_diff_pi_pi0_rho(const double s, const double t, const double m_rho);
  double xs_diff_pi0_rho0_pi0(const double s, const double t, const double m_rho);
  double xs_diff_pi_rho0_pi(const double s, const double t, const double m_rho);
  
  double xs_diff_pi_rho_pi0(const double s, const double t, const double m_rho);
  double xs_diff_pi_rho_pi0_rho_mediated(const double s, const double t, const double m_rho);
  double xs_diff_pi_rho_pi0_omega_mediated(const double s, const double t, const double m_rho);
  
  double xs_diff_pi0_rho_pi(const double s, const double t, const double m_rho);
  double xs_diff_pi0_rho_pi_rho_mediated(const double s, const double t, const double m_rho);
  double xs_diff_pi0_rho_pi_omega_mediated(const double s, const double t, const double m_rho);
  
};

template <>
class PhotonCrossSection<ComputationMethod::Parametrized> {};

} // namespace smash

#endif
