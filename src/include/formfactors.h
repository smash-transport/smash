/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_FORMFACTORS_H_
#define SRC_INCLUDE_FORMFACTORS_H_

#include <string>

#include "constants.h"
#include "pdgcode_constants.h"

namespace Smash {


/**
 * Returns the squared Blatt-Weisskopf functions, which influence the mass
 * dependence of the decay widths. See e.g. \iref{Effenberger:1999wlg}, page 28.
 *
 * \param p_ab Momentum of outgoing particles A and B in center-of-mass frame.
 * \param L Angular momentum of outgoing particles A and B.
 *
 * This is used as a standard form factor for all hadronic decays. Note that all
 * the Blatt-Weisskopf functions approach one for large p_ab and behave like
 * p_ab**L for small p_ab. They are increasing monotonically with p_ab.
 */
inline float blatt_weisskopf_sqr(const float p_ab, const int L)
#ifdef NDEBUG
    noexcept
#endif
{
  const float R = 1. / hbarc;  /* interaction radius = 1 fm */
  const auto x = p_ab * R;
  const auto x2 = x * x;
  const auto x4 = x2 * x2;
  switch (L) {
    case 0:
      return 1.f;
    case 1:
      return x2 / (1.f + x2);
    case 2:
      return x4 / (9.f + 3.f * x2 + x4);
    case 3:
      return x4 * x2 / (225.f + 45.f * x2 + 6.f * x4 + x4 * x2);
    case 4:
      return x4 * x4 /
             (11025.f + 1575.f * x2 + 135.f * x4 + 10.f * x2 * x4 + x4 * x4);
      // See also input sanitization in load_decaymodes in decaymodes.cc.
#ifndef NDEBUG
    default:
      throw std::invalid_argument(
          std::string("Wrong angular momentum in BlattWeisskopf: ") +
          std::to_string(L));
#endif
  }
  return 0.f;
}


/**
 * An additional form factor for unstable final states as used in GiBUU,
 * according to M. Post, see eq. (174) in \iref{Buss:2011mx} or eq. (13) in
 * \iref{Post:2003hu}.
 *
 * \param m Actual mass of the decaying resonance [GeV].
 * \param M0 Pole mass of the decaying resonance [GeV].
 * \param srts0 Threshold of the reaction, i.e. minimum possible sqrt(s) [GeV].
 * \param L Lambda parameter of the form factor [GeV]. This is a cut-off
 * parameter that can be different for baryons and mesons.
 *
 * \return The squared value of the form factor (dimensionless).
 *
 * \note This form factor is equal to one at m=M0 and m=srts0. For decreasing
 * values of L, the form factor results in a stronger and stronger suppression
 * of the high-mass tail (m > M0) and a corresponding enhancement of the
 * low-mass tail (m < M0).
 */
inline float post_ff_sqr(float m, float M0, float srts0, float L) {
  const auto L4 = L*L*L*L;
  const auto M2 = M0*M0;
  const auto s0 = srts0*srts0;
  const auto sminus = (s0-M2)*0.5;
  const auto splus = m*m - (s0+M2)*0.5;
  const auto FF = (L4 + sminus*sminus) / (L4 + splus*splus);
  return FF*FF;
}


// electromagnetic transition form factors for the dilepton dalitz decays

/** Electromagnetic transition form factor for P → γ e⁺ e⁻, with a
 * pseudoscalar meson P = π⁰,η,η', as a function of the dilepton mass.
 * For the π⁰ see \iref{Landsberg:1986fd}. For the η the Lambda parameter is
 * fitted to NA60 data, see \iref{Arnaldi:2009aa}. */
inline float em_form_factor_ps(PdgCode pdg, float mass) {
  switch (pdg.code()) {
  case pdg::pi_z:
    return 1. + 5.5*mass*mass;
  case pdg::eta: {
    const float lambda_eta = 0.716;
    const float m_over_eta = mass / lambda_eta;
    return 1. / (1. - m_over_eta*m_over_eta);
  }
  default:  /* η' etc */
    return 1.;  // use QED approximation
  }
}

/** Squared electromagnetic transition form factor for V → π⁰ e⁺ e⁻, with a
 * vector meson V = ω,φ, as a function of the dilepton mass.
 * For the ω see \iref{Bratkovskaya:1996qe}. */
inline float em_form_factor_sqr_vec(PdgCode pdg, float mass) {
  switch (pdg.code()) {
  case pdg::omega: {
    constexpr float lambda = 0.65;
    constexpr float gamma = 0.075;
    constexpr float lambda_sqr = lambda * lambda;
    constexpr float gamma_sqr = gamma * gamma;
    const float tmp = lambda_sqr - mass*mass;
    const float denom = tmp*tmp + lambda_sqr*gamma_sqr;
    return lambda_sqr * lambda_sqr / denom;
  }
  default:  /* φ etc */
    return 1.;  // use QED approximation
  }
}

/** Electromagnetic transition form factor for Delta -> N e+ e-
 * as a function of the dilepton mass. Currently assumed to be constant,
 * normalized at the real-photon point. */
inline float form_factor_delta(float) {
  return 3.12;
}


}  // namespace Smash

#endif  // SRC_INCLUDE_FORMFACTORS_H_
