/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/decaytype.h"

#include <math.h>

#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/resonances.h"

namespace Smash {

/**
 * Returns the squared Blatt-Weisskopf functions,
 * which influence the mass dependence of the decay widths.
 * See e.g. Effenberger's thesis, page 28.
 *
 * \param p_ab Momentum of outgoing particles A and B in center-of-mass frame.
 * \param L Angular momentum of outgoing particles A and B.
 */
static float BlattWeisskopf(const float p_ab, const int L)
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
    case 2: {
      return x4 / (9.f + 3.f * x2 + x4);
    case 3:
      return x4 * x2 / (225.f + 45.f * x2 + 6.f * x4 + x4 * x2);
    }
    /* The following lines should be correct. But since nothing in SMASH uses
     * L > 3, this code is untested and dead. Therefore we only keep it as a
     * reference for later.
     * See also input sanitization in load_decaymodes in decaymodes.cc.
    case 4:
      return x4 * x4 /
             (11025.f + 1575.f * x2 + 135.f * x4 + 10.f * x2 * x4 + x4 * x4);
    */
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
 * according to M. Post. Reference: \iref{Buss:2011mx}, eq. (174).
 * The function returns the squared value of the form factor.
 * \param m Actual mass of the decaying resonance [GeV].
 * \param M0 Pole mass of the decaying resonance [GeV].
 * \param srts0 Threshold of the reaction, i.e. minimum possible sqrt(s) [GeV].
 * \param L Lambda parameter of the form factor [GeV]. This is a cut-off
 * parameter that can be different for baryons and mesons.
 */
static float Post_FF_sqr(float m, float M0, float srts0, float L) {
  const auto L4 = L*L*L*L;
  const auto m2 = m*m;
  const auto M2 = M0*M0;
  const auto s0 = srts0*srts0;
  const float FF = (L4 + (s0-M2)*(s0-M2)/4.) /
                   (L4 + (m2-(s0+M2)/2.) * (m2-(s0+M2)/2.));
  return FF*FF;
}


// TwoBodyDecay

TwoBodyDecay::TwoBodyDecay(ParticleTypePtrList part_types, int l)
                          : DecayType(part_types, l) {
  if (part_types.size() != 2) {
    throw std::runtime_error(
      "Wrong number of particles in TwoBodyDecay constructor: " +
      std::to_string(part_types.size()));
  }
}

int TwoBodyDecay::particle_number() const {
  return 2;
}

bool TwoBodyDecay::has_particles(const ParticleType &t_a,
                                 const ParticleType &t_b) const {
  return (*particle_types_[0] == t_a && *particle_types_[1] == t_b) ||
         (*particle_types_[0] == t_b && *particle_types_[1] == t_a);
}


// TwoBodyDecayStable

TwoBodyDecayStable::TwoBodyDecayStable(ParticleTypePtrList part_types, int l)
                                      : TwoBodyDecay(part_types, l) {
  if (!(part_types[0]->is_stable() && part_types[1]->is_stable())) {
    throw std::runtime_error(
      "Error: Unstable particle in TwoBodyDecayStable constructor: " +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayStable::rho(float m) const {
  // Determine momentum of outgoing particles in rest frame of resonance
  const float p_ab = pCM(m, particle_types_[0]->mass(),
                            particle_types_[1]->mass());
  // determine rho(m)
  return p_ab / m * BlattWeisskopf(p_ab, L_);
}

float TwoBodyDecayStable::width(float m0, float G0, float m) const {
  if (m <= particle_types_[0]->mass() + particle_types_[1]->mass()) {
    return 0;
  } else {
    return G0 * rho(m) / rho(m0);
  }
}

float TwoBodyDecayStable::in_width(float m0, float G0, float m,
                                   float, float) const {
  // in-width = out-width
  return width(m0, G0, m);
}


// TwoBodyDecaySemistable

static float integrand_rho_Manley(float mass, float srts, float stable_mass,
                                  ParticleTypePtr type, int L) {
  if (srts <= mass + stable_mass) {
    return 0.;
  }

  /* center-of-mass momentum of final state particles */
  const float p_f = pCM(srts, stable_mass, mass);

  return p_f/srts * BlattWeisskopf(p_f, L) * 2.*srts
         * spectral_function(mass, type->mass(), type->total_width(srts));
}

static ParticleTypePtrList arrange_particles(ParticleTypePtrList part_types) {
  // re-arrange the particle list such that the first particle is the stable one
  if (part_types[1]->is_stable()) {
    std::swap(part_types[0], part_types[1]);
  }
  /* verify that this is really a "semi-stable" decay,
   * i.e. the first particle is stable and the second unstable */
  if (!part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error(
      "Error in TwoBodyDecaySemistable constructor: " +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
  return part_types;
}

TwoBodyDecaySemistable::TwoBodyDecaySemistable(ParticleTypePtrList part_types,
                                               int l)
  : TwoBodyDecay(arrange_particles(part_types), l),
    Lambda_((particle_types_[1]->baryon_number() != 0) ? 2.0 : 1.6),
    tabulation_(particle_types_[0]->mass() + particle_types_[1]->minimum_mass(),
                1.f, 50,
                [&](float srts) {
                  Integrator integrate;
                  return integrate(particle_types_[1]->minimum_mass(),
                                    srts - particle_types_[0]->mass(),
                                    [&](float m) {
                                      return integrand_rho_Manley(m, srts,
                                                  particle_types_[0]->mass(),
                                                  particle_types_[1], L_);
                                    });
                })
{}

float TwoBodyDecaySemistable::rho(float m) const {
  return tabulation_.get_value_linear(m);
}

float TwoBodyDecaySemistable::width(float m0, float G0, float m) const {
  return G0 * rho(m) / rho(m0)
         * Post_FF_sqr(m, m0, particle_types_[0]->mass()
                              + particle_types_[1]->minimum_mass(), Lambda_);
}

float TwoBodyDecaySemistable::in_width(float m0, float G0, float m,
                                       float m1, float m2) const {
  const float p_f = pCM(m, m1, m2);

  return G0 * p_f * BlattWeisskopf(p_f, L_)
         * Post_FF_sqr(m, m0, particle_types_[0]->mass()
                              + particle_types_[1]->minimum_mass(), Lambda_)
         / (m * rho(m0));
}


// TwoBodyDecayUnstable

TwoBodyDecayUnstable::TwoBodyDecayUnstable(ParticleTypePtrList part_types,
                                           int l)
                                          : TwoBodyDecay(part_types, l) {
  if (part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error(
      "Error: Stable particle in TwoBodyDecayUnstable constructor: " +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayUnstable::width(float, float G0, float) const {
  return G0;  // use on-shell width
}

float TwoBodyDecayUnstable::in_width(float, float G0, float,
                                     float, float) const {
  return G0;  // use on-shell width
}

// TwoBodyDecayDilepton

TwoBodyDecayDilepton::TwoBodyDecayDilepton(ParticleTypePtrList part_types,
                                           int l)
                                      : TwoBodyDecayStable(part_types, l) {
  if (!is_dilepton(particle_types_[0]->pdgcode(),
                  particle_types_[1]->pdgcode())) {
    throw std::runtime_error(
      "Error: No dilepton in TwoBodyDecayDilepton constructor: " +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayDilepton::width(float m0, float G0, float m) const {
  if (m <= particle_types_[0]->mass() + particle_types_[1]->mass()) {
    return 0;
  } else {
    /// dilepton decays: use width from \iref{Li:1996mi}, equation (19)
    const float ml = particle_types_[0]->mass();  // lepton mass
    const float ml_to_m_sqr = (ml/m) * (ml/m);
    const float m0_to_m_cubed = (m0/m) * (m0/m) * (m0/m);
    return G0 * m0_to_m_cubed * std::sqrt(1.0f - 4.0f * ml_to_m_sqr) *
           (1.0f + 2.0f * ml_to_m_sqr);
  }
}

// ThreeBodyDecay

ThreeBodyDecay::ThreeBodyDecay(ParticleTypePtrList part_types, int l)
                               : DecayType(part_types, l) {
  if (part_types.size() != 3) {
    throw std::runtime_error(
      "Wrong number of particles in ThreeBodyDecay constructor: " +
      std::to_string(part_types.size()));
  }
}

int ThreeBodyDecay::particle_number() const {
  return 3;
}

bool ThreeBodyDecay::has_particles(const ParticleType &,
                                   const ParticleType &) const {
  return false;
}

float ThreeBodyDecay::width(float, float G0, float) const {
  return G0;  // use on-shell width
}

float ThreeBodyDecay::in_width(float, float G0, float, float, float) const {
  return G0;  // use on-shell width
}

// ThreeBodyDecayDilepton

ThreeBodyDecayDilepton::ThreeBodyDecayDilepton(ParticleTypePtrList part_types,
                                           int l)
                                      : ThreeBodyDecay(part_types, l) {
  if (!has_lepton_pair(particle_types_[0]->pdgcode(),
                       particle_types_[1]->pdgcode(),
                       particle_types_[2]->pdgcode())) {
   throw std::runtime_error(
     "Error: No dilepton in ThreeBodyDecayDilepton constructor: " +
     part_types[0]->pdgcode().string() + " " +
     part_types[1]->pdgcode().string() + " " +
     part_types[2]->pdgcode().string());
  }
}

// Little helper functions for the form factor of the dilepton dalitz decays

float form_factor_pi(float mass) {
  return 1.+5.5*mass*mass;
}

float form_factor_eta(float mass) {
  return 1./(1.-(mass*mass/0.676));
}

float form_factor_sqr_omega(float mass) {
  float lambda = 0.65;
  float gamma_w = 0.075;
  float n = pow(lambda*lambda - mass*mass, 2.) +
             lambda*gamma_w * lambda*gamma_w;
  return pow(lambda, 4.) / n;
}

float form_factor_sqr_delta(float) {
  return 1.0;  // no form factor implemented
}


float ThreeBodyDecayDilepton::diff_width(float m_par, float m_dil,
                                             float m_other, PdgCode pdg) {
  if (m_par < m_dil + m_other) {
    return 0;
  } else {
    float gamma = 0.0;
    // abbreviations
    float m_dil_sqr = m_dil * m_dil;
    float m_par_sqr = m_par * m_par;
    float m_par_cubed = m_par * m_par*m_par;
    float m_other_sqr = m_other*m_other;

    switch (pdg.get_decimal()) {
      case 111: /*pi0*/ {
        gamma = 7.8e-9;
        float ff = form_factor_pi(m_dil);
        return (alpha*4./(3.*M_PI)) * gamma/m_dil *
                                   pow(1.-m_dil/m_par*m_dil/m_par,3.) * ff*ff;
      }
      case 221: /*eta*/ {
        gamma = 46e-8;
        float ff = form_factor_eta(m_dil);
        return (4.*alpha/(3.*M_PI)) * gamma/m_dil *
                                   pow(1.-m_dil/m_par*m_dil/m_par,3.) * ff*ff;
      }
      case 223: /*omega*/ {
        gamma = 0.703e-3;
        float n1 = (m_par_sqr - m_other_sqr);
        float n2 = ((m_par_sqr -m_other_sqr)*(m_par_sqr -m_other_sqr));
        float rad = pow(1+ m_dil_sqr / n1, 2.)  -  4*m_par_sqr*m_dil_sqr / n2;
        return (2.*alpha/(3.*M_PI))  *  gamma/m_dil  *   pow(sqrt(rad), 3.) *
                                                   form_factor_sqr_omega(m_dil);
      }
      case 2214: case 2114: /* Delta+ and Delta0 */ {
        float t1 = alpha/16. *
                   (m_par+m_other)*(m_par+m_other)/(m_par_cubed*m_other_sqr) *
                   std::sqrt((m_par+m_other)*(m_par+m_other) - m_dil_sqr);
        float t2 =
              pow(std::sqrt((m_par-m_other)*(m_par-m_other) - m_dil_sqr), 3.0);
        float gamma_vi = t1 * t2 * form_factor_sqr_delta(m_dil);
        return 2.*alpha/(3.*M_PI) * gamma_vi/m_dil;
      }
      default:
        throw std::runtime_error("Error in ThreeBodyDecayDilepton");
    }
  }  // else
}

}  // namespace Smash
