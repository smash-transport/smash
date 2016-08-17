/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "include/decaytype.h"

#include <algorithm>
#include <math.h>

#include "include/constants.h"
#include "include/cxx14compat.h"
#include "include/formfactors.h"
#include "include/integrate.h"
#include "include/kinematics.h"
#include "include/pdgcode_constants.h"

namespace Smash {

// auxiliary functions

static float integrand_rho_Manley_1res(float sqrts, float mass,
                              float stable_mass, ParticleTypePtr type, int L) {
  if (sqrts <= mass + stable_mass) {
    return 0.;
  }

  /* center-of-mass momentum of final state particles */
  const float p_f = pCM(sqrts, stable_mass, mass);

  return p_f/sqrts * blatt_weisskopf_sqr(p_f, L)
         * type->spectral_function(mass);
}

static float integrand_rho_Manley_2res(float sqrts, float m1, float m2,
                                ParticleTypePtr t1, ParticleTypePtr t2, int L) {
  if (sqrts <= m1 + m2) {
    return 0.;
  }

  /* center-of-mass momentum of final state particles */
  const float p_f = pCM(sqrts, m1, m2);

  return p_f/sqrts * blatt_weisskopf_sqr(p_f, L)
         * t1->spectral_function(m1) * t2->spectral_function(m2);
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

unsigned int TwoBodyDecay::particle_number() const {
  return 2;
}

bool TwoBodyDecay::has_particles(ParticleTypePtrList list) const {
  if (list.size() != particle_number()) {
    return false;
  }
  return (particle_types_[0] == list[0] && particle_types_[1] == list[1]) ||
         (particle_types_[0] == list[1] && particle_types_[1] == list[0]);
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
  return p_ab / m * blatt_weisskopf_sqr(p_ab, L_);
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

/// re-arrange the particle list such that the first particle is the stable one
static ParticleTypePtrList arrange_particles(ParticleTypePtrList part_types) {
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
                                Lambda_(get_Lambda()),
                                tabulation_(nullptr)
{}

float TwoBodyDecaySemistable::get_Lambda() {
  // "semi-stable" decays (first daughter is stable and second one unstable)
  if (particle_types_[1]->baryon_number() != 0) {
    return 2.;  // unstable baryons
  } else if (particle_types_[1]->pdgcode().is_rho() &&
             particle_types_[0]->pdgcode().is_pion()) {
    return 0.8;  // ρ+π
  } else {
    return 1.6;  // other unstable mesons
  }
}

// number of tabulation points
constexpr int num_tab_pts = 200;
static thread_local Integrator integrate;

float TwoBodyDecaySemistable::rho(float mass) const {
  if (tabulation_ == nullptr) {
    /* TODO(weil): Move this lazy init to a global initialization function,
      * in order to avoid race conditions in multi-threading. */
    tabulation_ = make_unique<Tabulation>(
                particle_types_[0]->mass() + particle_types_[1]->minimum_mass(),
                10*particle_types_[1]->width_at_pole(), num_tab_pts,
                [&](float sqrts) {
                  return integrate(particle_types_[1]->minimum_mass(),
                                    sqrts - particle_types_[0]->mass(),
                                    [&](float m) {
                                      return integrand_rho_Manley_1res(sqrts, m,
                                                  particle_types_[0]->mass(),
                                                  particle_types_[1], L_);
                                    });
                });
  }
  return tabulation_->get_value_linear(mass);
}

float TwoBodyDecaySemistable::width(float m0, float G0, float m) const {
  return G0 * rho(m) / rho(m0)
         * post_ff_sqr(m, m0, particle_types_[0]->mass()
                              + particle_types_[1]->minimum_mass(), Lambda_);
}

float TwoBodyDecaySemistable::in_width(float m0, float G0, float m,
                                       float m1, float m2) const {
  const float p_f = pCM(m, m1, m2);

  return G0 * p_f * blatt_weisskopf_sqr(p_f, L_)
         * post_ff_sqr(m, m0, particle_types_[0]->mass()
                              + particle_types_[1]->minimum_mass(), Lambda_)
         / (m * rho(m0));
}


// TwoBodyDecayUnstable

TwoBodyDecayUnstable::TwoBodyDecayUnstable(ParticleTypePtrList part_types,
                                           int l)
                                          : TwoBodyDecay(part_types, l),
                                            Lambda_(get_Lambda()),
                                            tabulation_(nullptr) {
  if (part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error(
      "Error: Stable particle in TwoBodyDecayUnstable constructor: " +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayUnstable::get_Lambda() {
  // for now: use the same value for all unstable decays (fixed on f₂ → ρ ρ)
  return 0.6;
}

static thread_local Integrator2d integrate2d(1E4);

float TwoBodyDecayUnstable::rho(float mass) const {
  if (tabulation_ == nullptr) {
    /* TODO(weil): Move this lazy init to a global initialization function,
      * in order to avoid race conditions in multi-threading. */
    const float m1_min = particle_types_[0]->minimum_mass();
    const float m2_min = particle_types_[1]->minimum_mass();
    const float sum_gamma = particle_types_[0]->width_at_pole()
                          + particle_types_[1]->width_at_pole();
    tabulation_
          = make_unique<Tabulation>(m1_min + m2_min, 10*sum_gamma, num_tab_pts,
            [&](float sqrts) {
              return integrate2d(m1_min, sqrts - m2_min, m2_min, sqrts - m1_min,
                                [&](float m1, float m2) {
                                  return integrand_rho_Manley_2res(sqrts,
                                                    m1, m2, particle_types_[0],
                                                    particle_types_[1], L_);
                                });
            });
  }
  return tabulation_->get_value_linear(mass);
}

float TwoBodyDecayUnstable::width(float m0, float G0, float m) const {
  return G0 * rho(m) / rho(m0)
            * post_ff_sqr(m, m0, particle_types_[0]->minimum_mass()
                               + particle_types_[1]->minimum_mass(), Lambda_);
}

float TwoBodyDecayUnstable::in_width(float m0, float G0, float m,
                                     float m1, float m2) const {
  const float p_f = pCM(m, m1, m2);

  return G0 * p_f * blatt_weisskopf_sqr(p_f, L_)
         * post_ff_sqr(m, m0, particle_types_[0]->minimum_mass()
                            + particle_types_[1]->minimum_mass(), Lambda_)
         / (m * rho(m0));
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

/// sort the particle list
static ParticleTypePtrList sort_particles(ParticleTypePtrList part_types) {
  std::sort(part_types.begin(), part_types.end());
  return part_types;
}

ThreeBodyDecay::ThreeBodyDecay(ParticleTypePtrList part_types, int l)
                               : DecayType(sort_particles(part_types), l) {
  if (part_types.size() != 3) {
    throw std::runtime_error(
      "Wrong number of particles in ThreeBodyDecay constructor: " +
      std::to_string(part_types.size()));
  }
}

unsigned int ThreeBodyDecay::particle_number() const {
  return 3;
}

bool ThreeBodyDecay::has_particles(ParticleTypePtrList list) const {
  if (list.size() != particle_number()) {
    return false;
  }
  // compare sorted vectors (particle_types_ is already sorted)
  std::sort(list.begin(), list.end());
  return particle_types_[0] == list[0] && particle_types_[1] == list[1] &&
         particle_types_[2] == list[2];
}

float ThreeBodyDecay::width(float, float G0, float) const {
  return G0;  // use on-shell width
}

float ThreeBodyDecay::in_width(float, float G0, float, float, float) const {
  return G0;  // use on-shell width
}


// ThreeBodyDecayDilepton

ThreeBodyDecayDilepton::ThreeBodyDecayDilepton(ParticleTypePtr mother,
                                          ParticleTypePtrList part_types, int l)
                         : ThreeBodyDecay(part_types, l), tabulation_(nullptr),
                           mother_(mother) {
  if (!has_lepton_pair(particle_types_[0]->pdgcode(),
                       particle_types_[1]->pdgcode(),
                       particle_types_[2]->pdgcode())) {
    throw std::runtime_error(
     "Error: No dilepton in ThreeBodyDecayDilepton constructor: " +
     part_types[0]->pdgcode().string() + " " +
     part_types[1]->pdgcode().string() + " " +
     part_types[2]->pdgcode().string());
  }

  int non_lepton_position = -1;
  for (int i = 0; i < 3; ++i) {
    if (!particle_types_[i]->is_lepton()) {
      non_lepton_position = i;
      break;
    }
  }

  if (mother->pdgcode() == pdg::invalid || non_lepton_position == -1) {
    throw std::runtime_error("Error: Unsupported dilepton Dalitz decay!");
  }
}

bool ThreeBodyDecayDilepton::has_mother(ParticleTypePtr mother) const {
  return mother == mother_;
}

float ThreeBodyDecayDilepton::diff_width(float m_par, float m_dil,
                                         float m_other, ParticleTypePtr t) {
  // check threshold
  if (m_par < m_dil + m_other) {
    return 0;
  }

  // abbreviations
  const float m_dil_sqr = m_dil * m_dil;
  const float m_par_sqr = m_par * m_par;
  const float m_par_cubed = m_par * m_par*m_par;
  const float m_other_sqr = m_other*m_other;

  PdgCode pdg = t->pdgcode();
  if (pdg.is_meson()) {
    const ParticleType &photon = ParticleType::find(pdg::photon);
    const ParticleType &pi0 = ParticleType::find(pdg::pi_z);
    switch (pdg.spin()) {
    case 0:  /* pseudoscalars: π⁰, η, η' */ {
      // width for decay into 2γ
      const float gamma_2g = t->get_partial_width(m_par, photon, photon);
      float ff = em_form_factor_ps(pdg, m_dil);  // form factor
      /// see \iref{Landsberg:1986fd}, equation (3.8)
      return (4.*alpha/(3.*M_PI)) * gamma_2g/m_dil
                                  * pow(1.-m_dil/m_par*m_dil/m_par, 3.) * ff*ff;
    }
    case 2: /* vectors: ω, φ */ {
      // width for decay into π⁰γ
      const float gamma_pig = t->get_partial_width(m_par, pi0, photon);
      float ff_sqr = em_form_factor_sqr_vec(pdg, m_dil);  // form factor squared
      /// see \iref{Landsberg:1986fd}, equation (3.4)
      const float n1 = m_par_sqr - m_other_sqr;
      const float rad = pow(1. + m_dil_sqr/n1, 2)
                        - 4.*m_par_sqr*m_dil_sqr/(n1*n1);
      if (rad < 0.) {
        assert(rad > -1E-5);
        return 0.;
      } else {
        return (2.*alpha/(3.*M_PI)) * gamma_pig/m_dil
               * pow(rad, 3./2.) * ff_sqr;
      }
    }
    default:
      throw std::runtime_error("Bad meson in ThreeBodyDecayDilepton: "
                               + pdg.string());
    }
  } else if (pdg.is_baryon()) {
    switch (pdg.code()) {
      case pdg::Delta_p: case -pdg::Delta_p:
      case pdg::Delta_z: case -pdg::Delta_z:  /* Δ⁺, Δ⁰ (and antiparticles) */ {
      /// see \iref{Krivoruchenko:2001hs}
      const float rad1 = (m_par+m_other)*(m_par+m_other) - m_dil_sqr;
      const float rad2 = (m_par-m_other)*(m_par-m_other) - m_dil_sqr;
      const float t1 = alpha/16. *
                  (m_par+m_other)*(m_par+m_other)/(m_par_cubed*m_other_sqr) *
                  std::sqrt(rad1);
      const float t2 = pow(std::sqrt(rad2), 3.0);
      const float ff = form_factor_delta(m_dil);
      const float gamma_vi = t1 * t2 * ff*ff;
      return 2.*alpha/(3.*M_PI) * gamma_vi/m_dil;
    }
    default:
      throw std::runtime_error("Bad baryon in ThreeBodyDecayDilepton: "
                               + pdg.string());
    }
  } else {
    throw std::runtime_error("Non-hadron in ThreeBodyDecayDilepton: "
                             + pdg.string());
  }
}


float ThreeBodyDecayDilepton::width(float, float G0, float m) const {
  if (mother_->is_stable()) {
    return G0;
  }

  if (!tabulation_) {
    int non_lepton_position = -1;
    for (int i = 0; i < 3; ++i) {
      if (!particle_types_[i]->is_lepton()) {
        non_lepton_position = i;
        break;
      }
    }
    // lepton mass
    const float m_l = particle_types_[(non_lepton_position+1)%3]->mass();
    // mass of non-leptonic particle in final state
    const float m_other = particle_types_[non_lepton_position]->mass();

    // integrate differential width to obtain partial width
    float M0 = mother_->mass();
    float G0tot = mother_->width_at_pole();
    tabulation_
          = make_unique<Tabulation>(m_other+2*m_l, M0 + 10*G0tot, num_tab_pts,
              [&](float m_parent) {
                const float bottom = 2*m_l;
                const float top = m_parent-m_other;
                if (top < bottom) {  // numerical problems at lower bound
                  return 0.;
                }
                return integrate(bottom, top,
                                 [&](float m_dil) {
                                    return diff_width(m_parent, m_dil, m_other,
                                                      mother_);
                                 }).value();
                });
  }

  return tabulation_->get_value_linear(m);
}

}  // namespace Smash
