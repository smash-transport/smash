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
#include "include/formfactors.h"
#include "include/integrate.h"
#include "include/kinematics.h"

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

/**
 * Determine the cutoff parameter Λ for semistable decays,
 * given the types of the daughter particles.
 *
 * For the values used in GiBUU, see \iref{Buss:2011mx}, eq. (175).
 * For the original values used by M. Post, see table 1 in \iref{Post:2003hu}.
 *
 * We mostly stick to the GiBUU values, but use a different value for the ρπ
 * decay, in order to avoid secondary bumps in the ω spectral function and
 * achieve a better normalization. In contrast to Smash, GiBUU does not have
 * an ω → ρ π decay.
 */
static float get_Lambda(const ParticleTypePtr type_stable,
                        const ParticleTypePtr type_unstable) {
  if (type_unstable->baryon_number() != 0) {
    return 2.;  // unstable baryons
  } else if (type_unstable->pdgcode().is_rho() &&
             type_stable->pdgcode().is_pion()) {
    return 0.8;  // ρ+π
  } else {
    return 1.6;  // other unstable mesons
  }
}

TwoBodyDecaySemistable::TwoBodyDecaySemistable(ParticleTypePtrList part_types,
                                               int l)
  : TwoBodyDecay(arrange_particles(part_types), l),
    Lambda_(get_Lambda(particle_types_[0], particle_types_[1])),
    tabulation_(nullptr)
{}

float TwoBodyDecaySemistable::rho(float mass) const {
  if (tabulation_ == nullptr) {
    const_cast<TwoBodyDecaySemistable*>(this)->tabulation_
        = make_unique<Tabulation>(
                particle_types_[0]->mass() + particle_types_[1]->minimum_mass(),
                10*particle_types_[1]->width_at_pole(), 60,
                [&](float sqrts) {
                  Integrator integrate;
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
                                            Lambda_(2.0), tabulation_(nullptr) {
  if (part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error(
      "Error: Stable particle in TwoBodyDecayUnstable constructor: " +
      part_types[0]->pdgcode().string() + " " +
      part_types[1]->pdgcode().string());
  }
}

float TwoBodyDecayUnstable::rho(float mass) const {
  if (tabulation_ == nullptr) {
    const float m1_min = particle_types_[0]->minimum_mass();
    const float m2_min = particle_types_[1]->minimum_mass();
    const float sum_gamma = particle_types_[0]->width_at_pole()
                          + particle_types_[1]->width_at_pole();
    const_cast<TwoBodyDecayUnstable*>(this)->tabulation_
          = make_unique<Tabulation>(m1_min + m2_min, 10*sum_gamma, 60,
            [&](float sqrts) {
              Integrator2d integrate(1E4);
              return integrate(m1_min, sqrts - m2_min, m2_min, sqrts - m1_min,
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
                         : ThreeBodyDecay(part_types, l), tabulation_(nullptr) {
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


float ThreeBodyDecayDilepton::diff_width(float m_par, float m_dil,
                                             float m_other, PdgCode pdg) {
  if (m_par < m_dil + m_other) {
    return 0;
  } else {
    float gamma = 0.0;
    // abbreviations
    const float m_dil_sqr = m_dil * m_dil;
    const float m_par_sqr = m_par * m_par;
    const float m_par_cubed = m_par * m_par*m_par;
    const float m_other_sqr = m_other*m_other;

    switch (pdg.code()) {
      case 0x111: /* pi0 */ {
        /// see \iref{Landsberg:1986fd}, equation (3.8)
        gamma = 7.6e-9;
        const float ff = form_factor_pi(m_dil);
        return (4.*alpha/(3.*M_PI)) * gamma/m_dil *
                                    pow(1.-m_dil/m_par*m_dil/m_par, 3.) * ff*ff;
      }
      case 0x221: /* eta */ {
        /// see \iref{Landsberg:1986fd}, equation (3.8)
        gamma = 52e-8;
        const float ff = form_factor_eta(m_dil);
        return (4.*alpha/(3.*M_PI)) * gamma/m_dil *
                                    pow(1.-m_dil/m_par*m_dil/m_par, 3.) * ff*ff;
      }
      case 0x223: /* omega */ {
        /// see \iref{Landsberg:1986fd}, equation (3.4)
        gamma = 0.703e-3;
        const float n1 = (m_par_sqr - m_other_sqr);
        const float n2 = ((m_par_sqr -m_other_sqr)*(m_par_sqr -m_other_sqr));
        const float rad = pow(1+ m_dil_sqr / n1, 2.)  -
                          4*m_par_sqr*m_dil_sqr / n2;
        return (2.*alpha/(3.*M_PI))  *  gamma/m_dil  *   pow(sqrt(rad), 3.) *
                                                   form_factor_sqr_omega(m_dil);
      }
      case 0x2214: case 0x2114: /* Delta+ and Delta0 */ {
        /// see \iref{Krivoruchenko:2001hs}
        const float rad1 = (m_par+m_other)*(m_par+m_other) - m_dil_sqr;
        const float rad2 = (m_par-m_other)*(m_par-m_other) - m_dil_sqr;
        const float t1 = alpha/16. *
                   (m_par+m_other)*(m_par+m_other)/(m_par_cubed*m_other_sqr) *
                   std::sqrt(rad1);
        const float t2 =
              pow(std::sqrt(rad2), 3.0);
        const float ff = form_factor_delta(m_dil);
        const float gamma_vi = t1 * t2 * ff*ff;
        return 2.*alpha/(3.*M_PI) * gamma_vi/m_dil;
      }
      default:
        throw std::runtime_error("Error in ThreeBodyDecayDilepton");
    }
  }  // else
}

float ThreeBodyDecayDilepton::width(float, float G0, float m) const {
  PdgCode pdg_par;
  int non_lepton_position = -1;

  for (int i = 0; i < 3; ++i) {
    if (particle_types_[i]->pdgcode() == 0x111) {
      pdg_par = 0x223;  // only omega decays into a lepton pair and a pi0
      non_lepton_position = i;
      break;
    }
    if (particle_types_[i]->pdgcode() == 0x2212) {
      pdg_par = 0x2214;  // only Delta+ decays into a lepton pair and a proton
      non_lepton_position = i;
      break;
    }
    if (particle_types_[i]->pdgcode() == 0x2112) {
      pdg_par = 0x2114;  // only Delta0 decays into a lepton pair and a neutron
      non_lepton_position = i;
      break;
    }
    if (particle_types_[i]->pdgcode() == 0x22) {
      // Only eta and pi0 decay into lepton pair and a photon. We assume here
      // that their width is on-shell.
      return G0;
    }
  }

  if (pdg_par == 0x0 || non_lepton_position == -1) {
    throw std::runtime_error("Error unsupported Dalitz Dilepton Decay");
  }

  // lepton mass
  const float m_l = particle_types_[(non_lepton_position+1)%3]->mass();
  // mass of non-leptonic particle in final state
  const float m_other = particle_types_[non_lepton_position]->mass();

  // integrate differential width to obtain partial width
  if (tabulation_ == nullptr) {
    const_cast<ThreeBodyDecayDilepton*>(this)->tabulation_
          = make_unique<Tabulation>(m_other+2*m_l, 10*G0, 60,
              [&](float m_parent) {
                const float bottom = 2*m_l;
                const float top = m_parent-m_other;
                if (top < bottom) {  // numerical problems at lower bound
                  return 0.0;
                }
                Integrator integrate;
                return integrate(bottom, top,
                               [&](float sqrts) {
                                  return diff_width(m_parent, sqrts,
                                                    m_other, pdg_par);
                                }).value();
                });
  }
  return tabulation_->get_value_linear(m);
}

}  // namespace Smash
