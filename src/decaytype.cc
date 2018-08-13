/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#include "smash/decaytype.h"

#include <algorithm>
#include <cmath>

#include "smash/constants.h"
#include "smash/cxx14compat.h"
#include "smash/formfactors.h"
#include "smash/integrate.h"
#include "smash/kinematics.h"
#include "smash/pdgcode_constants.h"
#include "smash/pow.h"

namespace smash {

// auxiliary functions

static double integrand_rho_Manley_1res(double sqrts, double mass,
                                        double stable_mass,
                                        ParticleTypePtr type, int L) {
  if (sqrts <= mass + stable_mass) {
    return 0.;
  }

  /* center-of-mass momentum of final state particles */
  const double p_f = pCM(sqrts, stable_mass, mass);

  return p_f / sqrts * blatt_weisskopf_sqr(p_f, L) *
         type->spectral_function(mass);
}

static double integrand_rho_Manley_2res(double sqrts, double m1, double m2,
                                        ParticleTypePtr t1, ParticleTypePtr t2,
                                        int L) {
  if (sqrts <= m1 + m2) {
    return 0.;
  }

  /* center-of-mass momentum of final state particles */
  const double p_f = pCM(sqrts, m1, m2);

  return p_f / sqrts * blatt_weisskopf_sqr(p_f, L) * t1->spectral_function(m1) *
         t2->spectral_function(m2);
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

unsigned int TwoBodyDecay::particle_number() const { return 2; }

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

double TwoBodyDecayStable::rho(double m) const {
  // Determine momentum of outgoing particles in rest frame of the resonance
  const double p_ab =
      pCM(m, particle_types_[0]->mass(), particle_types_[1]->mass());
  // determine rho(m)
  return p_ab / m * blatt_weisskopf_sqr(p_ab, L_);
}

double TwoBodyDecayStable::width(double m0, double G0, double m) const {
  assert(rho(m0) != 0);
  return (m <= threshold()) ? 0. : G0 * rho(m) / rho(m0);
}

double TwoBodyDecayStable::in_width(double m0, double G0, double m, double,
                                    double) const {
  // in-width = out-width
  return width(m0, G0, m);
}

// TwoBodyDecaySemistable

/**
 * Rearrange the particle list such that the first particle is the stable one.
 *
 * \param[inout] part_types Particle list to be rearranged.
 * \return Reference to rearranged particle list.
 */
static ParticleTypePtrList& arrange_particles(ParticleTypePtrList& part_types) {
  if (part_types[1]->is_stable()) {
    std::swap(part_types[0], part_types[1]);
  }
  /* verify that this is really a "semi-stable" decay,
   * i.e. the first particle is stable and the second unstable */
  if (!part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error("Error in TwoBodyDecaySemistable constructor: " +
                             part_types[0]->pdgcode().string() + " " +
                             part_types[1]->pdgcode().string());
  }
  return part_types;
}

TwoBodyDecaySemistable::TwoBodyDecaySemistable(ParticleTypePtrList part_types,
                                               int l)
    : TwoBodyDecay(arrange_particles(part_types), l),
      Lambda_(get_Lambda()),
      tabulation_(nullptr) {}

double TwoBodyDecaySemistable::get_Lambda() {
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
static /*thread_local (see #3075)*/ Integrator integrate;

double TwoBodyDecaySemistable::rho(double mass) const {
  if (tabulation_ == nullptr) {
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    const ParticleTypePtr res = particle_types_[1];
    const double tabulation_interval = std::max(2., 10. * res->width_at_pole());
    const double m_stable = particle_types_[0]->mass();
    const double mres_min = res->min_mass_kinematic();

    tabulation_ = make_unique<Tabulation>(
        threshold(), tabulation_interval, num_tab_pts, [&](double sqrts) {
          const double mres_max = sqrts - m_stable;
          return integrate(mres_min, mres_max, [&](double m) {
            return integrand_rho_Manley_1res(sqrts, m, m_stable, res, L_);
          });
        });
  }
  return tabulation_->get_value_linear(mass);
}

double TwoBodyDecaySemistable::width(double m0, double G0, double m) const {
  assert(rho(m0) != 0);
  return G0 * rho(m) / rho(m0) * post_ff_sqr(m, m0, threshold(), Lambda_);
}

double TwoBodyDecaySemistable::in_width(double m0, double G0, double m,
                                        double m1, double m2) const {
  assert(rho(m0) != 0);
  const double p_f = pCM(m, m1, m2);

  return G0 * p_f * blatt_weisskopf_sqr(p_f, L_) *
         post_ff_sqr(m, m0, threshold(), Lambda_) / (m * rho(m0));
}

// TwoBodyDecayUnstable

TwoBodyDecayUnstable::TwoBodyDecayUnstable(ParticleTypePtrList part_types,
                                           int l)
    : TwoBodyDecay(part_types, l), Lambda_(get_Lambda()), tabulation_(nullptr) {
  if (part_types[0]->is_stable() || part_types[1]->is_stable()) {
    throw std::runtime_error(
        "Error: Stable particle in TwoBodyDecayUnstable constructor: " +
        part_types[0]->pdgcode().string() + " " +
        part_types[1]->pdgcode().string());
  }
}

double TwoBodyDecayUnstable::get_Lambda() {
  // for now: use the same value for all unstable decays (fixed on f₂ → ρ ρ)
  return 0.6;
}

static /*thread_local*/ Integrator2dCuhre integrate2d(1E7);

double TwoBodyDecayUnstable::rho(double mass) const {
  if (tabulation_ == nullptr) {
    /* TODO(weil): Move this lazy init to a global initialization function,
     * in order to avoid race conditions in multi-threading. */
    const ParticleTypePtr r1 = particle_types_[0];
    const ParticleTypePtr r2 = particle_types_[1];
    const double m1_min = r1->min_mass_kinematic();
    const double m2_min = r2->min_mass_kinematic();
    const double sum_gamma = r1->width_at_pole() + r2->width_at_pole();
    const double tab_interval = std::max(2., 10. * sum_gamma);

    tabulation_ = make_unique<Tabulation>(
        m1_min + m2_min, tab_interval, num_tab_pts, [&](double sqrts) {
          const double m1_max = sqrts - m2_min;
          const double m2_max = sqrts - m1_min;

          const double result = integrate2d(m1_min, m1_max, m2_min, m2_max,
                                            [&](double m1, double m2) {
                                              return integrand_rho_Manley_2res(
                                                  sqrts, m1, m2, r1, r2, L_);
                                            })
                                    .value();
          return result;
        });
  }
  return tabulation_->get_value_linear(mass);
}

double TwoBodyDecayUnstable::width(double m0, double G0, double m) const {
  return G0 * rho(m) / rho(m0) * post_ff_sqr(m, m0, threshold(), Lambda_);
}

double TwoBodyDecayUnstable::in_width(double m0, double G0, double m, double m1,
                                      double m2) const {
  const double p_f = pCM(m, m1, m2);

  return G0 * p_f * blatt_weisskopf_sqr(p_f, L_) *
         post_ff_sqr(m, m0, threshold(), Lambda_) / (m * rho(m0));
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

double TwoBodyDecayDilepton::width(double m0, double G0, double m) const {
  if (m <= threshold()) {
    return 0;
  } else {
    /// dilepton decays: use width from \iref{Li:1996mi}, equation (19)
    const double ml = particle_types_[0]->mass();  // lepton mass
    const double ml_to_m_sqr = (ml / m) * (ml / m);
    const double m0_to_m_cubed = (m0 / m) * (m0 / m) * (m0 / m);
    return G0 * m0_to_m_cubed * std::sqrt(1. - 4. * ml_to_m_sqr) *
           (1. + 2. * ml_to_m_sqr);
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

unsigned int ThreeBodyDecay::particle_number() const { return 3; }

bool ThreeBodyDecay::has_particles(ParticleTypePtrList list) const {
  if (list.size() != particle_number()) {
    return false;
  }
  // compare sorted vectors (particle_types_ is already sorted)
  std::sort(list.begin(), list.end());
  return particle_types_[0] == list[0] && particle_types_[1] == list[1] &&
         particle_types_[2] == list[2];
}

double ThreeBodyDecay::width(double, double G0, double) const {
  return G0;  // use on-shell width
}

double ThreeBodyDecay::in_width(double, double G0, double, double,
                                double) const {
  return G0;  // use on-shell width
}

// ThreeBodyDecayDilepton

ThreeBodyDecayDilepton::ThreeBodyDecayDilepton(ParticleTypePtr mother,
                                               ParticleTypePtrList part_types,
                                               int l)
    : ThreeBodyDecay(part_types, l), tabulation_(nullptr), mother_(mother) {
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

double ThreeBodyDecayDilepton::diff_width(double m_par, double m_l,
                                          double m_dil, double m_other,
                                          ParticleTypePtr other,
                                          ParticleTypePtr t) {
  // check threshold
  if (m_par < m_dil + m_other) {
    return 0;
  }

  // abbreviations
  const double m_dil_sqr = m_dil * m_dil;
  const double m_par_sqr = m_par * m_par;
  const double m_par_cubed = m_par * m_par * m_par;
  const double m_other_sqr = m_other * m_other;
  const double ph_sp_factor = std::sqrt(1. - 4. * m_l * m_l / m_dil_sqr) *
                              (1. + 2. * m_l * m_l / m_dil_sqr);

  PdgCode pdg = t->pdgcode();
  if (pdg.is_meson()) {
    const ParticleType& photon = ParticleType::find(pdg::photon);
    switch (pdg.spin()) {
      case 0: /* pseudoscalars: π⁰, η, η' */ {
        // width for decay into 2γ
        const double gamma_2g = t->get_partial_width(m_par, photon, photon);
        double ff = em_form_factor_ps(pdg, m_dil);  // form factor
        /// see \iref{Landsberg:1986fd}, equation (3.8)
        return (4. * alpha / (3. * M_PI)) * gamma_2g / m_dil *
               pow_int(1. - m_dil / m_par * m_dil / m_par, 3) * ff * ff *
               ph_sp_factor;
      }
      case 2: /* vectors: ω, φ */ {
        // width for decay into Pγ with P = π,η
        const double gamma_pg = t->get_partial_width(m_par, *other, photon);
        double ff_sqr =
            em_form_factor_sqr_vec(pdg, m_dil);  // form factor squared
        /// see \iref{Landsberg:1986fd}, equation (3.4)
        const double n1 = m_par_sqr - m_other_sqr;
        const double rad = pow_int(1. + m_dil_sqr / n1, 2) -
                           4. * m_par_sqr * m_dil_sqr / (n1 * n1);
        if (rad < 0.) {
          assert(rad > -1E-5);
          return 0.;
        } else {
          return (2. * alpha / (3. * M_PI)) * gamma_pg / m_dil *
                 std::pow(rad, 3. / 2.) * ff_sqr * ph_sp_factor;
        }
      }
      default:
        throw std::runtime_error("Bad meson in ThreeBodyDecayDilepton: " +
                                 pdg.string());
    }
  } else if (pdg.is_baryon()) {
    switch (pdg.code()) {
      case pdg::Delta_p:
      case -pdg::Delta_p:
      case pdg::Delta_z:
      case -pdg::Delta_z: /* Δ⁺, Δ⁰ (and antiparticles) */ {
        /// see \iref{Krivoruchenko:2001hs}
        const double rad1 = (m_par + m_other) * (m_par + m_other) - m_dil_sqr;
        const double rad2 = (m_par - m_other) * (m_par - m_other) - m_dil_sqr;
        if (rad1 < 0.) {
          assert(rad1 > -1E-5);
          return 0.;
        } else if (rad2 < 0.) {
          assert(rad2 > -1E-5);
          return 0.;
        } else {
          const double t1 = alpha / 16. * (m_par + m_other) *
                            (m_par + m_other) / (m_par_cubed * m_other_sqr) *
                            std::sqrt(rad1);
          const double t2 = pow_int(std::sqrt(rad2), 3);
          const double ff = form_factor_delta(m_dil);
          const double gamma_vi = t1 * t2 * ff * ff;
          return 2. * alpha / (3. * M_PI) * gamma_vi / m_dil * ph_sp_factor;
        }
      }
      default:
        throw std::runtime_error("Bad baryon in ThreeBodyDecayDilepton: " +
                                 pdg.string());
    }
  } else {
    throw std::runtime_error("Non-hadron in ThreeBodyDecayDilepton: " +
                             pdg.string());
  }
}

double ThreeBodyDecayDilepton::width(double, double G0, double m) const {
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
    const double m_l = particle_types_[(non_lepton_position + 1) % 3]->mass();
    // mass of non-leptonic particle in final state
    const double m_other = particle_types_[non_lepton_position]->mass();

    // integrate differential width to obtain partial width
    double M0 = mother_->mass();
    double G0tot = mother_->width_at_pole();
    tabulation_ = make_unique<Tabulation>(
        m_other + 2 * m_l, M0 + 10 * G0tot, num_tab_pts, [&](double m_parent) {
          const double bottom = 2 * m_l;
          const double top = m_parent - m_other;
          if (top < bottom) {  // numerical problems at lower bound
            return 0.;
          }
          return integrate(bottom, top,
                           [&](double m_dil) {
                             return diff_width(
                                 m_parent, m_l, m_dil, m_other,
                                 particle_types_[non_lepton_position], mother_);
                           })
              .value();
        });
  }

  return tabulation_->get_value_linear(m, Extrapolation::Const);
}

}  // namespace smash
