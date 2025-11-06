/*
 *
 *    Copyright (c) 2018-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/crosssections.h"

#include "smash/clebschgordan.h"
#include "smash/constants.h"
#include "smash/logging.h"
#include "smash/parametrizations.h"
#include "smash/pow.h"

namespace smash {
static constexpr int LCrossSections = LogArea::CrossSections::id;
static constexpr int LScatterAction = LogArea::ScatterAction::id;

/**
 * Helper function:
 * Calculate the detailed balance factor R such that
 * \f[ R = \sigma(AB \to CD) / \sigma(CD \to AB) \f]
 * where \f$ A, B, C, D \f$ are stable.
 */
static double detailed_balance_factor_stable(double s, const ParticleType& a,
                                             const ParticleType& b,
                                             const ParticleType& c,
                                             const ParticleType& d) {
  double spin_factor = (c.spin() + 1) * (d.spin() + 1);
  spin_factor /= (a.spin() + 1) * (b.spin() + 1);
  double symmetry_factor = (1 + (a == b));
  symmetry_factor /= (1 + (c == d));
  const double momentum_factor = pCM_sqr_from_s(s, c.mass(), d.mass()) /
                                 pCM_sqr_from_s(s, a.mass(), b.mass());
  return spin_factor * symmetry_factor * momentum_factor;
}

/**
 * Helper function:
 * Calculate the detailed balance factor R such that
 * \f[ R = \sigma(AB \to CD) / \sigma(CD \to AB) \f]
 * where \f$A\f$ is unstable, \f$B\f$ is a kaon and \f$C, D\f$ are stable.
 */
static double detailed_balance_factor_RK(double sqrts, double pcm,
                                         const ParticleType& a,
                                         const ParticleType& b,
                                         const ParticleType& c,
                                         const ParticleType& d) {
  assert(!a.is_stable());
  assert(b.pdgcode().is_kaon());
  double spin_factor = (c.spin() + 1) * (d.spin() + 1);
  spin_factor /= (a.spin() + 1) * (b.spin() + 1);
  double symmetry_factor = (1 + (a == b));
  symmetry_factor /= (1 + (c == d));
  const double momentum_factor =
      pCM_sqr(sqrts, c.mass(), d.mass()) /
      (pcm * a.iso_multiplet()->get_integral_RK(sqrts));
  return spin_factor * symmetry_factor * momentum_factor;
}

/**
 * Helper function:
 * Calculate the detailed balance factor R such that
 * \f[ R = \sigma(AB \to CD) / \sigma(CD \to AB) \f]
 * where \f$A\f$ and \f$B\f$ are unstable, and \f$C\f$ and \f$D\f$ are stable.
 */
static double detailed_balance_factor_RR(double sqrts, double pcm,
                                         const ParticleType& a,
                                         const ParticleType& b,
                                         const ParticleType& c,
                                         const ParticleType& d) {
  assert(!a.is_stable());
  assert(!b.is_stable());
  double spin_factor = (c.spin() + 1) * (d.spin() + 1);
  spin_factor /= (a.spin() + 1) * (b.spin() + 1);
  double symmetry_factor = (1 + (a == b));
  symmetry_factor /= (1 + (c == d));
  const double momentum_factor =
      pCM_sqr(sqrts, c.mass(), d.mass()) /
      (pcm * a.iso_multiplet()->get_integral_RR(b.iso_multiplet(), sqrts));
  return spin_factor * symmetry_factor * momentum_factor;
}

/**
 * Helper function:
 * Append a list of processes to another (main) list of processes.
 */
static void append_list(CollisionBranchList& main_list,
                        CollisionBranchList in_list, double weight = 1.) {
  main_list.reserve(main_list.size() + in_list.size());
  for (auto& proc : in_list) {
    proc->set_weight(proc->weight() * weight);
    main_list.emplace_back(std::move(proc));
  }
}

/**
 * Helper function:
 * Shift the energy of a collision for AQM rescaled cross sections.
 *
 * \param[in] mandelstam_s the rest frame total energy squared
 * \param[in] m1 effective mass of incoming first particle
 * \param[in] m2 effective mass of incoming second particle
 * \param[in] m1_ref mass of the first AQM reference
 * \param[in] m2_ref mass of the second AQM reference
 * \return the shifted center of mass energy squared
 */
static double effective_AQM_s(double mandelstam_s, double m1, double m2,
                              double m1_ref, double m2_ref) {
  const double eff_sqrt_s = std::sqrt(mandelstam_s) - m1 - m2 + m1_ref + m2_ref;
  return eff_sqrt_s * eff_sqrt_s;
}

CrossSections::CrossSections(const ParticleList& incoming_particles,
                             const double sqrt_s,
                             const std::pair<FourVector, FourVector> potentials)
    : incoming_particles_(incoming_particles),
      sqrt_s_(sqrt_s),
      potentials_(potentials),
      is_BBbar_pair_(incoming_particles_[0].type().is_baryon() &&
                     incoming_particles_[1].type().is_baryon() &&
                     incoming_particles_[0].type().antiparticle_sign() ==
                         -incoming_particles_[1].type().antiparticle_sign()),
      is_NNbar_pair_(
          incoming_particles_[0].type().is_nucleon() &&
          incoming_particles_[1].pdgcode() ==
              incoming_particles_[0].type().get_antiparticle()->pdgcode()) {}

CollisionBranchList CrossSections::generate_collision_list(
    const ScatterActionsFinderParameters& finder_parameters,
    StringProcess* string_process) const {
  CollisionBranchList process_list;
  const ParticleType& t1 = incoming_particles_[0].type();
  const ParticleType& t2 = incoming_particles_[1].type();

  double p_pythia = 0.;
  if (finder_parameters.strings_with_probability) {
    p_pythia = string_probability(finder_parameters);
  }

  /* Elastic collisions between two nucleons with sqrt_s below
   * low_snn_cut can not happen. */
  const bool reject_by_nucleon_elastic_cutoff =
      t1.is_nucleon() && t2.is_nucleon() &&
      t1.antiparticle_sign() == t2.antiparticle_sign() &&
      sqrt_s_ < finder_parameters.low_snn_cut;
  bool incl_elastic =
      finder_parameters.included_2to2[IncludedReactions::Elastic];
  if (incl_elastic && !reject_by_nucleon_elastic_cutoff) {
    process_list.emplace_back(elastic(finder_parameters));
  }
  if (incoming_particles_[0].is_core() != incoming_particles_[1].is_core()) {
    return process_list;
  }
  if (p_pythia > 0.) {
    /* String-excitation cross section =
     * Parametrized total cross - the contributions
     * from all other present channels. */
    const double sig_current = sum_xs_of(process_list);
    const double sig_string = std::max(
        0., finder_parameters.scale_xs * high_energy(finder_parameters) -
                sig_current);
    append_list(
        process_list,
        string_excitation(sig_string, string_process, finder_parameters),
        p_pythia);
    append_list(process_list, rare_two_to_two(),
                p_pythia * finder_parameters.scale_xs);
  }
  if (p_pythia < 1.) {
    if (finder_parameters.two_to_one) {
      // resonance formation (2->1)
      append_list(process_list, two_to_one(),
                  (1. - p_pythia) * finder_parameters.scale_xs);
    }
    if (finder_parameters.included_2to2.any()) {
      // 2->2 (inelastic)
      append_list(
          process_list,
          two_to_two(finder_parameters.included_2to2,
                     finder_parameters.transition_high_energy.KN_offset),
          (1. - p_pythia) * finder_parameters.scale_xs);
    }
    if (finder_parameters
            .included_multi[IncludedMultiParticleReactions::Deuteron_3to2] ==
        1) {
      // 2->3 (deuterons only 2-to-3 reaction at the moment)
      append_list(process_list, two_to_three(),
                  (1. - p_pythia) * finder_parameters.scale_xs);
    }
    if (finder_parameters
            .included_multi[IncludedMultiParticleReactions::A3_Nuclei_4to2] ==
        1) {
      // 2->4
      append_list(process_list, two_to_four(),
                  (1. - p_pythia) * finder_parameters.scale_xs);
    }
  }
  if (finder_parameters.nnbar_treatment == NNbarTreatment::TwoToFive &&
      is_NNbar_pair_) {
    // NNbar directly to 5 pions (2-to-5)
    process_list.emplace_back(NNbar_to_5pi(finder_parameters.scale_xs));
  }

  /* NNbar annihilation thru NNbar → ρh₁(1170); combined with the decays
   * ρ → ππ and h₁(1170) → πρ, this gives a final state of 5 pions.
   * Only use in cases when detailed balance MUST happen, i.e. in a box! */
  if (finder_parameters.nnbar_treatment == NNbarTreatment::Resonances) {
    if (is_NNbar_pair_) {
      /* Has to be called after the other processes are already determined,
       * so that the sum of the cross sections includes all other processes. */
      process_list.emplace_back(NNbar_annihilation(sum_xs_of(process_list),
                                                   finder_parameters.scale_xs));
    } else {
      append_list(process_list, NNbar_creation(), finder_parameters.scale_xs);
    }
  }
  return process_list;
}

double CrossSections::parametrized_total(
    const ScatterActionsFinderParameters& finder_parameters) const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();
  double total_xs = 0.;
  if (pdg_a.is_baryon() && pdg_b.is_baryon() &&
      sqrt_s_ > finder_parameters.low_snn_cut) {
    if (pdg_a.antiparticle_sign() == pdg_b.antiparticle_sign()) {
      // NN
      total_xs = (pdg_a == pdg_b) ? pp_total(sqrt_s_ * sqrt_s_)
                                  : np_total(sqrt_s_ * sqrt_s_);
    } else {
      // NNbar
      total_xs = ppbar_total(sqrt_s_ * sqrt_s_);
    }
    total_xs *= finder_parameters.AQM_scaling_factor(pdg_a) *
                finder_parameters.AQM_scaling_factor(pdg_b);
  } else if ((pdg_a.is_baryon() && pdg_b.is_meson()) ||
             (pdg_a.is_meson() && pdg_b.is_baryon())) {
    const PdgCode& meson = pdg_a.is_meson() ? pdg_a : pdg_b;
    const PdgCode& baryon = pdg_a.is_meson() ? pdg_b : pdg_a;
    if (meson.is_kaon() && baryon.is_nucleon()) {
      if ((meson.code() == pdg::K_p && baryon.code() == pdg::p) ||
          (meson.code() == pdg::K_z && baryon.code() == pdg::n) ||
          (meson.code() == pdg::K_m && baryon.code() == -pdg::p) ||
          (meson.code() == pdg::Kbar_z && baryon.code() == -pdg::n)) {
        // K⁺p, K⁰n, and anti-processes
        total_xs = kplusp_total(sqrt_s_ * sqrt_s_);
      } else if ((meson.code() == pdg::K_p && baryon.code() == -pdg::p) ||
                 (meson.code() == pdg::K_z && baryon.code() == -pdg::n) ||
                 (meson.code() == pdg::K_m && baryon.code() == pdg::p) ||
                 (meson.code() == pdg::Kbar_z && baryon.code() == pdg::n)) {
        // K⁻p, K̅⁰n, and anti-processes
        total_xs = kminusp_total(sqrt_s_ * sqrt_s_);
      } else if ((meson.code() == pdg::K_p && baryon.code() == pdg::n) ||
                 (meson.code() == pdg::K_z && baryon.code() == pdg::p) ||
                 (meson.code() == pdg::K_m && baryon.code() == -pdg::n) ||
                 (meson.code() == pdg::Kbar_z && baryon.code() == -pdg::p)) {
        // K⁺n, K⁰p, and anti-processes
        total_xs = kplusn_total(sqrt_s_ * sqrt_s_);
      } else if ((meson.code() == pdg::K_p && baryon.code() == -pdg::n) ||
                 (meson.code() == pdg::K_z && baryon.code() == -pdg::p) ||
                 (meson.code() == pdg::K_m && baryon.code() == pdg::n) ||
                 (meson.code() == pdg::Kbar_z && baryon.code() == pdg::p)) {
        // K⁻n, K̅⁰p and anti-processes
        total_xs = kminusn_total(sqrt_s_ * sqrt_s_);
      }
    } else if (meson.is_pion() && baryon.is_nucleon()) {
      // π⁺(p,nbar), π⁻(n,pbar)
      if ((meson.code() == pdg::pi_p &&
           (baryon.code() == pdg::p || baryon.code() == -pdg::n)) ||
          (meson.code() == pdg::pi_m &&
           (baryon.code() == pdg::n || baryon.code() == -pdg::p))) {
        total_xs = piplusp_total(sqrt_s_);
      } else if (meson.code() == pdg::pi_z) {
        // π⁰N
        total_xs = 0.5 * (piplusp_total(sqrt_s_) + piminusp_total(sqrt_s_));
      } else {
        // π⁻(p,nbar), π⁺(n,pbar)
        total_xs = piminusp_total(sqrt_s_);
      }
    } else {
      // M*+B* goes to AQM high energy π⁻p
      total_xs = piminusp_high_energy(sqrt_s_ * sqrt_s_) *
                 finder_parameters.AQM_scaling_factor(pdg_a) *
                 finder_parameters.AQM_scaling_factor(pdg_b);
    }
  } else if (pdg_a.is_meson() && pdg_b.is_meson()) {
    if (pdg_a.is_pion() && pdg_b.is_pion()) {
      switch (pdg_a.isospin3() * pdg_b.isospin3() / 4) {
        // π⁺π⁻
        case -1:
          total_xs = pipluspiminus_total(sqrt_s_);
          break;
        case 0:
          // π⁰π⁰
          if (pdg_a.isospin3() + pdg_b.isospin3() == 0) {
            total_xs = pizeropizero_total(sqrt_s_);
          } else {
            // π⁺π⁰: similar to π⁺π⁻
            total_xs = pipluspiminus_total(sqrt_s_);
          }
          break;
        // π⁺π⁺ goes to π⁻p AQM
        case 1:
          total_xs = (2. / 3.) * piminusp_high_energy(sqrt_s_ * sqrt_s_);
          break;
        default:
          throw std::runtime_error("wrong isospin in ππ scattering");
      }
    } else {
      // M*+M* goes to AQM high energy π⁻p
      total_xs = (2. / 3.) * piminusp_high_energy(sqrt_s_ * sqrt_s_) *
                 finder_parameters.AQM_scaling_factor(pdg_a) *
                 finder_parameters.AQM_scaling_factor(pdg_b);
    }
  }
  return (total_xs + finder_parameters.additional_el_xs) *
         finder_parameters.scale_xs;
}

CollisionBranchPtr CrossSections::elastic(
    const ScatterActionsFinderParameters& finder_parameters) const {
  double elastic_xs = 0.;

  if (finder_parameters.elastic_parameter >= 0.) {
    // use constant elastic cross section from config file
    elastic_xs = finder_parameters.elastic_parameter;
  } else {
    // use parametrization
    elastic_xs = elastic_parametrization(finder_parameters);
  }
  /* when using a factor to scale the cross section and an additional
   * contribution to the elastic cross section, the contribution is added first
   * and then everything is scaled */
  return std::make_unique<CollisionBranch>(
      incoming_particles_[0].type(), incoming_particles_[1].type(),
      (elastic_xs + finder_parameters.additional_el_xs) *
          finder_parameters.scale_xs,
      ProcessType::Elastic);
}

CollisionBranchList CrossSections::rare_two_to_two() const {
  CollisionBranchList process_list;
  const ParticleData& data_a = incoming_particles_[0];
  const ParticleData& data_b = incoming_particles_[1];
  const auto& pdg_a = data_a.pdgcode();
  const auto& pdg_b = data_b.pdgcode();
  if ((pdg_a.is_nucleon() && pdg_b.is_pion()) ||
      (pdg_b.is_nucleon() && pdg_a.is_pion())) {
    process_list = npi_yk();
  }
  return process_list;
}

double CrossSections::elastic_parametrization(
    const ScatterActionsFinderParameters& finder_parameters) const {
  const bool use_AQM = finder_parameters.use_AQM;
  const double pipi_offset =
      finder_parameters.transition_high_energy.pipi_offset;
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();
  double elastic_xs = 0.0;
  if ((pdg_a.is_nucleon() && pdg_b.is_pion()) ||
      (pdg_b.is_nucleon() && pdg_a.is_pion())) {
    // Elastic Nucleon Pion Scattering
    elastic_xs = npi_el();
  } else if ((pdg_a.is_nucleon() && pdg_b.is_kaon()) ||
             (pdg_b.is_nucleon() && pdg_a.is_kaon())) {
    // Elastic Nucleon Kaon Scattering
    elastic_xs = nk_el();
  } else if (pdg_a.is_nucleon() && pdg_b.is_nucleon() &&
             pdg_a.antiparticle_sign() == pdg_b.antiparticle_sign()) {
    // Elastic Nucleon Nucleon Scattering
    elastic_xs = nn_el();
  } else if (pdg_a.is_nucleon() && pdg_b.is_nucleon() &&
             pdg_a.antiparticle_sign() == -pdg_b.antiparticle_sign()) {
    // Elastic Nucleon anti-Nucleon Scattering
    elastic_xs = ppbar_elastic(sqrt_s_ * sqrt_s_);
  } else if (pdg_a.is_nucleus() || pdg_b.is_nucleus()) {
    const PdgCode& pdg_nucleus = pdg_a.is_nucleus() ? pdg_a : pdg_b;
    const PdgCode& pdg_other = pdg_a.is_nucleus() ? pdg_b : pdg_a;
    const bool is_deuteron = pdg_nucleus.is_deuteron();  // d or anti-d
    if (is_deuteron && pdg_other.is_pion()) {
      // Elastic (Anti-)deuteron Pion Scattering
      elastic_xs = deuteron_pion_elastic(sqrt_s_ * sqrt_s_);
    } else if (is_deuteron && pdg_other.is_nucleon()) {
      // Elastic (Anti-)deuteron (Anti-)Nucleon Scattering
      elastic_xs = deuteron_nucleon_elastic(sqrt_s_ * sqrt_s_);
    }
  } else if (use_AQM) {
    const double m1 = incoming_particles_[0].effective_mass();
    const double m2 = incoming_particles_[1].effective_mass();
    const double s = sqrt_s_ * sqrt_s_;
    if (pdg_a.is_baryon() && pdg_b.is_baryon()) {
      elastic_xs = nn_el();  // valid also for annihilation
    } else if ((pdg_a.is_meson() && pdg_b.is_baryon()) ||
               (pdg_b.is_meson() && pdg_a.is_baryon())) {
      elastic_xs = piplusp_elastic_high_energy(s, m1, m2);
    } else if (pdg_a.is_meson() && pdg_b.is_meson()) {
      /* Special case: the pi+pi- elastic cross-section goes through resonances
       * at low sqrt_s, so we turn it off for this region so as not to destroy
       * the agreement with experimental data; this does not
       * apply to other pi pi cross-sections, which do not have any data */
      if (((pdg_a == pdg::pi_p && pdg_b == pdg::pi_m) ||
           (pdg_a == pdg::pi_m && pdg_b == pdg::pi_p)) &&
          (m1 + m2 + pipi_offset) > sqrt_s_) {
        elastic_xs = 0.0;
      } else {
        // meson-meson goes through scaling from π+p parametrization
        elastic_xs = 2. / 3. * piplusp_elastic_AQM(s, m1, m2);
      }
    }
    elastic_xs *= finder_parameters.AQM_scaling_factor(pdg_a) *
                  finder_parameters.AQM_scaling_factor(pdg_b);
  }
  return elastic_xs;
}

double CrossSections::nn_el() const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  // Use parametrized cross sections.
  double sig_el = -1.;
  const double s = sqrt_s_ * sqrt_s_;
  const bool is_NN_pair = pdg_a.is_nucleon() && pdg_b.is_nucleon();
  if (is_NN_pair) {
    if (is_BBbar_pair_) {
      // npbar and ppbar
      sig_el = ppbar_elastic(s);
    } else {
      sig_el = (pdg_a == pdg_b) ? pp_elastic(s) : np_elastic(s);
    }
  } else {
    // AQM - Additive Quark Model
    const double m1 = incoming_particles_[0].effective_mass();
    const double m2 = incoming_particles_[1].effective_mass();
    if (is_BBbar_pair_) {
      sig_el =
          ppbar_elastic(effective_AQM_s(s, m1, m2, nucleon_mass, nucleon_mass));
    } else {
      sig_el = pp_elastic_high_energy(s, m1, m2);
    }
  }

  if (sig_el > 0.) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a << " b=" << name_b
       << " j_a=" << pdg_a.spin() << " j_b=" << pdg_b.spin()
       << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

double CrossSections::npi_el() const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  const PdgCode& nucleon = pdg_a.is_nucleon() ? pdg_a : pdg_b;
  const PdgCode& pion = pdg_a.is_nucleon() ? pdg_b : pdg_a;
  assert(pion != nucleon);

  const double s = sqrt_s_ * sqrt_s_;

  double sig_el = 0.;
  switch (nucleon.code()) {
    case pdg::p:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    case pdg::n:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    case -pdg::p:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    case -pdg::n:
      switch (pion.code()) {
        case pdg::pi_p:
          sig_el = piplusp_elastic(s);
          break;
        case pdg::pi_m:
          sig_el = piminusp_elastic(s);
          break;
        case pdg::pi_z:
          sig_el = 0.5 * (piplusp_elastic(s) + piminusp_elastic(s));
          break;
      }
      break;
    default:
      throw std::runtime_error(
          "only the elastic cross section for proton-pion "
          "is implemented");
  }

  if (sig_el > 0) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a << " b=" << name_b
       << " j_a=" << pdg_a.spin() << " j_b=" << pdg_b.spin()
       << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

CollisionBranchList CrossSections::npi_yk() const {
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_nucleon = a.pdgcode().is_nucleon() ? a : b;
  const ParticleType& type_pion = a.pdgcode().is_nucleon() ? b : a;

  const auto pdg_nucleon = type_nucleon.pdgcode().code();
  const auto pdg_pion = type_pion.pdgcode().code();

  const double s = sqrt_s_ * sqrt_s_;

  /* The cross sections are paramectrized for four isospin channels. The
   * cross sections of the rest isospin channels are obtained using
   * Clebsch-Gordan coefficients */

  CollisionBranchList process_list;
  switch (pdg_nucleon) {
    case pdg::p: {
      switch (pdg_pion) {
        case pdg::pi_p: {
          const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          add_channel(
              process_list, [&] { return piplusp_sigmapluskplus_pdg(s); },
              sqrt_s_, type_K_p, type_Sigma_p);
          break;
        }
        case pdg::pi_m: {
          const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
          const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
          const auto& type_Lambda = ParticleType::find(pdg::Lambda);
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          add_channel(
              process_list, [&] { return piminusp_sigmaminuskplus_pdg(s); },
              sqrt_s_, type_K_p, type_Sigma_m);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_K_z, type_Sigma_z);
          add_channel(
              process_list, [&] { return piminusp_lambdak0_pdg(s); }, sqrt_s_,
              type_K_z, type_Lambda);
          break;
        }
        case pdg::pi_z: {
          const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
          const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
          const auto& type_Lambda = ParticleType::find(pdg::Lambda);
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          add_channel(
              process_list,
              [&] {
                return 0.5 * (piplusp_sigmapluskplus_pdg(s) -
                              piminusp_sigma0k0_res(s) +
                              piminusp_sigmaminuskplus_pdg(s));
              },
              sqrt_s_, type_K_p, type_Sigma_z);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_K_z, type_Sigma_p);
          add_channel(
              process_list, [&] { return 0.5 * piminusp_lambdak0_pdg(s); },
              sqrt_s_, type_K_p, type_Lambda);
          break;
        }
      }
      break;
    }
    case pdg::n: {
      switch (pdg_pion) {
        case pdg::pi_p: {
          const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
          const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
          const auto& type_Lambda = ParticleType::find(pdg::Lambda);
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          add_channel(
              process_list, [&] { return piminusp_sigmaminuskplus_pdg(s); },
              sqrt_s_, type_K_z, type_Sigma_p);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_K_p, type_Sigma_z);
          add_channel(
              process_list, [&] { return piminusp_lambdak0_pdg(s); }, sqrt_s_,
              type_K_p, type_Lambda);
          break;
        }
        case pdg::pi_m: {
          const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          add_channel(
              process_list, [&] { return piplusp_sigmapluskplus_pdg(s); },
              sqrt_s_, type_K_z, type_Sigma_m);
          break;
        }
        case pdg::pi_z: {
          const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
          const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
          const auto& type_Lambda = ParticleType::find(pdg::Lambda);
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          add_channel(
              process_list,
              [&] {
                return 0.5 * (piplusp_sigmapluskplus_pdg(s) -
                              piminusp_sigma0k0_res(s) +
                              piminusp_sigmaminuskplus_pdg(s));
              },
              sqrt_s_, type_K_z, type_Sigma_z);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_K_p, type_Sigma_m);
          add_channel(
              process_list, [&] { return 0.5 * piminusp_lambdak0_pdg(s); },
              sqrt_s_, type_K_z, type_Lambda);
          break;
        }
      }
      break;
    }
    case -pdg::p: {
      switch (pdg_pion) {
        case pdg::pi_p: {
          const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
          const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
          const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
          const auto& type_K_m = ParticleType::find(-pdg::K_p);
          const auto& type_Kbar_z = ParticleType::find(-pdg::K_z);
          add_channel(
              process_list, [&] { return piminusp_sigmaminuskplus_pdg(s); },
              sqrt_s_, type_K_m, type_Sigma_m_bar);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_Kbar_z, type_Sigma_z_bar);
          add_channel(
              process_list, [&] { return piminusp_lambdak0_pdg(s); }, sqrt_s_,
              type_Kbar_z, type_Lambda_bar);
          break;
        }
        case pdg::pi_m: {
          const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
          const auto& type_K_m = ParticleType::find(-pdg::K_p);
          add_channel(
              process_list, [&] { return piplusp_sigmapluskplus_pdg(s); },
              sqrt_s_, type_K_m, type_Sigma_p_bar);
          break;
        }
        case pdg::pi_z: {
          const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
          const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
          const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
          const auto& type_K_m = ParticleType::find(-pdg::K_p);
          const auto& type_Kbar_z = ParticleType::find(-pdg::K_z);
          add_channel(
              process_list,
              [&] {
                return 0.5 * (piplusp_sigmapluskplus_pdg(s) -
                              piminusp_sigma0k0_res(s) +
                              piminusp_sigmaminuskplus_pdg(s));
              },
              sqrt_s_, type_K_m, type_Sigma_z_bar);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_Kbar_z, type_Sigma_p_bar);
          add_channel(
              process_list, [&] { return 0.5 * piminusp_lambdak0_pdg(s); },
              sqrt_s_, type_K_m, type_Lambda_bar);
          break;
        }
      }
      break;
    }
    case -pdg::n: {
      switch (pdg_pion) {
        case pdg::pi_p: {
          const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
          const auto& type_Kbar_z = ParticleType::find(-pdg::K_z);
          add_channel(
              process_list, [&] { return piplusp_sigmapluskplus_pdg(s); },
              sqrt_s_, type_Kbar_z, type_Sigma_m_bar);
          break;
        }
        case pdg::pi_m: {
          const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
          const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
          const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
          const auto& type_K_m = ParticleType::find(-pdg::K_p);
          const auto& type_Kbar_z = ParticleType::find(-pdg::K_z);
          add_channel(
              process_list, [&] { return piminusp_sigmaminuskplus_pdg(s); },
              sqrt_s_, type_Kbar_z, type_Sigma_p_bar);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_K_m, type_Sigma_z_bar);
          add_channel(
              process_list, [&] { return piminusp_lambdak0_pdg(s); }, sqrt_s_,
              type_K_m, type_Lambda_bar);
          break;
        }
        case pdg::pi_z: {
          const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
          const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
          const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
          const auto& type_K_m = ParticleType::find(-pdg::K_p);
          const auto& type_Kbar_z = ParticleType::find(-pdg::K_z);
          add_channel(
              process_list,
              [&] {
                return 0.5 * (piplusp_sigmapluskplus_pdg(s) -
                              piminusp_sigma0k0_res(s) +
                              piminusp_sigmaminuskplus_pdg(s));
              },
              sqrt_s_, type_Kbar_z, type_Sigma_z_bar);
          add_channel(
              process_list, [&] { return piminusp_sigma0k0_res(s); }, sqrt_s_,
              type_K_m, type_Sigma_m_bar);
          add_channel(
              process_list, [&] { return 0.5 * piminusp_lambdak0_pdg(s); },
              sqrt_s_, type_Kbar_z, type_Lambda_bar);
          break;
        }
      }
      break;
    }
  }

  return process_list;
}

double CrossSections::nk_el() const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  const PdgCode& nucleon = pdg_a.is_nucleon() ? pdg_a : pdg_b;
  const PdgCode& kaon = pdg_a.is_nucleon() ? pdg_b : pdg_a;
  assert(kaon != nucleon);

  const double s = sqrt_s_ * sqrt_s_;

  double sig_el = 0.;
  switch (nucleon.code()) {
    case pdg::p:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kplusp_elastic_background(s);
          break;
        case pdg::K_m:
          sig_el = kminusp_elastic_background(s);
          break;
        case pdg::K_z:
          sig_el = k0p_elastic_background(s);
          break;
        case pdg::Kbar_z:
          sig_el = kbar0p_elastic_background(s);
          break;
      }
      break;
    case pdg::n:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kplusn_elastic_background(s);
          break;
        case pdg::K_m:
          sig_el = kminusn_elastic_background(s);
          break;
        case pdg::K_z:
          sig_el = k0n_elastic_background(s);
          break;
        case pdg::Kbar_z:
          sig_el = kbar0n_elastic_background(s);
          break;
      }
      break;
    case -pdg::p:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kminusp_elastic_background(s);
          break;
        case pdg::K_m:
          sig_el = kplusp_elastic_background(s);
          break;
        case pdg::K_z:
          sig_el = kbar0p_elastic_background(s);
          break;
        case pdg::Kbar_z:
          sig_el = k0p_elastic_background(s);
          break;
      }
      break;
    case -pdg::n:
      switch (kaon.code()) {
        case pdg::K_p:
          sig_el = kminusn_elastic_background(s);
          break;
        case pdg::K_m:
          sig_el = kplusn_elastic_background(s);
          break;
        case pdg::K_z:
          sig_el = kbar0n_elastic_background(s);
          break;
        case pdg::Kbar_z:
          sig_el = k0n_elastic_background(s);
          break;
      }
      break;
    default:
      throw std::runtime_error(
          "elastic cross section for antinucleon-kaon "
          "not implemented");
  }

  if (sig_el > 0) {
    return sig_el;
  } else {
    std::stringstream ss;
    const auto name_a = incoming_particles_[0].type().name();
    const auto name_b = incoming_particles_[1].type().name();
    ss << "problem in CrossSections::elastic: a=" << name_a << " b=" << name_b
       << " j_a=" << pdg_a.spin() << " j_b=" << pdg_b.spin()
       << " sigma=" << sig_el << " s=" << s;
    throw std::runtime_error(ss.str());
  }
}

CollisionBranchList CrossSections::two_to_one() const {
  CollisionBranchList resonance_process_list;
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  const double p_cm_sqr = pCM_sqr(sqrt_s_, m1, m2);

  ParticleTypePtrList possible_resonances =
      list_possible_resonances(&type_particle_a, &type_particle_b);

  // Find all the possible resonances
  for (const ParticleTypePtr type_resonance : possible_resonances) {
    double resonance_xsection = formation(*type_resonance, p_cm_sqr);

    // If cross section is non-negligible, add resonance to the list
    if (resonance_xsection > really_small) {
      resonance_process_list.push_back(std::make_unique<CollisionBranch>(
          *type_resonance, resonance_xsection, ProcessType::TwoToOne));
      logg[LCrossSections].debug("Found resonance: ", *type_resonance);
      logg[LCrossSections].debug(type_particle_a.name(), type_particle_b.name(),
                                 "->", type_resonance->name(),
                                 " at sqrt(s)[GeV] = ", sqrt_s_,
                                 " with xs[mb] = ", resonance_xsection);
    }
  }
  return resonance_process_list;
}

double CrossSections::formation(const ParticleType& type_resonance,
                                double cm_momentum_sqr) const {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  // Calculate partial in-width.
  const double partial_width = type_resonance.get_partial_in_width(
      sqrt_s_, incoming_particles_[0], incoming_particles_[1]);
  if (partial_width <= 0.) {
    return 0.;
  }

  assert(type_resonance.charge() ==
         type_particle_a.charge() + type_particle_b.charge());
  assert(type_resonance.baryon_number() ==
         type_particle_a.baryon_number() + type_particle_b.baryon_number());

  // Calculate spin factor
  const double spinfactor =
      static_cast<double>(type_resonance.spin() + 1) /
      ((type_particle_a.spin() + 1) * (type_particle_b.spin() + 1));
  const int sym_factor =
      (type_particle_a.pdgcode() == type_particle_b.pdgcode()) ? 2 : 1;
  /** Calculate resonance production cross section
   * using the Breit-Wigner distribution as probability amplitude.
   * See Eq. (176) in \iref{Buss:2011mx}. */
  return spinfactor * sym_factor * 2. * M_PI * M_PI / cm_momentum_sqr *
         type_resonance.spectral_function(sqrt_s_) * partial_width * hbarc *
         hbarc / fm2_mb;
}

CollisionBranchList CrossSections::two_to_two(
    const ReactionsBitSet& included_2to2, const double KN_offset) const {
  CollisionBranchList process_list;
  const ParticleData& data_a = incoming_particles_[0];
  const ParticleData& data_b = incoming_particles_[1];
  const ParticleType& type_a = data_a.type();
  const ParticleType& type_b = data_b.type();
  const auto& pdg_a = data_a.pdgcode();
  const auto& pdg_b = data_b.pdgcode();

  if (data_a.is_baryon() && data_b.is_baryon()) {
    if (pdg_a.is_nucleon() && pdg_b.is_nucleon() &&
        pdg_a.antiparticle_sign() == pdg_b.antiparticle_sign()) {
      // Nucleon Nucleon Scattering
      process_list = nn_xx(included_2to2);
    } else {
      // Baryon Baryon Scattering
      process_list = bb_xx_except_nn(included_2to2);
    }
  } else if ((type_a.is_baryon() && type_b.is_meson()) ||
             (type_a.is_meson() && type_b.is_baryon())) {
    // Baryon Meson Scattering
    if ((pdg_a.is_nucleon() && pdg_b.is_kaon()) ||
        (pdg_b.is_nucleon() && pdg_a.is_kaon())) {
      // Nucleon Kaon Scattering
      process_list = nk_xx(included_2to2, KN_offset);
    } else if ((pdg_a.is_hyperon() && pdg_b.is_pion()) ||
               (pdg_b.is_hyperon() && pdg_a.is_pion())) {
      // Hyperon Pion Scattering
      process_list = ypi_xx(included_2to2);
    } else if ((pdg_a.is_Delta() && pdg_b.is_kaon()) ||
               (pdg_b.is_Delta() && pdg_a.is_kaon())) {
      // Delta Kaon Scattering
      process_list = deltak_xx(included_2to2);
    }
  } else if (type_a.is_nucleus() || type_b.is_nucleus()) {
    if ((type_a.is_nucleon() && type_b.is_nucleus()) ||
        (type_b.is_nucleon() && type_a.is_nucleus())) {
      // Nucleon Deuteron and Nucleon d' Scattering
      process_list = dn_xx(included_2to2);
    } else if (((type_a.is_deuteron() || type_a.is_dprime()) &&
                pdg_b.is_pion()) ||
               ((type_b.is_deuteron() || type_b.is_dprime()) &&
                pdg_a.is_pion())) {
      // Pion Deuteron and Pion d' Scattering
      process_list = dpi_xx(included_2to2);
    }
  }
  return process_list;
}

CollisionBranchList CrossSections::two_to_three() const {
  CollisionBranchList process_list;
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();

  if ((type_a.is_deuteron() && type_b.pdgcode().is_pion()) ||
      (type_b.is_deuteron() && type_a.pdgcode().is_pion())) {
    const ParticleType& type_pi = type_a.pdgcode().is_pion() ? type_a : type_b;
    const ParticleType& type_nucleus = type_a.is_nucleus() ? type_a : type_b;

    if (type_nucleus.baryon_number() > 0) {
      // πd → πpn
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_n = ParticleType::find(pdg::n);

      process_list.push_back(std::make_unique<CollisionBranch>(
          type_pi, type_p, type_n, two_to_three_xs(type_a, type_b, sqrt_s_),
          ProcessType::TwoToThree));
    } else {
      // πd̅ → πp̅n̅
      const auto& type_anti_p = ParticleType::find(-pdg::p);
      const auto& type_anti_n = ParticleType::find(-pdg::n);

      process_list.push_back(std::make_unique<CollisionBranch>(
          type_pi, type_anti_p, type_anti_n,
          two_to_three_xs(type_a, type_b, sqrt_s_), ProcessType::TwoToThree));
    }
  }

  if ((type_a.is_nucleon() && type_b.is_deuteron()) ||
      (type_b.is_nucleon() && type_a.is_deuteron())) {
    const ParticleType& type_N = type_a.is_nucleon() ? type_a : type_b;
    const ParticleType& type_nucleus = type_a.is_deuteron() ? type_a : type_b;

    if (type_nucleus.baryon_number() > 0) {
      // Nd → Nnp, N̅d → N̅np
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_n = ParticleType::find(pdg::n);

      process_list.push_back(std::make_unique<CollisionBranch>(
          type_N, type_p, type_n, two_to_three_xs(type_a, type_b, sqrt_s_),
          ProcessType::TwoToThree));
    } else {
      // Nd̅ → Np̅n̅, N̅d̅ → N̅p̅n̅
      const auto& type_anti_p = ParticleType::find(-pdg::p);
      const auto& type_anti_n = ParticleType::find(-pdg::n);

      process_list.push_back(std::make_unique<CollisionBranch>(
          type_N, type_anti_p, type_anti_n,
          two_to_three_xs(type_a, type_b, sqrt_s_), ProcessType::TwoToThree));
    }
  }
  return process_list;
}

CollisionBranchList CrossSections::two_to_four() const {
  CollisionBranchList process_list;
  ParticleTypePtr type_nucleus = &(incoming_particles_[0].type());
  ParticleTypePtr type_catalyzer = &(incoming_particles_[1].type());
  if (!type_nucleus->is_nucleus()) {
    type_nucleus = &(incoming_particles_[1].type());
    type_catalyzer = &(incoming_particles_[0].type());
  }

  if (type_nucleus->is_nucleus() &&
      std::abs(type_nucleus->baryon_number()) == 3 &&
      (type_catalyzer->is_pion() || type_catalyzer->is_nucleon())) {
    const ParticleTypePtr type_p = ParticleType::try_find(pdg::p);
    const ParticleTypePtr type_n = ParticleType::try_find(pdg::n);
    const ParticleTypePtr type_anti_p = ParticleType::try_find(-pdg::p);
    const ParticleTypePtr type_anti_n = ParticleType::try_find(-pdg::n);
    const ParticleTypePtr type_la = ParticleType::try_find(pdg::Lambda);
    const ParticleTypePtr type_anti_la = ParticleType::try_find(-pdg::Lambda);

    // Find nucleus components
    ParticleTypePtrList components;
    components.reserve(3);
    const PdgCode nucleus_pdg = type_nucleus->pdgcode();
    for (int i = 0; i < nucleus_pdg.nucleus_p(); i++) {
      components.push_back(type_p);
    }
    for (int i = 0; i < nucleus_pdg.nucleus_n(); i++) {
      components.push_back(type_n);
    }
    for (int i = 0; i < nucleus_pdg.nucleus_ap(); i++) {
      components.push_back(type_anti_p);
    }
    for (int i = 0; i < nucleus_pdg.nucleus_an(); i++) {
      components.push_back(type_anti_n);
    }
    for (int i = 0; i < nucleus_pdg.nucleus_La(); i++) {
      components.push_back(type_la);
    }
    for (int i = 0; i < nucleus_pdg.nucleus_aLa(); i++) {
      components.push_back(type_anti_la);
    }
    if (sqrt_s_ > type_catalyzer->mass() + components[0]->mass() +
                      components[1]->mass() + components[2]->mass()) {
      process_list.push_back(std::make_unique<CollisionBranch>(
          *type_catalyzer, *(components[0]), *(components[1]), *(components[2]),
          two_to_four_xs(*type_nucleus, *type_catalyzer, sqrt_s_),
          ProcessType::TwoToFour));
    }
  }
  return process_list;
}

double CrossSections::d_pi_inelastic_xs(double pion_kinetic_energy) {
  const double x = pion_kinetic_energy;
  return x * (4.3 + 10.0 * x) / ((x - 0.16) * (x - 0.16) + 0.007);
}

double CrossSections::d_N_inelastic_xs(double N_kinetic_energy) {
  const double x = N_kinetic_energy;
  return x * (1.0 + 50 * x) / (x * x + 0.01) +
         4 * x / ((x - 0.008) * (x - 0.008) + 0.0004);
}

double CrossSections::d_aN_inelastic_xs(double aN_kinetic_energy) {
  return 55.0 / (aN_kinetic_energy + 0.17);
}

double CrossSections::two_to_three_xs(const ParticleType& type_a,
                                      const ParticleType& type_b,
                                      double sqrts) {
  double xs = 0.0;
  ParticleTypePtr type_nucleus = &type_a, type_catalyzer = &type_b;
  if (!type_nucleus->is_nucleus()) {
    type_nucleus = &type_b;
    type_catalyzer = &type_a;
  }

  bool nonzero_xs = type_nucleus->is_nucleus() &&
                    (type_catalyzer->is_pion() || type_catalyzer->is_nucleon());
  if (!nonzero_xs) {
    return 0.0;
  }

  const double md = type_nucleus->mass(), mcat = type_catalyzer->mass();
  const double Tkin = (sqrts * sqrts - (md + mcat) * (md + mcat)) / (2.0 * md);

  // Should normally never happen, but may be a useful safeguard
  if (Tkin <= 0.0) {
    return 0.0;
  }

  if (type_catalyzer->is_pion()) {
    xs = d_pi_inelastic_xs(Tkin);
  } else if (type_catalyzer->is_nucleon()) {
    if (type_nucleus->pdgcode().antiparticle_sign() ==
        type_catalyzer->pdgcode().antiparticle_sign()) {
      // Nd and N̅d̅
      xs = d_N_inelastic_xs(Tkin);
    } else {
      // N̅d and Nd̅
      xs = d_aN_inelastic_xs(Tkin);
    }
  }
  return xs;
}

double CrossSections::two_to_four_xs(const ParticleType& type_a,
                                     const ParticleType& type_b, double sqrts) {
  double xs = 0.0;
  ParticleTypePtr type_nucleus = &type_a, type_catalyzer = &type_b;
  if (!type_nucleus->is_nucleus()) {
    type_nucleus = &type_b;
    type_catalyzer = &type_a;
  }
  bool nonzero_xs = type_nucleus->is_nucleus() &&
                    (type_catalyzer->is_pion() || type_catalyzer->is_nucleon());
  if (!nonzero_xs) {
    return 0.0;
  }

  const double mA = type_nucleus->mass(), mcat = type_catalyzer->mass();
  const double Tkin = (sqrts * sqrts - (mA + mcat) * (mA + mcat)) / (2.0 * mA);
  const int A = type_nucleus->pdgcode().nucleus_A();
  // Should normally never happen, but may be a useful safeguard
  if (A != 3 || Tkin <= 0.0) {
    return 0.0;
  }

  if (type_catalyzer->is_pion()) {
    xs = A / 2. * d_pi_inelastic_xs(Tkin);
  } else if (type_catalyzer->is_nucleon()) {
    if (type_nucleus->pdgcode().antiparticle_sign() ==
        type_catalyzer->pdgcode().antiparticle_sign()) {
      // N + A, anti-N + anti-A
      xs = A / 2. * d_N_inelastic_xs(Tkin);
    } else {
      // N̅ + A and N + anti-A
      xs = A / 2. * d_aN_inelastic_xs(Tkin);
    }
  }
  return xs;
}

CollisionBranchList CrossSections::bb_xx_except_nn(
    const ReactionsBitSet& included_2to2) const {
  CollisionBranchList process_list;
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();

  bool same_sign = type_a.antiparticle_sign() == type_b.antiparticle_sign();
  bool any_nucleus = type_a.is_nucleus() || type_b.is_nucleus();
  if (!same_sign && !any_nucleus) {
    return process_list;
  }
  bool anti_particles = type_a.antiparticle_sign() == -1;
  if (type_a.is_nucleon() || type_b.is_nucleon()) {
    // N R → N N, N̅ R → N̅ N̅
    if (included_2to2[IncludedReactions::NN_to_NR] == 1) {
      process_list = bar_bar_to_nuc_nuc(anti_particles);
    }
  } else if (type_a.is_Delta() || type_b.is_Delta()) {
    // Δ R → N N, Δ̅ R → N̅ N̅
    if (included_2to2[IncludedReactions::NN_to_DR] == 1) {
      process_list = bar_bar_to_nuc_nuc(anti_particles);
    }
  }

  return process_list;
}

CollisionBranchList CrossSections::nn_xx(
    const ReactionsBitSet& included_2to2) const {
  CollisionBranchList process_list, channel_list;

  const double sqrts = sqrt_s_;

  /* Find whether colliding particles are nucleons or anti-nucleons;
   * adjust lists of produced particles. */
  bool both_antinucleons =
      (incoming_particles_[0].type().antiparticle_sign() == -1) &&
      (incoming_particles_[1].type().antiparticle_sign() == -1);
  const ParticleTypePtrList& nuc_or_anti_nuc =
      both_antinucleons ? ParticleType::list_anti_nucleons()
                        : ParticleType::list_nucleons();
  const ParticleTypePtrList& delta_or_anti_delta =
      both_antinucleons ? ParticleType::list_anti_Deltas()
                        : ParticleType::list_Deltas();
  // Find N N → N R channels.
  if (included_2to2[IncludedReactions::NN_to_NR] == 1) {
    channel_list = find_nn_xsection_from_type(
        ParticleType::list_baryon_resonances(), nuc_or_anti_nuc,
        [&sqrts](const ParticleType& type_res_1, const ParticleType&) {
          return type_res_1.iso_multiplet()->get_integral_NR(sqrts);
        });
    process_list.reserve(process_list.size() + channel_list.size());
    std::move(channel_list.begin(), channel_list.end(),
              std::inserter(process_list, process_list.end()));
    channel_list.clear();
  }

  // Find N N → Δ R channels.
  if (included_2to2[IncludedReactions::NN_to_DR] == 1) {
    channel_list = find_nn_xsection_from_type(
        ParticleType::list_baryon_resonances(), delta_or_anti_delta,
        [&sqrts](const ParticleType& type_res_1,
                 const ParticleType& type_res_2) {
          return type_res_1.iso_multiplet()->get_integral_RR(
              type_res_2.iso_multiplet(), sqrts);
        });
    process_list.reserve(process_list.size() + channel_list.size());
    std::move(channel_list.begin(), channel_list.end(),
              std::inserter(process_list, process_list.end()));
    channel_list.clear();
  }

  // Find N N → dπ and N̅ N̅→ d̅π channels.
  ParticleTypePtr deuteron = ParticleType::try_find(pdg::deuteron);
  ParticleTypePtr antideutron = ParticleType::try_find(pdg::antideuteron);
  ParticleTypePtr pim = ParticleType::try_find(pdg::pi_m);
  ParticleTypePtr pi0 = ParticleType::try_find(pdg::pi_z);
  ParticleTypePtr pip = ParticleType::try_find(pdg::pi_p);
  // Make sure all the necessary particle types are found
  if (deuteron && antideutron && pim && pi0 && pip &&
      included_2to2[IncludedReactions::PiDeuteron_to_NN] == 1) {
    const ParticleTypePtrList deutron_list = {deuteron};
    const ParticleTypePtrList antideutron_list = {antideutron};
    const ParticleTypePtrList pion_list = {pim, pi0, pip};
    channel_list = find_nn_xsection_from_type(
        (both_antinucleons ? antideutron_list : deutron_list), pion_list,
        [&sqrts](const ParticleType& type_res_1,
                 const ParticleType& type_res_2) {
          return pCM(sqrts, type_res_1.mass(), type_res_2.mass());
        });
    process_list.reserve(process_list.size() + channel_list.size());
    std::move(channel_list.begin(), channel_list.end(),
              std::inserter(process_list, process_list.end()));
    channel_list.clear();
  }

  return process_list;
}

CollisionBranchList CrossSections::nk_xx(const ReactionsBitSet& included_2to2,
                                         const double KN_offset) const {
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_nucleon = a.pdgcode().is_nucleon() ? a : b;
  const ParticleType& type_kaon = a.pdgcode().is_nucleon() ? b : a;

  const auto pdg_nucleon = type_nucleon.pdgcode().code();
  const auto pdg_kaon = type_kaon.pdgcode().code();

  const double s = sqrt_s_ * sqrt_s_;

  // Some variable declarations for frequently used quantities
  const auto sigma_kplusp = kplusp_inelastic_background(s);
  const auto sigma_kplusn = kplusn_inelastic_background(s);

  /* At high energy, the parametrization we use diverges from experimental
   * data. This cutoff represents the point where the AQM cross section
   * becomes smaller than this parametrization, so we cut it here, and fully
   * switch to AQM beyond this point. */
  const double KN_to_KDelta_cutoff = KN_offset +
                                     incoming_particles_[0].pole_mass() +
                                     incoming_particles_[1].pole_mass();

  bool incl_KN_to_KN = included_2to2[IncludedReactions::KN_to_KN] == 1;
  bool incl_KN_to_KDelta =
      included_2to2[IncludedReactions::KN_to_KDelta] == 1 &&
      sqrt_s_ < KN_to_KDelta_cutoff;
  bool incl_Strangeness_exchange =
      included_2to2[IncludedReactions::Strangeness_exchange] == 1;

  CollisionBranchList process_list;
  switch (pdg_kaon) {
    case pdg::K_m: {
      /* All inelastic K- N channels here are strangeness exchange, plus one
       * charge exchange. */
      switch (pdg_nucleon) {
        case pdg::p: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_m = ParticleType::find(pdg::pi_m);
            const auto& type_pi_p = ParticleType::find(pdg::pi_p);
            const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
            const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
            const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
            const auto& type_Lambda = ParticleType::find(pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusp_piminussigmaplus(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Sigma_p);
            add_channel(
                process_list, [&] { return kminusp_piplussigmaminus(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Sigma_m);
            add_channel(
                process_list, [&] { return kminusp_pi0sigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_z);
            add_channel(
                process_list, [&] { return kminusp_pi0lambda(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Lambda);
          }
          if (incl_KN_to_KN) {
            const auto& type_n = ParticleType::find(pdg::n);
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            add_channel(
                process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                type_Kbar_z, type_n);
          }
          break;
        }
        case pdg::n: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_m = ParticleType::find(pdg::pi_m);
            const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
            const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
            const auto& type_Lambda = ParticleType::find(pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Sigma_z);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_m);
            add_channel(
                process_list, [&] { return kminusn_piminuslambda(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Lambda);
          }
          break;
        }
        case -pdg::p: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            const auto& type_Delta_pp_bar = ParticleType::find(-pdg::Delta_pp);
            const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon,
                                            type_Kbar_z, type_Delta_pp_bar);
                },
                sqrt_s_, type_Kbar_z, type_Delta_pp_bar);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon, type_K_m,
                                            type_Delta_p_bar);
                },
                sqrt_s_, type_K_m, type_Delta_p_bar);
          }
          break;
        }
        case -pdg::n: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
            const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon,
                                            type_Kbar_z, type_Delta_p_bar);
                },
                sqrt_s_, type_Kbar_z, type_Delta_p_bar);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon, type_K_m,
                                            type_Delta_z_bar);
                },
                sqrt_s_, type_K_m, type_Delta_z_bar);
          }
          if (incl_KN_to_KN) {
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            const auto& type_p_bar = ParticleType::find(-pdg::p);
            add_channel(
                process_list, [&] { return kplusn_k0p(s); }, sqrt_s_,
                type_Kbar_z, type_p_bar);
          }
          break;
        }
      }
      break;
    }
    case pdg::K_p: {
      /* All inelastic channels are K+ N -> K Delta -> K pi N, with identical
       * cross section, weighted by the isospin factor. */
      switch (pdg_nucleon) {
        case pdg::p: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            const auto& type_Delta_pp = ParticleType::find(pdg::Delta_pp);
            const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_pp);
                },
                sqrt_s_, type_K_z, type_Delta_pp);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_p);
                },
                sqrt_s_, type_K_p, type_Delta_p);
          }
          break;
        }
        case pdg::n: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
            const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_p);
                },
                sqrt_s_, type_K_z, type_Delta_p);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_z);
                },
                sqrt_s_, type_K_p, type_Delta_z);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            const auto& type_p = ParticleType::find(pdg::p);
            add_channel(
                process_list, [&] { return kplusn_k0p(s); }, sqrt_s_, type_K_z,
                type_p);
          }
          break;
        }
        case -pdg::p: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_m = ParticleType::find(pdg::pi_m);
            const auto& type_pi_p = ParticleType::find(pdg::pi_p);
            const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
            const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
            const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
            const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusp_piminussigmaplus(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Sigma_p_bar);
            add_channel(
                process_list, [&] { return kminusp_piplussigmaminus(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Sigma_m_bar);
            add_channel(
                process_list, [&] { return kminusp_pi0sigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_z_bar);
            add_channel(
                process_list, [&] { return kminusp_pi0lambda(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Lambda_bar);
          }
          if (incl_KN_to_KN) {
            const auto& type_n_bar = ParticleType::find(-pdg::n);
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            add_channel(
                process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                type_K_z, type_n_bar);
          }
          break;
        }
        case -pdg::n: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_p = ParticleType::find(pdg::pi_p);
            const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
            const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
            const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Sigma_z_bar);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_m_bar);
            add_channel(
                process_list, [&] { return kminusn_piminuslambda(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Lambda_bar);
          }
          break;
        }
      }
      break;
    }
    case pdg::K_z: {
      // K+ and K0 have the same mass and spin, so their cross sections are
      // assumed to only differ in isospin factors. For the initial state, we
      // assume that K0 p is equivalent to K+ n and K0 n equivalent to K+ p,
      // like for the elastic background.

      switch (pdg_nucleon) {
        case pdg::p: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
            const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_p);
                },
                sqrt_s_, type_K_z, type_Delta_p);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_z);
                },
                sqrt_s_, type_K_p, type_Delta_z);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_n = ParticleType::find(pdg::n);
            add_channel(
                process_list,
                [&] {
                  // The isospin factor is 1, see the parametrizations
                  // tests.
                  return kplusn_k0p(s);
                },
                sqrt_s_, type_K_p, type_n);
          }
          break;
        }
        case pdg::n: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
            const auto& type_Delta_m = ParticleType::find(pdg::Delta_m);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_z);
                },
                sqrt_s_, type_K_z, type_Delta_z);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp *
                         kaon_nucleon_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_m);
                },
                sqrt_s_, type_K_p, type_Delta_m);
          }
          break;
        }
        case -pdg::p: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_m = ParticleType::find(pdg::pi_m);
            const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
            const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
            const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Sigma_z_bar);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_p_bar);
            add_channel(
                process_list, [&] { return kminusn_piminuslambda(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Lambda_bar);
          }
          break;
        }
        case -pdg::n: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_m = ParticleType::find(pdg::pi_m);
            const auto& type_pi_p = ParticleType::find(pdg::pi_p);
            const auto& type_Sigma_p_bar = ParticleType::find(-pdg::Sigma_p);
            const auto& type_Sigma_m_bar = ParticleType::find(-pdg::Sigma_m);
            const auto& type_Sigma_z_bar = ParticleType::find(-pdg::Sigma_z);
            const auto& type_Lambda_bar = ParticleType::find(-pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusp_piminussigmaplus(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Sigma_m_bar);
            add_channel(
                process_list, [&] { return kminusp_piplussigmaminus(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Sigma_p_bar);
            add_channel(
                process_list, [&] { return kminusp_pi0sigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_z_bar);
            add_channel(
                process_list, [&] { return kminusp_pi0lambda(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Lambda_bar);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_p_bar = ParticleType::find(-pdg::p);
            add_channel(
                process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                type_K_p, type_p_bar);
          }
          break;
        }
      }
      break;
    }
    case pdg::Kbar_z:
      switch (pdg_nucleon) {
        case pdg::p: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_p = ParticleType::find(pdg::pi_p);
            const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
            const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
            const auto& type_Lambda = ParticleType::find(pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_p);
            add_channel(
                process_list, [&] { return kminusn_piminussigma0(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Sigma_z);
            add_channel(
                process_list, [&] { return kminusn_piminuslambda(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Lambda);
          }
          break;
        }
        case pdg::n: {
          if (incl_Strangeness_exchange) {
            const auto& type_pi_z = ParticleType::find(pdg::pi_z);
            const auto& type_pi_m = ParticleType::find(pdg::pi_m);
            const auto& type_pi_p = ParticleType::find(pdg::pi_p);
            const auto& type_Sigma_p = ParticleType::find(pdg::Sigma_p);
            const auto& type_Sigma_m = ParticleType::find(pdg::Sigma_m);
            const auto& type_Sigma_z = ParticleType::find(pdg::Sigma_z);
            const auto& type_Lambda = ParticleType::find(pdg::Lambda);
            add_channel(
                process_list, [&] { return kminusp_piminussigmaplus(sqrt_s_); },
                sqrt_s_, type_pi_p, type_Sigma_m);
            add_channel(
                process_list, [&] { return kminusp_piplussigmaminus(sqrt_s_); },
                sqrt_s_, type_pi_m, type_Sigma_p);
            add_channel(
                process_list, [&] { return kminusp_pi0sigma0(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Sigma_z);
            add_channel(
                process_list, [&] { return kminusp_pi0lambda(sqrt_s_); },
                sqrt_s_, type_pi_z, type_Lambda);
          }
          if (incl_KN_to_KN) {
            const auto& type_p = ParticleType::find(pdg::p);
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            add_channel(
                process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                type_K_m, type_p);
          }
          break;
        }
        case -pdg::p: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_Kbar_z = type_kaon;
            const auto& type_Delta_bar_m = ParticleType::find(-pdg::Delta_p);
            const auto& type_Delta_bar_z = ParticleType::find(-pdg::Delta_z);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon,
                                            type_Kbar_z, type_Delta_bar_m);
                },
                sqrt_s_, type_Kbar_z, type_Delta_bar_m);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon, type_K_m,
                                            type_Delta_bar_z);
                },
                sqrt_s_, type_K_m, type_Delta_bar_z);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_n_bar = ParticleType::find(-pdg::n);
            add_channel(
                process_list,
                [&] {
                  // The isospin factor is 1, see the parametrizations
                  // tests.
                  return kplusn_k0p(s);
                },
                sqrt_s_, type_K_m, type_n_bar);
          }
          break;
        }
        case -pdg::n: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
            const auto& type_Delta_m_bar = ParticleType::find(-pdg::Delta_m);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon,
                                            type_Kbar_z, type_Delta_z_bar);
                },
                sqrt_s_, type_Kbar_z, type_Delta_z_bar);
            add_channel(
                process_list,
                [&] {
                  return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                            type_nucleon, type_kaon, type_K_m,
                                            type_Delta_m_bar);
                },
                sqrt_s_, type_K_m, type_Delta_m_bar);
          }
          break;
        }
      }
      break;
  }

  return process_list;
}

CollisionBranchList CrossSections::deltak_xx(
    const ReactionsBitSet& included_2to2) const {
  CollisionBranchList process_list;
  if (included_2to2[IncludedReactions::KN_to_KDelta] == 0) {
    return process_list;
  }
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_delta = a.pdgcode().is_Delta() ? a : b;
  const ParticleType& type_kaon = a.pdgcode().is_Delta() ? b : a;

  const auto pdg_delta = type_delta.pdgcode().code();
  const auto pdg_kaon = type_kaon.pdgcode().code();

  const double s = sqrt_s_ * sqrt_s_;
  const double pcm = cm_momentum();
  /* The cross sections are determined from the backward reactions via
   * detailed balance. The same isospin factors as for the backward reaction
   * are used. */
  switch (pack(pdg_delta, pdg_kaon)) {
    case pack(pdg::Delta_pp, pdg::K_z):
    case pack(pdg::Delta_p, pdg::K_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_p, type_K_p) *
                   kaon_nucleon_ratios.get_ratio(type_p, type_K_p, type_kaon,
                                                 type_delta) *
                   kplusp_inelastic_background(s);
          },
          sqrt_s_, type_p, type_K_p);
      break;
    }
    case pack(-pdg::Delta_pp, pdg::Kbar_z):
    case pack(-pdg::Delta_p, pdg::K_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_p_bar, type_K_m) *
                   kaon_nucleon_ratios.get_ratio(type_p_bar, type_K_m,
                                                 type_kaon, type_delta) *
                   kplusp_inelastic_background(s);
          },
          sqrt_s_, type_p_bar, type_K_m);
      break;
    }
    case pack(pdg::Delta_p, pdg::K_z):
    case pack(pdg::Delta_z, pdg::K_p): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_n, type_K_p) *
                   kaon_nucleon_ratios.get_ratio(type_n, type_K_p, type_kaon,
                                                 type_delta) *
                   kplusn_inelastic_background(s);
          },
          sqrt_s_, type_n, type_K_p);

      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_p, type_K_z) *
                   kaon_nucleon_ratios.get_ratio(type_p, type_K_z, type_kaon,
                                                 type_delta) *
                   kplusn_inelastic_background(s);
          },
          sqrt_s_, type_p, type_K_z);
      break;
    }
    case pack(-pdg::Delta_p, pdg::Kbar_z):
    case pack(-pdg::Delta_z, pdg::K_m): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_n_bar, type_K_m) *
                   kaon_nucleon_ratios.get_ratio(type_n_bar, type_K_m,
                                                 type_kaon, type_delta) *
                   kplusn_inelastic_background(s);
          },
          sqrt_s_, type_n_bar, type_K_m);

      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_p_bar,
                                              type_Kbar_z) *
                   kaon_nucleon_ratios.get_ratio(type_p_bar, type_Kbar_z,
                                                 type_kaon, type_delta) *
                   kplusn_inelastic_background(s);
          },
          sqrt_s_, type_p_bar, type_Kbar_z);
      break;
    }
    case pack(pdg::Delta_z, pdg::K_z):
    case pack(pdg::Delta_m, pdg::K_p): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_n, type_K_z) *
                   kaon_nucleon_ratios.get_ratio(type_n, type_K_z, type_kaon,
                                                 type_delta) *
                   kplusp_inelastic_background(s);
          },
          sqrt_s_, type_n, type_K_z);
      break;
    }
    case pack(-pdg::Delta_z, pdg::Kbar_z):
    case pack(-pdg::Delta_m, pdg::K_m): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                              type_kaon, type_n_bar,
                                              type_Kbar_z) *
                   kaon_nucleon_ratios.get_ratio(type_n_bar, type_Kbar_z,
                                                 type_kaon, type_delta) *
                   kplusp_inelastic_background(s);
          },
          sqrt_s_, type_n_bar, type_Kbar_z);
      break;
    }
    default:
      break;
  }

  return process_list;
}

CollisionBranchList CrossSections::ypi_xx(
    const ReactionsBitSet& included_2to2) const {
  CollisionBranchList process_list;
  if (included_2to2[IncludedReactions::Strangeness_exchange] == 0) {
    return process_list;
  }
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_hyperon = a.pdgcode().is_hyperon() ? a : b;
  const ParticleType& type_pion = a.pdgcode().is_hyperon() ? b : a;

  const auto pdg_hyperon = type_hyperon.pdgcode().code();
  const auto pdg_pion = type_pion.pdgcode().code();

  const double s = sqrt_s_ * sqrt_s_;

  switch (pack(pdg_hyperon, pdg_pion)) {
    case pack(pdg::Sigma_z, pdg::pi_m): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_K_m) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_n, type_K_m);
      break;
    }
    case pack(pdg::Sigma_z, pdg::pi_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_Kbar_z) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_p, type_Kbar_z);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_p): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_p) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_p);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_z) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_z);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_z): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_K_m) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_n, type_K_m);
      break;
    }
    case pack(pdg::Sigma_p, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_Kbar_z) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_p, type_Kbar_z);
      break;
    }
    case pack(-pdg::Sigma_m, pdg::pi_z): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_p) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_p);
      break;
    }
    case pack(-pdg::Sigma_p, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_z) *
                   kminusn_piminussigma0(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_z);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_m): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_K_m) *
                   kminusn_piminuslambda(sqrt_s_);
          },
          sqrt_s_, type_n, type_K_m);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_Kbar_z) *
                   kminusn_piminuslambda(sqrt_s_);
          },
          sqrt_s_, type_p, type_Kbar_z);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_p): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_p) *
                   kminusn_piminuslambda(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_p);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_z) *
                   kminusn_piminuslambda(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_z);
      break;
    }
    case pack(pdg::Sigma_z, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_K_m) *
                   kminusp_pi0sigma0(sqrt_s_);
          },
          sqrt_s_, type_p, type_K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_Kbar_z) *
                   kminusp_pi0sigma0(sqrt_s_);
          },
          sqrt_s_, type_n, type_Kbar_z);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_p) *
                   kminusp_pi0sigma0(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_z) *
                   kminusp_pi0sigma0(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_z);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_K_m) *
                   kminusp_piplussigmaminus(sqrt_s_);
          },
          sqrt_s_, type_p, type_K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_Kbar_z) *
                   kminusp_piminussigmaplus(sqrt_s_);
          },
          sqrt_s_, type_n, type_Kbar_z);
      break;
    }
    case pack(-pdg::Sigma_m, pdg::pi_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_p) *
                   kminusp_piplussigmaminus(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_z) *
                   kminusp_piminussigmaplus(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_z);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_K_m) *
                   kminusp_pi0lambda(sqrt_s_);
          },
          sqrt_s_, type_p, type_K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_Kbar_z) *
                   kminusp_pi0lambda(sqrt_s_);
          },
          sqrt_s_, type_n, type_Kbar_z);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_p) *
                   kminusp_pi0lambda(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_z) *
                   kminusp_pi0lambda(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_z);
      break;
    }
    case pack(pdg::Sigma_p, pdg::pi_m): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p, type_K_m) *
                   kminusp_piminussigmaplus(sqrt_s_);
          },
          sqrt_s_, type_p, type_K_m);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n, type_Kbar_z) *
                   kminusp_piplussigmaminus(sqrt_s_);
          },
          sqrt_s_, type_n, type_Kbar_z);
      break;
    }
    case pack(-pdg::Sigma_p, pdg::pi_p): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_p_bar, type_K_p) *
                   kminusp_piminussigmaplus(sqrt_s_);
          },
          sqrt_s_, type_p_bar, type_K_p);
      add_channel(
          process_list,
          [&] {
            return detailed_balance_factor_stable(s, type_hyperon, type_pion,
                                                  type_n_bar, type_K_z) *
                   kminusp_piplussigmaminus(sqrt_s_);
          },
          sqrt_s_, type_n_bar, type_K_z);
      break;
    }
    default:
      break;
  }

  return process_list;
}

double CrossSections::xs_dpi_dprimepi(const double sqrts, const double cm_mom,
                                      ParticleTypePtr produced_nucleus,
                                      const ParticleType& type_pi) {
  const double s = sqrts * sqrts;
  // same matrix element for πd and πd̅
  const double tmp = sqrts - pion_mass - deuteron_mass;
  // Matrix element is fit to match the inelastic pi+ d -> pi+ n p
  // cross-section from the Fig. 5 of [\iref{Arndt:1994bs}].
  const double matrix_element =
      295.5 + 2.862 / (0.00283735 + pow_int(sqrts - 2.181, 2)) +
      0.0672 / pow_int(tmp, 2) - 6.61753 / tmp;

  const double spin_factor =
      (produced_nucleus->spin() + 1) * (type_pi.spin() + 1);
  /* Isospin factor is always the same, so it is included into the
   * matrix element.
   * Symmetry factor is always 1 here.
   * The (hbarc)^2/16 pi factor is absorbed into matrix element. */
  double xsection = matrix_element * spin_factor / (s * cm_mom);
  if (produced_nucleus->is_stable()) {
    xsection *= pCM_from_s(s, type_pi.mass(), produced_nucleus->mass());
  } else {
    const double resonance_integral =
        produced_nucleus->iso_multiplet()->get_integral_piR(sqrts);
    xsection *= resonance_integral;
    logg[LScatterAction].debug("Resonance integral ", resonance_integral,
                               ", matrix element: ", matrix_element,
                               ", cm_momentum: ", cm_mom);
  }
  return xsection;
}

CollisionBranchList CrossSections::dpi_xx(
    const ReactionsBitSet& included_2to2) const {
  CollisionBranchList process_list;
  const double sqrts = sqrt_s_;
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();

  // pi d -> N N
  bool is_pid = (type_a.is_deuteron() && type_b.pdgcode().is_pion()) ||
                (type_b.is_deuteron() && type_a.pdgcode().is_pion());
  if (is_pid && included_2to2[IncludedReactions::PiDeuteron_to_NN] == 1) {
    const int baryon_number = type_a.baryon_number() + type_b.baryon_number();
    ParticleTypePtrList nuc = (baryon_number > 0)
                                  ? ParticleType::list_nucleons()
                                  : ParticleType::list_anti_nucleons();
    const double s = sqrt_s_ * sqrt_s_;
    for (ParticleTypePtr nuc_a : nuc) {
      for (ParticleTypePtr nuc_b : nuc) {
        if (type_a.charge() + type_b.charge() !=
            nuc_a->charge() + nuc_b->charge()) {
          continue;
        }
        // loop over total isospin
        for (const int twoI : I_tot_range(*nuc_a, *nuc_b)) {
          const double isospin_factor = isospin_clebsch_gordan_sqr_2to2(
              type_a, type_b, *nuc_a, *nuc_b, twoI);
          // If Clebsch-Gordan coefficient = 0, don't bother with the rest.
          if (std::abs(isospin_factor) < really_small) {
            continue;
          }

          // Calculate matrix element for inverse process.
          const double matrix_element =
              nn_to_resonance_matrix_element(sqrts, type_a, type_b, twoI);
          if (matrix_element <= 0.) {
            continue;
          }

          const double spin_factor = (nuc_a->spin() + 1) * (nuc_b->spin() + 1);
          const int sym_fac_in =
              (type_a.iso_multiplet() == type_b.iso_multiplet()) ? 2 : 1;
          const int sym_fac_out =
              (nuc_a->iso_multiplet() == nuc_b->iso_multiplet()) ? 2 : 1;
          double p_cm_final = pCM_from_s(s, nuc_a->mass(), nuc_b->mass());
          const double xsection = isospin_factor * spin_factor * sym_fac_in /
                                  sym_fac_out * p_cm_final * matrix_element /
                                  (s * cm_momentum());

          if (xsection > really_small) {
            process_list.push_back(std::make_unique<CollisionBranch>(
                *nuc_a, *nuc_b, xsection, ProcessType::TwoToTwo));
            logg[LScatterAction].debug(type_a.name(), type_b.name(), "->",
                                       nuc_a->name(), nuc_b->name(),
                                       " at sqrts [GeV] = ", sqrts,
                                       " with cs[mb] = ", xsection);
          }
        }
      }
    }
  }

  // pi d -> pi d' (effectively pi d -> pi p n)  AND reverse, pi d' -> pi d
  bool is_pid_or_pidprime = ((type_a.is_deuteron() || type_a.is_dprime()) &&
                             type_b.pdgcode().is_pion()) ||
                            ((type_b.is_deuteron() || type_b.is_dprime()) &&
                             type_a.pdgcode().is_pion());
  if (is_pid_or_pidprime &&
      included_2to2[IncludedReactions::PiDeuteron_to_pidprime] == 1) {
    const ParticleType& type_pi = type_a.pdgcode().is_pion() ? type_a : type_b;
    const ParticleType& type_nucleus = type_a.is_nucleus() ? type_a : type_b;
    ParticleTypePtrList nuclei = ParticleType::list_light_nuclei();
    for (ParticleTypePtr produced_nucleus : nuclei) {
      // Elastic collisions are treated in a different function
      if (produced_nucleus == &type_nucleus ||
          produced_nucleus->charge() != type_nucleus.charge() ||
          produced_nucleus->baryon_number() != type_nucleus.baryon_number()) {
        continue;
      }
      const double xsection =
          xs_dpi_dprimepi(sqrts, cm_momentum(), produced_nucleus, type_pi);
      process_list.push_back(std::make_unique<CollisionBranch>(
          type_pi, *produced_nucleus, xsection, ProcessType::TwoToTwo));
      logg[LScatterAction].debug(type_pi.name(), type_nucleus.name(), "→ ",
                                 type_pi.name(), produced_nucleus->name(),
                                 " at ", sqrts, " GeV, xs[mb] = ", xsection);
    }
  }
  return process_list;
}

double CrossSections::xs_dn_dprimen(const double sqrts, const double cm_mom,
                                    ParticleTypePtr produced_nucleus,
                                    const ParticleType& type_nucleus,
                                    const ParticleType& type_N) {
  const double s = sqrts * sqrts;
  double matrix_element = 0.0;
  double tmp = sqrts - nucleon_mass - deuteron_mass;
  assert(tmp >= 0.0);
  if (std::signbit(type_N.baryon_number()) ==
      std::signbit(type_nucleus.baryon_number())) {
    /** Nd → Nd', N̅d̅→ N̅d̅' and reverse:
     * Fit to match experimental cross-section Nd -> Nnp from
     * \cite Carlson1973. */
    matrix_element = 79.0474 / std::pow(tmp, 0.7897) + 654.596 * tmp;
  } else {
    /** N̅d →  N̅d', Nd̅→ Nd̅' and reverse:
     * Fit to roughly match experimental cross-section N̅d -> N̅ np from
     * \iref{Bizzarri:1973sp}. */
    matrix_element = 342.572 / std::pow(tmp, 0.6);
  }
  const double spin_factor =
      (produced_nucleus->spin() + 1) * (type_N.spin() + 1);
  /* Isospin factor is always the same, so it is included into matrix element
   * Symmetry factor is always 1 here
   * Absorb (hbarc)^2/16 pi factor into matrix element */
  double xsection = matrix_element * spin_factor / (s * cm_mom);
  if (produced_nucleus->is_stable()) {
    assert(!type_nucleus.is_stable());
    xsection *= pCM_from_s(s, type_N.mass(), produced_nucleus->mass());
  } else {
    assert(type_nucleus.is_stable());
    const double resonance_integral =
        produced_nucleus->iso_multiplet()->get_integral_NR(sqrts);
    xsection *= resonance_integral;
  }
  return xsection;
}

CollisionBranchList CrossSections::dn_xx(
    const ReactionsBitSet& included_2to2) const {
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();
  const ParticleType& type_N = type_a.is_nucleon() ? type_a : type_b;
  const ParticleType& type_nucleus = type_a.is_nucleus() ? type_a : type_b;
  CollisionBranchList process_list;
  if (included_2to2[IncludedReactions::NDeuteron_to_Ndprime] == 0) {
    return process_list;
  }
  ParticleTypePtrList nuclei = ParticleType::list_light_nuclei();

  for (ParticleTypePtr produced_nucleus : nuclei) {
    // No elastic collisions for now, respect conservation laws
    if (produced_nucleus == &type_nucleus ||
        produced_nucleus->charge() != type_nucleus.charge() ||
        produced_nucleus->baryon_number() != type_nucleus.baryon_number()) {
      continue;
    }
    const double xsection = xs_dn_dprimen(
        sqrt_s_, cm_momentum(), produced_nucleus, type_nucleus, type_N);
    process_list.push_back(std::make_unique<CollisionBranch>(
        type_N, *produced_nucleus, xsection, ProcessType::TwoToTwo));
    logg[LScatterAction].debug(type_N.name(), type_nucleus.name(), "→ ",
                               type_N.name(), produced_nucleus->name(), " at ",
                               sqrt_s_, " GeV, xs[mb] = ", xsection);
  }
  return process_list;
}

CollisionBranchList CrossSections::string_excitation(
    double total_string_xs, StringProcess* string_process,
    const ScatterActionsFinderParameters& finder_parameters) const {
  if (!string_process) {
    throw std::runtime_error("string_process should be initialized.");
  }

  CollisionBranchList channel_list;
  if (total_string_xs <= 0.) {
    return channel_list;
  }

  double mandelstam_s = sqrt_s_ * sqrt_s_;
  /* Get mapped PDG id for evaluation of the parametrized cross sections
   * for diffractive processes.
   * This must be rescaled according to additive quark model
   * in the case of exotic hadrons.
   * Also calculate the multiplicative factor for AQM
   * based on the quark contents. */
  std::array<int, 2> pdgid;
  double AQM_scaling = 1.;
  for (int i = 0; i < 2; i++) {
    PdgCode pdg = incoming_particles_[i].type().pdgcode();
    pdgid[i] = StringProcess::pdg_map_for_pythia(pdg);
    AQM_scaling *= finder_parameters.AQM_scaling_factor(pdg);
  }

  /* Determine if the initial state is a baryon-antibaryon pair,
   * which can annihilate. */
  bool can_annihilate = false;
  if (is_BBbar_pair_) {
    int n_q_types = 5;  // u, d, s, c, b
    for (int iq = 1; iq <= n_q_types; iq++) {
      std::array<int, 2> nquark;
      for (int i = 0; i < 2; i++) {
        nquark[i] =
            incoming_particles_[i].type().pdgcode().net_quark_number(iq);
      }
      if (nquark[0] != 0 && nquark[1] != 0) {
        can_annihilate = true;
        break;
      }
    }
  }

  /* The case for baryon/anti-baryon annihilation is treated separately,
   * as in this case we use only one way to break up the particles, namely
   * into 2 mesonic strings of equal mass after annihilating one quark-
   * anti-quark pair. See StringProcess::next_BBbarAnn() */
  double sig_annihilation = 0.0;
  if (can_annihilate) {
    /* In the case of baryon-antibaryon pair,
     * the parametrized cross section for annihilation will be added.
     * See xs_ppbar_annihilation(). */
    mandelstam_s = effective_AQM_s(
        mandelstam_s, incoming_particles_[0].effective_mass(),
        incoming_particles_[1].effective_mass(), nucleon_mass, nucleon_mass);
    double xs_param = xs_ppbar_annihilation(mandelstam_s);
    if (finder_parameters.use_AQM) {
      xs_param *= AQM_scaling;
    }
    sig_annihilation = std::min(total_string_xs, xs_param);
  }

  /* Total parametrized cross-section (I) and pythia-produced total
   * cross-section (II) do not necessarily coincide. If I > II then
   * non-diffractive cross-section is reinforced to get I == II.
   * If I < II then partial cross-sections are drained one-by-one
   * to reduce II until I == II:
   * first non-diffractive, then double-diffractive, then
   * single-diffractive AB->AX and AB->XB in equal proportion.
   * The way it is done here is not unique. I (ryu) think that at high energy
   * collision this is not an issue, but at sqrt_s < 10 GeV it may
   * matter. */
  std::array<double, 3> xs = string_process->cross_sections_diffractive(
      pdgid[0], pdgid[1], std::sqrt(mandelstam_s));
  if (finder_parameters.use_AQM) {
    for (int ip = 0; ip < 3; ip++) {
      xs[ip] *= AQM_scaling;
    }
  }
  double single_diffr_AX = xs[0], single_diffr_XB = xs[1], double_diffr = xs[2];
  double single_diffr = single_diffr_AX + single_diffr_XB;
  double diffractive = single_diffr + double_diffr;

  const double nondiffractive_all =
      std::max(0., total_string_xs - sig_annihilation - diffractive);
  diffractive = total_string_xs - sig_annihilation - nondiffractive_all;
  double_diffr = std::max(0., diffractive - single_diffr);
  const double a = (diffractive - double_diffr) / single_diffr;
  single_diffr_AX *= a;
  single_diffr_XB *= a;
  assert(std::abs(single_diffr_AX + single_diffr_XB + double_diffr +
                  sig_annihilation + nondiffractive_all - total_string_xs) <
         1.e-6);

  double nondiffractive_soft = 0.;
  double nondiffractive_hard = 0.;
  if (nondiffractive_all > 0.) {
    /* Hard string process is added by hard cross section
     * in conjunction with multipartion interaction picture
     * \iref{Sjostrand:1987su}. */
    const double hard_xsec = AQM_scaling * string_hard_cross_section();
    nondiffractive_soft =
        nondiffractive_all * std::exp(-hard_xsec / nondiffractive_all);
    nondiffractive_hard = nondiffractive_all - nondiffractive_soft;
  }
  logg[LCrossSections].debug("String cross sections [mb] are");
  logg[LCrossSections].debug("Single-diffractive AB->AX: ", single_diffr_AX);
  logg[LCrossSections].debug("Single-diffractive AB->XB: ", single_diffr_XB);
  logg[LCrossSections].debug("Double-diffractive AB->XX: ", double_diffr);
  logg[LCrossSections].debug("Soft non-diffractive: ", nondiffractive_soft);
  logg[LCrossSections].debug("Hard non-diffractive: ", nondiffractive_hard);
  logg[LCrossSections].debug("B-Bbar annihilation: ", sig_annihilation);

  /* cross section of soft string excitation including annihilation */
  const double sig_string_soft = total_string_xs - nondiffractive_hard;

  /* fill the list of process channels */
  if (sig_string_soft > 0.) {
    channel_list.push_back(std::make_unique<CollisionBranch>(
        single_diffr_AX, ProcessType::StringSoftSingleDiffractiveAX));
    channel_list.push_back(std::make_unique<CollisionBranch>(
        single_diffr_XB, ProcessType::StringSoftSingleDiffractiveXB));
    channel_list.push_back(std::make_unique<CollisionBranch>(
        double_diffr, ProcessType::StringSoftDoubleDiffractive));
    channel_list.push_back(std::make_unique<CollisionBranch>(
        nondiffractive_soft, ProcessType::StringSoftNonDiffractive));
    if (can_annihilate) {
      channel_list.push_back(std::make_unique<CollisionBranch>(
          sig_annihilation, ProcessType::StringSoftAnnihilation));
    }
  }
  if (nondiffractive_hard > 0.) {
    channel_list.push_back(std::make_unique<CollisionBranch>(
        nondiffractive_hard, ProcessType::StringHard));
  }
  return channel_list;
}

double CrossSections::high_energy(
    const ScatterActionsFinderParameters& finder_parameters) const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  const double s = sqrt_s_ * sqrt_s_;
  double xs = 0.;

  // Currently all BB collisions use the nucleon-nucleon parametrizations.
  if (pdg_a.is_baryon() && pdg_b.is_baryon()) {
    const double eff_s = effective_AQM_s(
        s, incoming_particles_[0].effective_mass(),
        incoming_particles_[1].effective_mass(), nucleon_mass, nucleon_mass);
    if (pdg_a == pdg_b) {
      xs = pp_high_energy(eff_s);  // pp, nn
    } else if (pdg_a.antiparticle_sign() * pdg_b.antiparticle_sign() == 1) {
      xs = np_high_energy(eff_s);  // np, nbarpbar
    } else if (pdg_a.antiparticle_sign() * pdg_b.antiparticle_sign() == -1) {
      /* In the case of baryon-antibaryon interactions,
       * the low-energy cross section must be involved
       * due to annihilation processes (via strings). */
      double xs_l = ppbar_total(eff_s);
      double xs_h = 0.;
      if (pdg_a.is_antiparticle_of(pdg_b)) {
        xs_h = ppbar_high_energy(eff_s);  // ppbar, nnbar
      } else {
        xs_h = npbar_high_energy(eff_s);  // npbar, nbarp
      }
      /* Transition between low and high energy is set to be consistent with
       * that defined in string_probability(). */
      auto [region_lower, region_upper] =
          finder_parameters.transition_high_energy.sqrts_range_NN;
      double prob_high = probability_transit_high(region_lower, region_upper);
      xs = xs_l * (1. - prob_high) + xs_h * prob_high;
    }
  }

  // Pion nucleon interaction / baryon-meson
  if ((pdg_a == pdg::pi_p && pdg_b == pdg::p) ||
      (pdg_b == pdg::pi_p && pdg_a == pdg::p) ||
      (pdg_a == pdg::pi_m && pdg_b == pdg::n) ||
      (pdg_b == pdg::pi_m && pdg_a == pdg::n)) {
    xs = piplusp_high_energy(s);  // pi+ p, pi- n
  } else if ((pdg_a == pdg::pi_m && pdg_b == pdg::p) ||
             (pdg_b == pdg::pi_m && pdg_a == pdg::p) ||
             (pdg_a == pdg::pi_p && pdg_b == pdg::n) ||
             (pdg_b == pdg::pi_p && pdg_a == pdg::n)) {
    xs = piminusp_high_energy(s);  // pi- p, pi+ n
  } else if ((pdg_a.is_meson() && pdg_b.is_baryon()) ||
             (pdg_b.is_meson() && pdg_a.is_baryon())) {
    xs = piminusp_high_energy(s);  // default for baryon-meson
  }

  /* Meson-meson interaction goes through AQM from pi+p,
   * see user guide "Use_AQM"*/
  if (pdg_a.is_meson() && pdg_b.is_meson()) {
    /* 2/3 factor since difference of 1 meson between meson-meson
     * and baryon-meson */
    xs = 2. / 3. * piplusp_high_energy(s);
  }

  // AQM scaling for cross-sections
  xs *= finder_parameters.AQM_scaling_factor(pdg_a) *
        finder_parameters.AQM_scaling_factor(pdg_b);

  return xs;
}

double CrossSections::string_hard_cross_section() const {
  double cross_sec = 0.;
  /* Hard strings can only be excited if the lower cutoff by
   * Pythia is fulfilled */
  if (sqrt_s_ <= minimum_sqrts_pythia_can_handle) {
    return cross_sec;
  }
  const ParticleData& data_a = incoming_particles_[0];
  const ParticleData& data_b = incoming_particles_[1];

  if (data_a.is_baryon() && data_b.is_baryon()) {
    // Nucleon-nucleon cross section is used for all baryon-baryon cases.
    const double eff_s =
        effective_AQM_s(sqrt_s_ * sqrt_s_, data_a.effective_mass(),
                        data_b.effective_mass(), nucleon_mass, nucleon_mass);
    cross_sec = NN_string_hard(eff_s);
  } else if (data_a.is_baryon() || data_b.is_baryon()) {
    // Nucleon-pion cross section is used for all baryon-meson cases.
    cross_sec = Npi_string_hard(sqrt_s_ * sqrt_s_);
  } else {
    // Pion-pion cross section is used for all meson-meson cases.
    cross_sec = pipi_string_hard(sqrt_s_ * sqrt_s_);
  }

  return cross_sec;
}

CollisionBranchPtr CrossSections::NNbar_to_5pi(const double scale_xs) const {
  const double s = sqrt_s_ * sqrt_s_;
  /* Use difference between total and elastic in order to conserve detailed
   * balance for all inelastoc NNbar processes. */
  const double nnbar_xsec = std::max(0., ppbar_total(s) - ppbar_elastic(s));
  logg[LCrossSections].debug("NNbar cross section for 2-to-5 is: ", nnbar_xsec);

  /* Make collision channel NNbar -> 5π (with same final state as resonance
   * approach). */
  const auto& type_piz = ParticleType::find(pdg::pi_z);
  const auto& type_pip = ParticleType::find(pdg::pi_p);
  const auto& type_pim = ParticleType::find(pdg::pi_m);
  return std::make_unique<CollisionBranch>(
      type_pip, type_pim, type_pip, type_pim, type_piz, nnbar_xsec * scale_xs,
      ProcessType::TwoToFive);
}

CollisionBranchPtr CrossSections::NNbar_annihilation(
    const double current_xs, const double scale_xs) const {
  /* Calculate NNbar cross section:
   * Parametrized total minus all other present channels.*/
  const double s = sqrt_s_ * sqrt_s_;
  double nnbar_xsec = std::max(0., ppbar_total(s) * scale_xs - current_xs);
  logg[LCrossSections].debug("NNbar cross section is: ", nnbar_xsec);
  // Make collision channel NNbar -> ρh₁(1170); eventually decays into 5π
  return std::make_unique<CollisionBranch>(ParticleType::find(pdg::h1),
                                           ParticleType::find(pdg::rho_z),
                                           nnbar_xsec, ProcessType::TwoToTwo);
}

CollisionBranchList CrossSections::NNbar_creation() const {
  CollisionBranchList channel_list;
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();
  if ((type_a.pdgcode() == pdg::rho_z && type_b.pdgcode() == pdg::h1) ||
      (type_a.pdgcode() == pdg::h1 && type_b.pdgcode() == pdg::rho_z)) {
    /* Calculate NNbar reverse cross section:
     * from reverse reaction (see NNbar_annihilation_cross_section).*/
    const double s = sqrt_s_ * sqrt_s_;
    const double pcm = cm_momentum();

    const auto& type_N = ParticleType::find(pdg::p);
    const auto& type_Nbar = ParticleType::find(-pdg::p);

    // Check available energy
    if (sqrt_s_ - 2 * type_N.mass() < 0) {
      return channel_list;
    }

    double xsection = detailed_balance_factor_RR(sqrt_s_, pcm, type_a, type_b,
                                                 type_N, type_Nbar) *
                      std::max(0., ppbar_total(s) - ppbar_elastic(s));
    logg[LCrossSections].debug("NNbar reverse cross section is: ", xsection);
    channel_list.push_back(std::make_unique<CollisionBranch>(
        type_N, type_Nbar, xsection, ProcessType::TwoToTwo));
    channel_list.push_back(std::make_unique<CollisionBranch>(
        ParticleType::find(pdg::n), ParticleType::find(-pdg::n), xsection,
        ProcessType::TwoToTwo));
  }
  return channel_list;
}

CollisionBranchList CrossSections::bar_bar_to_nuc_nuc(
    const bool is_anti_particles) const {
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();
  CollisionBranchList process_list;

  const double s = sqrt_s_ * sqrt_s_;
  // CM momentum in final state
  double p_cm_final = std::sqrt(s - 4. * nucleon_mass * nucleon_mass) / 2.;

  ParticleTypePtrList nuc_or_anti_nuc;
  if (is_anti_particles) {
    nuc_or_anti_nuc = ParticleType::list_anti_nucleons();
  } else {
    nuc_or_anti_nuc = ParticleType::list_nucleons();
  }

  // Loop over all nucleon or anti-nucleon charge states.
  for (ParticleTypePtr nuc_a : nuc_or_anti_nuc) {
    for (ParticleTypePtr nuc_b : nuc_or_anti_nuc) {
      /* Check for charge conservation. */
      if (type_a.charge() + type_b.charge() !=
          nuc_a->charge() + nuc_b->charge()) {
        continue;
      }
      // loop over total isospin
      for (const int twoI : I_tot_range(*nuc_a, *nuc_b)) {
        const double isospin_factor = isospin_clebsch_gordan_sqr_2to2(
            type_a, type_b, *nuc_a, *nuc_b, twoI);
        // If Clebsch-Gordan coefficient is zero, don't bother with the rest
        if (std::abs(isospin_factor) < really_small) {
          continue;
        }

        // Calculate matrix element for inverse process.
        const double matrix_element =
            nn_to_resonance_matrix_element(sqrt_s_, type_a, type_b, twoI);
        if (matrix_element <= 0.) {
          continue;
        }

        /** Cross section for 2->2 resonance absorption, obtained via detailed
         * balance from the inverse reaction.
         * See eqs. (B.6), (B.9) and (181) in \iref{Buss:2011mx}.
         * There are factors for spin, isospin and symmetry involved. */
        const double spin_factor = (nuc_a->spin() + 1) * (nuc_b->spin() + 1);
        const int sym_fac_in =
            (type_a.iso_multiplet() == type_b.iso_multiplet()) ? 2 : 1;
        const int sym_fac_out =
            (nuc_a->iso_multiplet() == nuc_b->iso_multiplet()) ? 2 : 1;
        const double xsection = isospin_factor * spin_factor * sym_fac_in /
                                sym_fac_out * p_cm_final * matrix_element /
                                (s * cm_momentum());

        if (xsection > really_small) {
          process_list.push_back(std::make_unique<CollisionBranch>(
              *nuc_a, *nuc_b, xsection, ProcessType::TwoToTwo));
          logg[LCrossSections].debug(
              "2->2 absorption with original particles: ", type_a, type_b);
        }
      }
    }
  }
  return process_list;
}

double CrossSections::nn_to_resonance_matrix_element(double sqrts,
                                                     const ParticleType& type_a,
                                                     const ParticleType& type_b,
                                                     const int twoI) {
  const double m_a = type_a.mass();
  const double m_b = type_b.mass();
  const double msqr = 2. * (m_a * m_a + m_b * m_b);
  /* If the c.m. energy is larger than the sum of the pole masses of the
   * outgoing particles plus three times of the sum of the widths plus 3 GeV,
   * the collision will be neglected.
   *
   * This can be problematic for some final-state cross sections, but at
   * energies that high strings are used anyway.
   */
  const double w_a = type_a.width_at_pole();
  const double w_b = type_b.width_at_pole();
  const double uplmt = m_a + m_b + 3.0 * (w_a + w_b) + 3.0;
  if (sqrts > uplmt) {
    return 0.;
  }
  /// NN → NΔ: fit sqrt(s)-dependence to OBE model [\iref{Dmitriev:1986st}]
  if (((type_a.is_Delta() && type_b.is_nucleon()) ||
       (type_b.is_Delta() && type_a.is_nucleon())) &&
      (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    return 68. / std::pow(sqrts - 1.104, 1.951);
    /** All other processes use a constant matrix element,
     *  similar to \iref{Bass:1998ca}, equ. (3.35). */
  } else if (((type_a.is_Nstar() && type_b.is_nucleon()) ||
              (type_b.is_Nstar() && type_a.is_nucleon())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → NN*
    if (twoI == 2) {
      return 4.5 / msqr;
    } else if (twoI == 0) {
      const double parametrization = 14. / msqr;
      /** pn → pnη cross section is known to be larger than the corresponding
       * pp → ppη cross section by a factor of 6.5 [\iref{Calen:1998vh}].
       * Since the eta is mainly produced by an intermediate N*(1535) we
       * introduce an explicit isospin asymmetry for the production of N*(1535)
       * produced in pn vs. pp similar to [\iref{Teis:1996kx}], eq. 29. */
      if (type_a.is_Nstar1535() || type_b.is_Nstar1535()) {
        return 6.5 * parametrization;
      } else {
        return parametrization;
      }
    }
  } else if (((type_a.is_Deltastar() && type_b.is_nucleon()) ||
              (type_b.is_Deltastar() && type_a.is_nucleon())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → NΔ*
    return 15. / msqr;
  } else if ((type_a.is_Delta() && type_b.is_Delta()) &&
             (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    // NN → ΔΔ
    if (twoI == 2) {
      return 45. / msqr;
    } else if (twoI == 0) {
      return 120. / msqr;
    }
  } else if (((type_a.is_Nstar() && type_b.is_Delta()) ||
              (type_b.is_Nstar() && type_a.is_Delta())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → ΔN*
    return 7. / msqr;
  } else if (((type_a.is_Deltastar() && type_b.is_Delta()) ||
              (type_b.is_Deltastar() && type_a.is_Delta())) &&
             type_a.antiparticle_sign() == type_b.antiparticle_sign()) {
    // NN → ΔΔ*
    if (twoI == 2) {
      return 15. / msqr;
    } else if (twoI == 0) {
      return 25. / msqr;
    }
  } else if ((type_a.is_deuteron() && type_b.pdgcode().is_pion()) ||
             (type_b.is_deuteron() && type_a.pdgcode().is_pion())) {
    /* This parametrization is the result of fitting d+pi->NN cross-section.
     * Already Breit-Wigner-like part provides a good fit, exponential fixes
     * behaviour around the treshold. The d+pi experimental cross-section
     * was taken from Fig. 2 of [\iref{Tanabe:1987vg}]. */
    return 0.055 / (pow_int(sqrts - 2.145, 2) + pow_int(0.065, 2)) *
           (1.0 - std::exp(-(sqrts - 2.0) * 20.0));
  }

  // all cases not listed: zero!
  return 0.;
}

template <class IntegrationMethod>
CollisionBranchList CrossSections::find_nn_xsection_from_type(
    const ParticleTypePtrList& list_res_1,
    const ParticleTypePtrList& list_res_2,
    const IntegrationMethod integrator) const {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  CollisionBranchList channel_list;
  const double s = sqrt_s_ * sqrt_s_;

  // Loop over specified first resonance list
  for (ParticleTypePtr type_res_1 : list_res_1) {
    // Loop over specified second resonance list
    for (ParticleTypePtr type_res_2 : list_res_2) {
      // Check for charge conservation.
      if (type_res_1->charge() + type_res_2->charge() !=
          type_particle_a.charge() + type_particle_b.charge()) {
        continue;
      }

      // loop over total isospin
      for (const int twoI : I_tot_range(type_particle_a, type_particle_b)) {
        const double isospin_factor = isospin_clebsch_gordan_sqr_2to2(
            type_particle_a, type_particle_b, *type_res_1, *type_res_2, twoI);
        // If Clebsch-Gordan coefficient is zero, don't bother with the rest.
        if (std::abs(isospin_factor) < really_small) {
          continue;
        }

        // Integration limits.
        const double lower_limit = type_res_1->min_mass_kinematic();
        const double upper_limit = sqrt_s_ - type_res_2->mass();
        /* Check the available energy (requiring it to be a little above the
         * threshold, because the integration will not work if it's too close).
         */
        if (upper_limit - lower_limit < 1E-3) {
          continue;
        }

        // Calculate matrix element.
        const double matrix_element = nn_to_resonance_matrix_element(
            sqrt_s_, *type_res_1, *type_res_2, twoI);
        if (matrix_element <= 0.) {
          continue;
        }

        /* Calculate resonance production cross section
         * using the Breit-Wigner distribution as probability amplitude.
         * Integrate over the allowed resonance mass range. */
        const double resonance_integral = integrator(*type_res_1, *type_res_2);

        /** Cross section for 2->2 process with 1/2 resonance(s) in final state.
         * Based on Eq. (46) in \iref{Weil:2013mya} and Eq. (3.29) in
         * \iref{Bass:1998ca} */
        const double spin_factor =
            (type_res_1->spin() + 1) * (type_res_2->spin() + 1);
        const double xsection = isospin_factor * spin_factor * matrix_element *
                                resonance_integral / (s * cm_momentum());

        if (xsection > really_small) {
          channel_list.push_back(std::make_unique<CollisionBranch>(
              *type_res_1, *type_res_2, xsection, ProcessType::TwoToTwo));
          logg[LCrossSections].debug(
              "Found 2->2 creation process for resonance ", type_res_1, ", ",
              type_res_2);
          logg[LCrossSections].debug("2->2 with original particles: ",
                                     type_particle_a, type_particle_b);
        }
      }
    }
  }
  return channel_list;
}

double CrossSections::string_probability(
    const ScatterActionsFinderParameters& finder_parameters) const {
  /* string fragmentation is enabled when strings_switch is on and the process
   * is included in pythia. */
  if (!finder_parameters.strings_switch) {
    return 0.;
  }

  const ParticleType& t1 = incoming_particles_[0].type();
  const ParticleType& t2 = incoming_particles_[1].type();
  const bool treat_BBbar_with_strings =
      (finder_parameters.nnbar_treatment == NNbarTreatment::Strings);
  const bool is_NN_scattering =
      t1.is_nucleon() && t2.is_nucleon() &&
      t1.antiparticle_sign() == t2.antiparticle_sign();
  const bool is_BBbar_scattering =
      (treat_BBbar_with_strings && is_BBbar_pair_ &&
       finder_parameters.use_AQM) ||
      (t1.is_nucleon() && t2.is_nucleon() &&
       t1.antiparticle_sign() != t2.antiparticle_sign());
  const bool is_Npi_scattering = (t1.pdgcode().is_pion() && t2.is_nucleon()) ||
                                 (t1.is_nucleon() && t2.pdgcode().is_pion());
  /* True for baryon-baryon, anti-baryon-anti-baryon, baryon-meson,
   * anti-baryon-meson and meson-meson*/
  const bool is_AQM_scattering =
      finder_parameters.use_AQM &&
      ((t1.is_baryon() && t2.is_baryon() &&
        t1.antiparticle_sign() == t2.antiparticle_sign()) ||
       ((t1.is_baryon() && t2.is_meson()) ||
        (t2.is_baryon() && t1.is_meson())) ||
       (t1.is_meson() && t2.is_meson()));
  const double mass_sum =
      incoming_particles_[0].pole_mass() + incoming_particles_[1].pole_mass();

  if (!is_NN_scattering && !is_BBbar_scattering && !is_Npi_scattering &&
      !is_AQM_scattering) {
    return 0.;
  } else if (is_NNbar_pair_ && !treat_BBbar_with_strings) {
    return 0.;
  } else if (is_BBbar_scattering) {
    // BBbar only goes through strings, so there are no "window" considerations
    return 1.;
  } else {
    /* true for K+ p and K0 p (+ antiparticles), which have special treatment
     * to fit data */
    const PdgCode pdg1 = t1.pdgcode(), pdg2 = t2.pdgcode();
    const bool is_KplusP =
        ((pdg1 == pdg::K_p || pdg1 == pdg::K_z) && (pdg2 == pdg::p)) ||
        ((pdg2 == pdg::K_p || pdg2 == pdg::K_z) && (pdg1 == pdg::p)) ||
        ((pdg1 == -pdg::K_p || pdg1 == -pdg::K_z) && (pdg2 == -pdg::p)) ||
        ((pdg2 == -pdg::K_p || pdg2 == -pdg::K_z) && (pdg1 == -pdg::p));
    // where to start the AQM strings above mass sum
    double aqm_offset =
        finder_parameters.transition_high_energy.sqrts_add_lower;
    if (is_KplusP) {
      /* for this specific case we have data. This corresponds to the point
       * where the AQM parametrization is smaller than the current 2to2
       * parametrization, which starts growing and diverges from exp. data */
      aqm_offset = finder_parameters.transition_high_energy.KN_offset;
    } else if (pdg1.is_pion() && pdg2.is_pion()) {
      aqm_offset = finder_parameters.transition_high_energy.pipi_offset;
    }
    /* if we do not use the probability transition algorithm, this is always a
     * string contribution if the energy is large enough */
    if (!finder_parameters.strings_with_probability) {
      return static_cast<double>(sqrt_s_ > mass_sum + aqm_offset);
    }
    /* No strings at low energy, only strings at high energy and
     * a transition region in the middle. Determine transition region: */
    double region_lower, region_upper;
    if (is_Npi_scattering) {
      std::tie(region_lower, region_upper) =
          finder_parameters.transition_high_energy.sqrts_range_Npi;
    } else if (is_NN_scattering) {
      std::tie(region_lower, region_upper) =
          finder_parameters.transition_high_energy.sqrts_range_NN;
    } else {  // AQM - Additive Quark Model
      /* Transition region around 0.9 larger than the sum of pole masses;
       * highly arbitrary, feel free to improve */
      region_lower = mass_sum + aqm_offset;
      region_upper = mass_sum + aqm_offset +
                     finder_parameters.transition_high_energy.sqrts_range_width;
    }

    if (sqrt_s_ > region_upper) {
      return 1.;
    } else if (sqrt_s_ < region_lower) {
      return 0.;
    } else {
      // Rescale transition region to [-1, 1]
      return probability_transit_high(region_lower, region_upper);
    }
  }
}

double CrossSections::probability_transit_high(
    const double region_lower, const double region_upper) const {
  if (sqrt_s_ < region_lower) {
    return 0.;
  }

  if (sqrt_s_ > region_upper) {
    return 1.;
  }

  double x = (sqrt_s_ - 0.5 * (region_lower + region_upper)) /
             (region_upper - region_lower);
  assert(x >= -0.5 && x <= 0.5);
  double prob = 0.5 * (std::sin(M_PI * x) + 1.0);
  assert(prob >= 0. && prob <= 1.);

  return prob;
}

}  // namespace smash
