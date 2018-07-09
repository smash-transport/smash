/*
 *
 *    Copyright (c) 2013-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/crosssections.h"

#include "smash/clebschgordan.h"
#include "smash/constants.h"
#include "smash/kinematics.h"
#include "smash/logging.h"
#include "smash/parametrizations.h"
#include "smash/particletype.h"
#include "smash/pow.h"

namespace smash {

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
      (pcm * a.iso_multiplet()->get_integral_RR(b, sqrts));
  return spin_factor * symmetry_factor * momentum_factor;
}

/**
 * Helper function:
 * Add a 2-to-2 channel to a collision branch list given a cross section.
 *
 * The cross section is only calculated if there is enough energy
 * for the process. If the cross section is small, the branch is not added.
 */
template <typename F>
void add_channel(CollisionBranchList& process_list, F get_xsection,
                 double sqrts, const ParticleType& type_a,
                 const ParticleType& type_b) {
  const double sqrt_s_min =
      type_a.min_mass_spectral() + type_b.min_mass_spectral();
  if (sqrts <= sqrt_s_min) {
    return;
  }
  const auto xsection = get_xsection();
  if (xsection > really_small) {
    process_list.push_back(make_unique<CollisionBranch>(
        type_a, type_b, xsection, ProcessType::TwoToTwo));
  }
}

/**
 * Helper function:
 * Append a list of processes to another (main) list of processes.
 */
static void append_list(CollisionBranchList& main_list,
                        CollisionBranchList in_list) {
  main_list.reserve(main_list.size() + in_list.size());
  for (auto& proc : in_list) {
    main_list.emplace_back(std::move(proc));
  }
}

/**
 * Helper function:
 * Sum all cross sections of the given process list.
 */
static double sum_xs_of(CollisionBranchList& list) {
  double xs_sum = 0.0;
  for (auto& proc : list) {
    xs_sum += proc->weight();
  }
  return xs_sum;
}

CrossSections::CrossSections(const ParticleList& incoming_particles,
                             const double sqrt_s)
    : incoming_particles_(incoming_particles), sqrt_s_(sqrt_s),
      is_BBbar_pair_(incoming_particles_[0].type().is_baryon() &&
                     incoming_particles_[1].type().is_baryon() &&
                     incoming_particles_[0].type().antiparticle_sign() ==
                     -incoming_particles_[1].type().antiparticle_sign()) {
}

CollisionBranchList CrossSections::generate_collision_list(
    double elastic_parameter, bool two_to_one_switch,
    ReactionsBitSet included_2to2, double low_snn_cut, bool strings_switch,
    bool use_AQM, bool strings_with_probability, NNbarTreatment nnbar_treatment,
    StringProcess* string_process) {
  CollisionBranchList process_list;
  const ParticleType& t1 = incoming_particles_[0].type();
  const ParticleType& t2 = incoming_particles_[1].type();

  const bool is_pythia =
      strings_with_probability &&
      decide_string(strings_switch, strings_with_probability, use_AQM,
                    nnbar_treatment == NNbarTreatment::Strings);

  /* Elastic collisions between two nucleons with sqrt_s below
   * low_snn_cut can not happen. */
  const bool reject_by_nucleon_elastic_cutoff =
      t1.is_nucleon() && t2.is_nucleon() &&
      t1.antiparticle_sign() == t2.antiparticle_sign() && sqrt_s_ < low_snn_cut;
  bool incl_elastic = included_2to2[IncludedReactions::Elastic];
  if (incl_elastic && !reject_by_nucleon_elastic_cutoff) {
    process_list.emplace_back(elastic(elastic_parameter, use_AQM));
  }
  if (is_pythia) {
    /* String-excitation cross section =
     * Parametrized total cross - the contributions
     * from all other present channels. */
    const double sig_string =
        std::max(0., high_energy() - elastic_parametrization(use_AQM));
    append_list(process_list,
        string_excitation(sig_string, string_process, use_AQM));
  } else {
    if (two_to_one_switch) {
      // resonance formation (2->1)
      append_list(process_list, two_to_one());
    }
    if (included_2to2.any()) {
      // 2->2 (inelastic)
      append_list(process_list, two_to_two(included_2to2));
    }
  }
  /* NNbar annihilation thru NNbar → ρh₁(1170); combined with the decays
   * ρ → ππ and h₁(1170) → πρ, this gives a final state of 5 pions.
   * Only use in cases when detailed balance MUST happen, i.e. in a box! */
  if (nnbar_treatment == NNbarTreatment::Resonances) {
    if (t1.is_nucleon() && t2.pdgcode() == t1.get_antiparticle()->pdgcode()) {
      /* Has to be called after the other processes are already determined,
       * so that the sum of the cross sections includes all other processes. */
      process_list.emplace_back(NNbar_annihilation(sum_xs_of(process_list)));
    }
    if ((t1.pdgcode() == pdg::rho_z && t2.pdgcode() == pdg::h1) ||
        (t1.pdgcode() == pdg::h1 && t2.pdgcode() == pdg::rho_z)) {
      append_list(process_list, NNbar_creation());
    }
  }
  return process_list;
}

CollisionBranchPtr CrossSections::elastic(double elast_par, bool use_AQM) {
  double elastic_xs = 0.;
  if (elast_par >= 0.) {
    // use constant elastic cross section from config file
    elastic_xs = elast_par;
  } else {
    // use parametrization
    elastic_xs = elastic_parametrization(use_AQM);
  }
  return make_unique<CollisionBranch>(incoming_particles_[0].type(),
                                      incoming_particles_[1].type(), elastic_xs,
                                      ProcessType::Elastic);
}

double CrossSections::elastic_parametrization(bool use_AQM) {
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
  } else if (pdg_a.is_nucleus() || pdg_b.is_nucleus()) {
    const PdgCode& pdg_nucleus = pdg_a.is_nucleus() ? pdg_a : pdg_b;
    const PdgCode& pdg_other = pdg_a.is_nucleus() ? pdg_b : pdg_a;
    const bool is_deuteron =
        std::abs(pdg_nucleus.get_decimal()) == pdg::decimal_d;
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
      // meson-meson goes through scaling from π+p parametrization
      elastic_xs = 2. / 3. * piplusp_elastic_high_energy(s, m1, m2);
    }
    elastic_xs *=
        (1. - 0.4 * pdg_a.frac_strange()) * (1. - 0.4 * pdg_b.frac_strange());
  }
  return elastic_xs;
}

double CrossSections::nn_el() {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  const double s = sqrt_s_ * sqrt_s_;

  // Use parametrized cross sections.
  double sig_el = -1.;
  if (pdg_a.is_antiparticle_of(pdg_b)) {
    sig_el = ppbar_elastic(s);
  } else if (pdg_a.is_nucleon() && pdg_b.is_nucleon()) {
    sig_el = (pdg_a == pdg_b) ? pp_elastic(s)
                              : np_elastic(s);
  } else {
    // AQM - Additive Quark Model
    const double m1 = incoming_particles_[0].effective_mass();
    const double m2 = incoming_particles_[1].effective_mass();
    sig_el = pp_elastic_high_energy(s, m1, m2);
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

double CrossSections::npi_el() {
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

double CrossSections::nk_el() {
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

CollisionBranchList CrossSections::two_to_one() {
  const auto& log = logger<LogArea::CrossSections>();
  CollisionBranchList resonance_process_list;
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  const double p_cm_sqr = pCM_sqr(sqrt_s_, m1, m2);

  // Find all the possible resonances
  for (const ParticleType& type_resonance : ParticleType::list_all()) {
    /* Not a resonance, go to next type of particle */
    if (type_resonance.is_stable()) {
      continue;
    }

    // Same resonance as in the beginning, ignore
    if ((!type_particle_a.is_stable() &&
         type_resonance.pdgcode() == type_particle_a.pdgcode()) ||
        (!type_particle_b.is_stable() &&
         type_resonance.pdgcode() == type_particle_b.pdgcode())) {
      continue;
    }

    double resonance_xsection = formation(type_resonance, p_cm_sqr);

    // If cross section is non-negligible, add resonance to the list
    if (resonance_xsection > really_small) {
      resonance_process_list.push_back(make_unique<CollisionBranch>(
          type_resonance, resonance_xsection, ProcessType::TwoToOne));
      log.debug("Found resonance: ", type_resonance);
      log.debug(type_particle_a.name(), type_particle_b.name(), "->",
                type_resonance.name(), " at sqrt(s)[GeV] = ", sqrt_s_,
                " with xs[mb] = ", resonance_xsection);
    }
  }
  return resonance_process_list;
}

double CrossSections::formation(const ParticleType& type_resonance,
                                double cm_momentum_sqr) {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();
  // Check for charge conservation.
  if (type_resonance.charge() !=
      type_particle_a.charge() + type_particle_b.charge()) {
    return 0.;
  }

  // Check for baryon-number conservation.
  if (type_resonance.baryon_number() !=
      type_particle_a.baryon_number() + type_particle_b.baryon_number()) {
    return 0.;
  }

  // Calculate partial in-width.
  const double partial_width = type_resonance.get_partial_in_width(
      sqrt_s_, incoming_particles_[0], incoming_particles_[1]);
  if (partial_width <= 0.) {
    return 0.;
  }

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

CollisionBranchList CrossSections::two_to_two(ReactionsBitSet included_2to2) {
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
      process_list = nk_xx(included_2to2);
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

CollisionBranchList CrossSections::bb_xx_except_nn(
    ReactionsBitSet included_2to2) {
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

CollisionBranchList CrossSections::nn_xx(ReactionsBitSet included_2to2) {
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
          return type_res_1.iso_multiplet()->get_integral_RR(type_res_2, sqrts);
        });
    process_list.reserve(process_list.size() + channel_list.size());
    std::move(channel_list.begin(), channel_list.end(),
              std::inserter(process_list, process_list.end()));
    channel_list.clear();
  }

  // Find N N → dπ and N̅ N̅→ d̅π channels.
  ParticleTypePtr deutron =
      ParticleType::try_find(PdgCode::from_decimal(pdg::decimal_d));
  ParticleTypePtr antideutron =
      ParticleType::try_find(PdgCode::from_decimal(pdg::decimal_antid));
  ParticleTypePtr pim = ParticleType::try_find(pdg::pi_m);
  ParticleTypePtr pi0 = ParticleType::try_find(pdg::pi_z);
  ParticleTypePtr pip = ParticleType::try_find(pdg::pi_p);
  // Make sure all the necessary particle types are found
  if (deutron && antideutron && pim && pi0 && pip) {
    const ParticleTypePtrList deutron_list = {deutron};
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

CollisionBranchList CrossSections::nk_xx(ReactionsBitSet included_2to2) {
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

  bool incl_KN_to_KN = included_2to2[IncludedReactions::KN_to_KN] == 1;
  bool incl_KN_to_KDelta = included_2to2[IncludedReactions::KN_to_KDelta] == 1;
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
            add_channel(process_list,
                        [&] { return kminusp_piminussigmaplus(sqrt_s_); },
                        sqrt_s_, type_pi_m, type_Sigma_p);
            add_channel(process_list,
                        [&] { return kminusp_piplussigmaminus(sqrt_s_); },
                        sqrt_s_, type_pi_p, type_Sigma_m);
            add_channel(process_list,
                        [&] { return kminusp_pi0sigma0(sqrt_s_); }, sqrt_s_,
                        type_pi_z, type_Sigma_z);
            add_channel(process_list,
                        [&] { return kminusp_pi0lambda(sqrt_s_); }, sqrt_s_,
                        type_pi_z, type_Lambda);
          }
          if (incl_KN_to_KN) {
            const auto& type_n = ParticleType::find(pdg::n);
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            add_channel(process_list, [&] { return kminusp_kbar0n(s); },
                        sqrt_s_, type_Kbar_z, type_n);
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
            add_channel(process_list,
                        [&] { return kminusn_piminussigma0(sqrt_s_); }, sqrt_s_,
                        type_pi_m, type_Sigma_z);
            add_channel(process_list,
                        [&] { return kminusn_pi0sigmaminus(sqrt_s_); }, sqrt_s_,
                        type_pi_z, type_Sigma_m);
            add_channel(process_list,
                        [&] { return kminusn_piminuslambda(sqrt_s_); }, sqrt_s_,
                        type_pi_m, type_Lambda);
          }
          break;
        }
        case -pdg::p: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            const auto& type_Delta_pp_bar = ParticleType::find(-pdg::Delta_pp);
            const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_Kbar_z,
                                                    type_Delta_pp_bar);
                        },
                        sqrt_s_, type_Kbar_z, type_Delta_pp_bar);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_m, type_Delta_p_bar);
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
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_Kbar_z,
                                                    type_Delta_p_bar);
                        },
                        sqrt_s_, type_Kbar_z, type_Delta_p_bar);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_m, type_Delta_z_bar);
                        },
                        sqrt_s_, type_K_m, type_Delta_z_bar);
          }
          if (incl_KN_to_KN) {
            const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
            const auto& type_p_bar = ParticleType::find(-pdg::p);
            add_channel(process_list, [&] { return kplusn_k0p(s); }, sqrt_s_,
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
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_z, type_Delta_pp);
                        },
                        sqrt_s_, type_K_z, type_Delta_pp);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
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
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_z, type_Delta_p);
                        },
                        sqrt_s_, type_K_z, type_Delta_p);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_p, type_Delta_z);
                        },
                        sqrt_s_, type_K_p, type_Delta_z);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            const auto& type_p = ParticleType::find(pdg::p);
            add_channel(process_list, [&] { return kplusn_k0p(s); }, sqrt_s_,
                        type_K_z, type_p);
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
            add_channel(process_list,
                        [&] { return kminusp_piminussigmaplus(sqrt_s_); },
                        sqrt_s_, type_pi_p, type_Sigma_p_bar);
            add_channel(process_list,
                        [&] { return kminusp_piplussigmaminus(sqrt_s_); },
                        sqrt_s_, type_pi_m, type_Sigma_m_bar);
            add_channel(process_list,
                        [&] { return kminusp_pi0sigma0(sqrt_s_); }, sqrt_s_,
                        type_pi_z, type_Sigma_z_bar);
            add_channel(process_list,
                        [&] { return kminusp_pi0lambda(sqrt_s_); }, sqrt_s_,
                        type_pi_z, type_Lambda_bar);
          }
          if (incl_KN_to_KN) {
            const auto& type_n_bar = ParticleType::find(-pdg::n);
            const auto& type_K_z = ParticleType::find(pdg::K_z);
            add_channel(process_list, [&] { return kminusp_kbar0n(s); },
                        sqrt_s_, type_K_z, type_n_bar);
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
            add_channel(process_list,
                        [&] { return kminusn_piminussigma0(sqrt_s_); }, sqrt_s_,
                        type_pi_p, type_Sigma_z_bar);
            add_channel(process_list,
                        [&] { return kminusn_pi0sigmaminus(sqrt_s_); }, sqrt_s_,
                        type_pi_z, type_Sigma_m_bar);
            add_channel(process_list,
                        [&] { return kminusn_piminuslambda(sqrt_s_); }, sqrt_s_,
                        type_pi_p, type_Lambda_bar);
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
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_z, type_Delta_p);
                        },
                        sqrt_s_, type_K_z, type_Delta_p);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_p, type_Delta_z);
                        },
                        sqrt_s_, type_K_p, type_Delta_z);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_n = ParticleType::find(pdg::n);
            add_channel(process_list,
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
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_z, type_Delta_z);
                        },
                        sqrt_s_, type_K_z, type_Delta_z);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_p, type_Delta_m);
                        },
                        sqrt_s_, type_K_p, type_Delta_m);
          }
          break;
        }
        case -pdg::n: {
          if (incl_KN_to_KN) {
            const auto& type_K_p = ParticleType::find(pdg::K_p);
            const auto& type_p_bar = ParticleType::find(-pdg::p);
            add_channel(process_list, [&] { return kminusp_kbar0n(s); },
                        sqrt_s_, type_K_p, type_p_bar);
          }
          break;
        }
      }
      break;
    }
    case pdg::Kbar_z:
      switch (pdg_nucleon) {
        case pdg::n: {
          if (incl_KN_to_KN) {
            const auto& type_p = ParticleType::find(pdg::p);
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            add_channel(process_list, [&] { return kminusp_kbar0n(s); },
                        sqrt_s_, type_K_m, type_p);
          }
          break;
        }
        case -pdg::p: {
          if (incl_KN_to_KDelta) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_Kbar_z = type_kaon;
            const auto& type_Delta_bar_m = ParticleType::find(-pdg::Delta_p);
            const auto& type_Delta_bar_z = ParticleType::find(-pdg::Delta_z);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_Kbar_z,
                                                    type_Delta_bar_m);
                        },
                        sqrt_s_, type_Kbar_z, type_Delta_bar_m);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusn * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_m, type_Delta_bar_z);
                        },
                        sqrt_s_, type_K_m, type_Delta_bar_z);
          }
          if (incl_KN_to_KN) {
            const auto& type_K_m = ParticleType::find(pdg::K_m);
            const auto& type_n_bar = ParticleType::find(-pdg::n);
            add_channel(process_list,
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
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_Kbar_z,
                                                    type_Delta_z_bar);
                        },
                        sqrt_s_, type_Kbar_z, type_Delta_z_bar);
            add_channel(process_list,
                        [&] {
                          return sigma_kplusp * kaon_nucleon_ratios.get_ratio(
                                                    type_nucleon, type_kaon,
                                                    type_K_m, type_Delta_m_bar);
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

CollisionBranchList CrossSections::deltak_xx(ReactionsBitSet included_2to2) {
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
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p,
                                                      type_K_p) *
                           kaon_nucleon_ratios.get_ratio(
                               type_p, type_K_p, type_kaon, type_delta) *
                           kplusp_inelastic_background(s);
                  },
                  sqrt_s_, type_p, type_K_p);
      break;
    }
    case pack(-pdg::Delta_pp, pdg::Kbar_z):
    case pack(-pdg::Delta_p, pdg::K_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p_bar,
                                                      type_K_m) *
                           kaon_nucleon_ratios.get_ratio(
                               type_p_bar, type_K_m, type_kaon, type_delta) *
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
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_n,
                                                      type_K_p) *
                           kaon_nucleon_ratios.get_ratio(
                               type_n, type_K_p, type_kaon, type_delta) *
                           kplusn_inelastic_background(s);
                  },
                  sqrt_s_, type_n, type_K_p);

      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p,
                                                      type_K_z) *
                           kaon_nucleon_ratios.get_ratio(
                               type_p, type_K_z, type_kaon, type_delta) *
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
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_n_bar,
                                                      type_K_m) *
                           kaon_nucleon_ratios.get_ratio(
                               type_n_bar, type_K_m, type_kaon, type_delta) *
                           kplusn_inelastic_background(s);
                  },
                  sqrt_s_, type_n_bar, type_K_m);

      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p_bar,
                                                      type_Kbar_z) *
                           kaon_nucleon_ratios.get_ratio(
                               type_p_bar, type_Kbar_z, type_kaon, type_delta) *
                           kplusn_inelastic_background(s);
                  },
                  sqrt_s_, type_p_bar, type_Kbar_z);
      break;
    }
    case pack(pdg::Delta_z, pdg::K_z):
    case pack(pdg::Delta_m, pdg::K_p): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_z = ParticleType::find(pdg::K_z);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_n,
                                                      type_K_z) *
                           kaon_nucleon_ratios.get_ratio(
                               type_n, type_K_z, type_kaon, type_delta) *
                           kplusp_inelastic_background(s);
                  },
                  sqrt_s_, type_n, type_K_z);
      break;
    }
    case pack(-pdg::Delta_z, pdg::Kbar_z):
    case pack(-pdg::Delta_m, pdg::K_m): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_n_bar,
                                                      type_Kbar_z) *
                           kaon_nucleon_ratios.get_ratio(
                               type_n_bar, type_Kbar_z, type_kaon, type_delta) *
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

CollisionBranchList CrossSections::ypi_xx(ReactionsBitSet included_2to2) {
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
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_n, type_K_m) *
                           kminusn_piminussigma0(sqrt_s_);
                  },
                  sqrt_s_, type_n, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_p): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_n_bar,
                                                          type_K_p) *
                           kminusn_piminussigma0(sqrt_s_);
                  },
                  sqrt_s_, type_n_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_z): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_n, type_K_m) *
                           kminusn_pi0sigmaminus(sqrt_s_);
                  },
                  sqrt_s_, type_n, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_m, pdg::pi_z): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_n_bar,
                                                          type_K_p) *
                           kminusn_pi0sigmaminus(sqrt_s_);
                  },
                  sqrt_s_, type_n_bar, type_K_p);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_m): {
      const auto& type_n = ParticleType::find(pdg::n);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_n, type_K_m) *
                           kminusn_piminuslambda(sqrt_s_);
                  },
                  sqrt_s_, type_n, type_K_m);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_p): {
      const auto& type_n_bar = ParticleType::find(-pdg::n);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_n_bar,
                                                          type_K_p) *
                           kminusn_piminuslambda(sqrt_s_);
                  },
                  sqrt_s_, type_n_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_z, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_p, type_K_m) *
                           kminusp_pi0sigma0(sqrt_s_);
                  },
                  sqrt_s_, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_z, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_p_bar,
                                                          type_K_p) *
                           kminusp_pi0sigma0(sqrt_s_);
                  },
                  sqrt_s_, type_p_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_m, pdg::pi_p): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_p, type_K_m) *
                           kminusp_piplussigmaminus(sqrt_s_);
                  },
                  sqrt_s_, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_m, pdg::pi_m): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_p_bar,
                                                          type_K_p) *
                           kminusp_piplussigmaminus(sqrt_s_);
                  },
                  sqrt_s_, type_p_bar, type_K_p);
      break;
    }
    case pack(pdg::Lambda, pdg::pi_z): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_p, type_K_m) *
                           kminusp_pi0lambda(sqrt_s_);
                  },
                  sqrt_s_, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Lambda, pdg::pi_z): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_p_bar,
                                                          type_K_p) *
                           kminusp_pi0lambda(sqrt_s_);
                  },
                  sqrt_s_, type_p_bar, type_K_p);
      break;
    }
    case pack(pdg::Sigma_p, pdg::pi_m): {
      const auto& type_p = ParticleType::find(pdg::p);
      const auto& type_K_m = ParticleType::find(pdg::K_m);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(
                               s, type_hyperon, type_pion, type_p, type_K_m) *
                           kminusp_piminussigmaplus(sqrt_s_);
                  },
                  sqrt_s_, type_p, type_K_m);
      break;
    }
    case pack(-pdg::Sigma_p, pdg::pi_p): {
      const auto& type_p_bar = ParticleType::find(-pdg::p);
      const auto& type_K_p = ParticleType::find(pdg::K_p);
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_stable(s, type_hyperon,
                                                          type_pion, type_p_bar,
                                                          type_K_p) *
                           kminusp_piminussigmaplus(sqrt_s_);
                  },
                  sqrt_s_, type_p_bar, type_K_p);
      break;
    }
    default:
      break;
  }

  return process_list;
}

CollisionBranchList CrossSections::dpi_xx(ReactionsBitSet
                                          /*included_2to2*/) {
  const auto& log = logger<LogArea::ScatterAction>();
  CollisionBranchList process_list;
  const double sqrts = sqrt_s_;
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();

  // pi d -> N N
  if ((type_a.is_deuteron() && type_b.pdgcode().is_pion()) ||
      (type_b.is_deuteron() && type_a.pdgcode().is_pion())) {
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
            process_list.push_back(make_unique<CollisionBranch>(
                *nuc_a, *nuc_b, xsection, ProcessType::TwoToTwo));
            log.debug(type_a.name(), type_b.name(), "->", nuc_a->name(),
                      nuc_b->name(), " at sqrts [GeV] = ", sqrts,
                      " with cs[mb] = ", xsection);
          }
        }
      }
    }
  }

  // pi d -> pi d' (effectively pi d -> pi p n)  AND reverse, pi d' -> pi d
  if (((type_a.is_deuteron() || type_a.is_dprime()) &&
       type_b.pdgcode().is_pion()) ||
      ((type_b.is_deuteron() || type_b.is_dprime()) &&
       type_a.pdgcode().is_pion())) {
    const ParticleType& type_pi = type_a.pdgcode().is_pion() ? type_a : type_b;
    const ParticleType& type_nucleus = type_a.is_nucleus() ? type_a : type_b;
    ParticleTypePtrList nuclei = ParticleType::list_light_nuclei();
    const double s = sqrt_s_ * sqrt_s_;
    for (ParticleTypePtr produced_nucleus : nuclei) {
      // Elastic collisions are treated in a different function
      if (produced_nucleus == &type_nucleus ||
          produced_nucleus->charge() != type_nucleus.charge() ||
          produced_nucleus->baryon_number() != type_nucleus.baryon_number()) {
        continue;
      }
      // same matrix element for πd and πd̅
      const double tmp =
          sqrts - type_a.min_mass_kinematic() - type_b.min_mass_kinematic();
      /* Matrix element is fit to match the inelastic pi+ d -> pi+ n p
       * cross-section from the Fig. 5 of [\iref{Arndt:1994bs}]. */
      const double matrix_element =
          295.5 + 2.862 / (0.00283735 + pow_int(sqrts - 2.181, 2)) +
          0.0672 / pow_int(tmp, 2) - 6.61753 / tmp;
      const double spin_factor =
          (produced_nucleus->spin() + 1) * (type_pi.spin() + 1);
      /* Isospin factor is always the same, so it is included into the
       * matrix element.
       * Symmetry factor is always 1 here.
       * The (hbarc)^2/16 pi factor is absorbed into matrix element. */
      double xsection = matrix_element * spin_factor / (s * cm_momentum());
      if (produced_nucleus->is_stable()) {
        assert(!type_nucleus.is_stable());
        xsection *= pCM_from_s(s, type_pi.mass(), produced_nucleus->mass());
      } else {
        assert(type_nucleus.is_stable());
        const double resonance_integral =
            produced_nucleus->iso_multiplet()->get_integral_piR(sqrts);
        xsection *= resonance_integral;
        log.debug("Resonance integral ", resonance_integral,
                  ", matrix element: ", matrix_element,
                  ", cm_momentum: ", cm_momentum());
      }
      process_list.push_back(make_unique<CollisionBranch>(
          type_pi, *produced_nucleus, xsection, ProcessType::TwoToTwo));
      log.debug(type_pi.name(), type_nucleus.name(), "→ ", type_pi.name(),
                produced_nucleus->name(), " at ", sqrts,
                " GeV, xs[mb] = ", xsection);
    }
  }
  return process_list;
}

CollisionBranchList CrossSections::dn_xx(ReactionsBitSet /*included_2to2*/) {
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();
  const ParticleType& type_N = type_a.is_nucleon() ? type_a : type_b;
  const ParticleType& type_nucleus = type_a.is_nucleus() ? type_a : type_b;
  CollisionBranchList process_list;
  ParticleTypePtrList nuclei = ParticleType::list_light_nuclei();
  const double s = sqrt_s_ * sqrt_s_;
  const double sqrts = sqrt_s_;

  for (ParticleTypePtr produced_nucleus : nuclei) {
    // No elastic collisions for now, respect conservation laws
    if (produced_nucleus == &type_nucleus ||
        produced_nucleus->charge() != type_nucleus.charge() ||
        produced_nucleus->baryon_number() != type_nucleus.baryon_number()) {
      continue;
    }
    double matrix_element = 0.0;
    if (std::signbit(type_N.baryon_number()) ==
        std::signbit(type_nucleus.baryon_number())) {
      // Nd → Nd', N̅d̅→ N̅d̅' and reverse
      const double tmp = sqrts - type_N.min_mass_kinematic() -
                         type_nucleus.min_mass_kinematic();
      assert(tmp >= 0.0);
      /* Fit to match experimental cross-section Nd -> Nnp from
       * [\iref{Carlson1973}] */
      matrix_element = 79.0474 / std::pow(tmp, 0.7897) + 654.596 * tmp;
    } else {
      /* N̅d →  N̅d', Nd̅→ Nd̅' and reverse
       * Fit to roughly match experimental cross-section N̅d -> N̅ np from
       * [\iref{Bizzarri:1973sp}]. */
      matrix_element = 681.4;
    }
    const double spin_factor =
        (produced_nucleus->spin() + 1) * (type_N.spin() + 1);
    /* Isospin factor is always the same, so it is included into matrix element
     * Symmetry factor is always 1 here
     * Absorb (hbarc)^2/16 pi factor into matrix element */
    double xsection = matrix_element * spin_factor / (s * cm_momentum());
    if (produced_nucleus->is_stable()) {
      assert(!type_nucleus.is_stable());
      xsection *= pCM_from_s(s, type_N.mass(), produced_nucleus->mass());
    } else {
      assert(type_nucleus.is_stable());
      const double resonance_integral =
          produced_nucleus->iso_multiplet()->get_integral_NR(sqrts);
      xsection *= resonance_integral;
    }
    process_list.push_back(make_unique<CollisionBranch>(
        type_N, *produced_nucleus, xsection, ProcessType::TwoToTwo));
    const auto& log = logger<LogArea::ScatterAction>();
    log.debug(type_N.name(), type_nucleus.name(), "→ ", type_N.name(),
              produced_nucleus->name(), " at ", sqrts,
              " GeV, xs[mb] = ", xsection);
  }
  return process_list;
}

CollisionBranchList CrossSections::string_excitation(
    double total_string_xs, StringProcess* string_process, bool use_AQM) {
  const auto& log = logger<LogArea::CrossSections>();

  if (!string_process) {
    throw std::runtime_error("string_process should be initialized.");
  }

  CollisionBranchList channel_list;
  if (total_string_xs <= 0.) {
    return channel_list;
  }

  /* Get PDG id for evaluation of the parametrized cross sections
   * for diffractive processes.
   * (anti-)proton is used for (anti-)baryons and
   * pion is used for mesons.
   * This must be rescaled according to additive quark model
   * in the case of exotic hadrons.
   * Also calculate the multiplicative factor for AQM
   * based on the quark contents. */
  std::array<int, 2> pdgid;
  double AQM_factor = 1.;
  for (int i = 0; i < 2; i++) {
    PdgCode pdg = incoming_particles_[i].type().pdgcode();
    // FIXME: What about charge?
    if (pdg.baryon_number() == 1) {
      pdgid[i] = pdg::p_decimal;
    } else if (pdg.baryon_number() == -1) {
      pdgid[i] = -pdg::p_decimal;
    } else {
      pdgid[i] = pdg::pi_p_decimal;
    }
    AQM_factor *= (1. - 0.4 * pdg.frac_strange());
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
  std::array<double, 3> xs =
      string_process->cross_sections_diffractive(pdgid[0], pdgid[1], sqrt_s_);
  if (use_AQM) {
    for (int ip = 0; ip < 3; ip++) {
      xs[ip] *= AQM_factor;
    }
  }
  double single_diffr_AX = xs[0], single_diffr_XB = xs[1],
         double_diffr = xs[2];
  double single_diffr = single_diffr_AX + single_diffr_XB;
  double diffractive = single_diffr + double_diffr;

  /* The case for baryon/anti-baryon annihilation is treated separately,
   * as in this case we use only one way to break up the particles, namely
   * into 2 mesonic strings of equal mass after annihilating one quark-
   * anti-quark pair. See StringProcess::next_BBbarAnn() */
  double sig_annihilation = 0.0;
  if (can_annihilate) {
    /* In the case of baryon-antibaryon pair,
     * the parametrized cross section for annihilation will be added.
     * See xs_ppbar_annihilation(). */
    double xs_param = xs_ppbar_annihilation(sqrt_s_ * sqrt_s_);
    if (use_AQM) {
      xs_param *= AQM_factor;
    }
    sig_annihilation = std::min(total_string_xs, xs_param);
  }

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
    const double hard_xsec = string_hard_cross_section();
    nondiffractive_soft =
        nondiffractive_all * std::exp(-hard_xsec / nondiffractive_all);
    nondiffractive_hard = nondiffractive_all - nondiffractive_soft;
  }
  log.debug("String cross sections [mb] are");
  log.debug("Single-diffractive AB->AX: ", single_diffr_AX);
  log.debug("Single-diffractive AB->XB: ", single_diffr_XB);
  log.debug("Double-diffractive AB->XX: ", double_diffr);
  log.debug("Soft non-diffractive: ", nondiffractive_soft);
  log.debug("Hard non-diffractive: ", nondiffractive_hard);
  log.debug("B-Bbar annihilation: ", sig_annihilation);

  /* cross section of soft string excitation including annihilation */
  const double sig_string_soft = total_string_xs - nondiffractive_hard;

  /* fill the list of process channels */
  if (sig_string_soft > 0.) {
    channel_list.push_back(make_unique<CollisionBranch>(
        single_diffr_AX,
        ProcessType::StringSoftSingleDiffractiveAX));
    channel_list.push_back(make_unique<CollisionBranch>(
        single_diffr_XB,
        ProcessType::StringSoftSingleDiffractiveXB));
    channel_list.push_back(make_unique<CollisionBranch>(
        double_diffr, ProcessType::StringSoftDoubleDiffractive));
    channel_list.push_back(make_unique<CollisionBranch>(
        nondiffractive_soft, ProcessType::StringSoftNonDiffractive));
    if (can_annihilate) {
      channel_list.push_back(make_unique<CollisionBranch>(
          sig_annihilation, ProcessType::StringSoftAnnihilation));
    }
  }
  if (nondiffractive_hard > 0.) {
    channel_list.push_back(make_unique<CollisionBranch>(
        nondiffractive_hard, ProcessType::StringHard));
  }
  return channel_list;
}

double CrossSections::high_energy() const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  // Disable AQM cross section for KN, because it destroys the total cross
  // section. This is only a temporary fix, see #6262.
  if (((pdg_a == pdg::K_p || pdg_a == pdg::K_z) && pdg_b.is_nucleon()) ||
      ((pdg_b == pdg::K_p || pdg_b == pdg::K_z) && pdg_a.is_nucleon()) ||
      ((pdg_a == -pdg::K_p || pdg_a == -pdg::K_z) &&
       pdg_b.get_antiparticle().is_nucleon()) ||
      ((pdg_b == -pdg::K_p || pdg_b == -pdg::K_z) &&
       pdg_a.get_antiparticle().is_nucleon())) {
    return 0.;
  }

  const double s = sqrt_s_ * sqrt_s_;
  double xs = 0.;

  // Currently all BB collisions use the nucleon-nucleon parametrizations.
  if (pdg_a.is_baryon() && pdg_b.is_baryon()) {
    if (pdg_a == pdg_b) {
      xs = pp_high_energy(s);  // pp, nn
    } else if (pdg_a.is_antiparticle_of(pdg_b)) {
      xs = ppbar_high_energy(s);  // ppbar, nnbar
    } else if (pdg_a.antiparticle_sign() * pdg_b.antiparticle_sign() == 1) {
      xs = np_high_energy(s);  // np, nbarpbar
    } else {
      xs = npbar_high_energy(s);  // npbar, nbarp
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
  xs *= (1. - 0.4 * pdg_a.frac_strange()) * (1. - 0.4 * pdg_b.frac_strange());

  return xs;
}

double CrossSections::string_hard_cross_section() const {
  double cross_sec = 0.;
  const ParticleData& data_a = incoming_particles_[0];
  const ParticleData& data_b = incoming_particles_[1];
  if (data_a.is_baryon() && data_b.is_baryon()) {
    /* Currently nucleon-nucleon cross section is used for all baryon-baryon
     * casees. This will be changed later by applying additive quark model.
     */
    cross_sec = NN_string_hard(sqrt_s_ * sqrt_s_);
  } else if (data_a.is_baryon() || data_b.is_baryon()) {
    /* Currently nucleon-pion cross section is used for all baryon-meson cases.
     * This will be changed later by applying additive quark model.
     */
    cross_sec = Npi_string_hard(sqrt_s_ * sqrt_s_);
  } else {
    /* Currently pion-pion cross section is used for all meson-meson cases.
     * This will be changed later by applying additive quark model.
     */
    cross_sec = pipi_string_hard(sqrt_s_ * sqrt_s_);
  }
  return cross_sec;
}

CollisionBranchPtr CrossSections::NNbar_annihilation(const double current_xs) {
  const auto& log = logger<LogArea::CrossSections>();
  /* Calculate NNbar cross section:
   * Parametrized total minus all other present channels.*/
  const double s = sqrt_s_ * sqrt_s_;
  double nnbar_xsec = std::max(0., ppbar_total(s) - current_xs);
  log.debug("NNbar cross section is: ", nnbar_xsec);
  // Make collision channel NNbar -> ρh₁(1170); eventually decays into 5π
  return make_unique<CollisionBranch>(ParticleType::find(pdg::h1),
                                      ParticleType::find(pdg::rho_z),
                                      nnbar_xsec, ProcessType::TwoToTwo);
}

CollisionBranchList CrossSections::NNbar_creation() {
  const auto& log = logger<LogArea::CrossSections>();
  CollisionBranchList channel_list;
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

  double xsection = detailed_balance_factor_RR(
                        sqrt_s_, pcm, incoming_particles_[0].type(),
                        incoming_particles_[1].type(), type_N, type_Nbar) *
                    std::max(0., ppbar_total(s) - ppbar_elastic(s));
  log.debug("NNbar reverse cross section is: ", xsection);
  channel_list.push_back(make_unique<CollisionBranch>(
      type_N, type_Nbar, xsection, ProcessType::TwoToTwo));
  channel_list.push_back(make_unique<CollisionBranch>(
      ParticleType::find(pdg::n), ParticleType::find(-pdg::n), xsection,
      ProcessType::TwoToTwo));
  return channel_list;
}

CollisionBranchList CrossSections::bar_bar_to_nuc_nuc(
    const bool is_anti_particles) {
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
          process_list.push_back(make_unique<CollisionBranch>(
              *nuc_a, *nuc_b, xsection, ProcessType::TwoToTwo));
          const auto& log = logger<LogArea::CrossSections>();
          log.debug("2->2 absorption with original particles: ", type_a,
                    type_b);
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
   * the collision will be neglected.*/
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
      return 7. / msqr;
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
    const ParticleTypePtrList& list_res_2, const IntegrationMethod integrator) {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const auto& log = logger<LogArea::CrossSections>();
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
          channel_list.push_back(make_unique<CollisionBranch>(
              *type_res_1, *type_res_2, xsection, ProcessType::TwoToTwo));
          log.debug("Found 2->2 creation process for resonance ", type_res_1,
                    ", ", type_res_2);
          log.debug("2->2 with original particles: ", type_particle_a,
                    type_particle_b);
        }
      }
    }
  }
  return channel_list;
}

bool CrossSections::decide_string(bool strings_switch,
                                  bool use_transition_probability, bool use_AQM,
                                  bool treat_BBbar_with_strings) const {
  /* string fragmentation is enabled when strings_switch is on and the process
   * is included in pythia. */
  if (!strings_switch) {
    return false;
  }

  const ParticleType& t1 = incoming_particles_[0].type();
  const ParticleType& t2 = incoming_particles_[1].type();

  const bool is_NN_scattering =
      t1.is_nucleon() && t2.is_nucleon() &&
      t1.antiparticle_sign() == t2.antiparticle_sign();
  const bool is_BBbar_scattering =
      treat_BBbar_with_strings && is_BBbar_pair_;
  const bool is_Npi_scattering = (t1.pdgcode().is_pion() && t2.is_nucleon()) ||
                                 (t1.is_nucleon() && t2.pdgcode().is_pion());
  /* True for baryon-baryon, anti-baryon-anti-baryon, baryon-meson,
   * anti-baryon-meson and meson-meson*/
  const bool is_AQM_scattering =
      use_AQM && ((t1.is_baryon() && t2.is_baryon() &&
                   t1.antiparticle_sign() == t2.antiparticle_sign()) ||
                  ((t1.is_baryon() && t2.is_meson()) ||
                   (t2.is_baryon() && t1.is_meson())) ||
                  (t1.is_meson() && t2.is_meson()));
  const double mass_sum =
      incoming_particles_[0].pole_mass() + incoming_particles_[1].pole_mass();

  if (!is_NN_scattering && !is_BBbar_scattering && !is_Npi_scattering &&
      !is_AQM_scattering) {
    return false;
  } else if (is_BBbar_scattering) {
    // BBbar only goes through strings, so there are no "window" considerations
    return true;
  } else {
    /* if we do not use the probability transition algorithm, this is always a
     * string contribution if the energy is large enough */
    if (!use_transition_probability) {
      return (sqrt_s_ > mass_sum + 0.9);
    }
    /* No strings at low energy, only strings at high energy and
     * a transition region in the middle. Determine transition region: */
    double region_lower, region_upper;
    if (is_Npi_scattering) {
      region_lower = 1.9;
      region_upper = 2.2;
    } else if (is_NN_scattering) {
      region_lower = 4.0;
      region_upper = 5.0;
    } else {  // AQM - Additive Quark Model
      /* Transition region around 0.9 larger than the sum of pole masses;
       * highly arbitrary, feel free to improve */
      region_lower = mass_sum + 0.9;
      region_upper = mass_sum + 1.9;
    }

    if (sqrt_s_ > region_upper) {
      return true;
    } else if (sqrt_s_ < region_lower) {
      return false;
    } else {
      // Rescale transition region to [-1, 1]
      double x = (sqrt_s_ - (region_lower + region_upper) / 2.) /
                 (region_upper - region_lower);
      assert(x >= -0.5 && x <= 0.5);
      double prob_pythia = 0.5 * (std::sin(M_PI * x) + 1.0);
      assert(prob_pythia >= 0. && prob_pythia <= 1.);
      return prob_pythia > random::uniform(0., 1.);
    }
  }
}

}  // namespace smash
