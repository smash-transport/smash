/*
 *
 *    Copyright (c) 2013-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "include/crosssections.h"

#include "include/clebschgordan.h"
#include "include/constants.h"
#include "include/kinematics.h"
#include "include/logging.h"
#include "include/parametrizations.h"
#include "include/particletype.h"

namespace smash {

cross_sections::cross_sections(const ParticleList& scat_particles,
                               const double sqrt_s)
    : incoming_particles_(scat_particles), sqrt_s_(sqrt_s) {}

CollisionBranchList cross_sections::generate_collision_list(
    double elastic_parameter, bool two_to_one_switch, bool two_to_two_switch,
    double low_snn_cut, bool strings_switch, NNbarTreatment nnbar_treatment,
    StringProcess* string_process) {
  CollisionBranchList process_list;
  const ParticleType& t1 = incoming_particles_[0].type();
  const ParticleType& t2 = incoming_particles_[1].type();
  const bool both_are_nucleons = t1.is_nucleon() && t2.is_nucleon();

  const bool is_pythia = decide_string(strings_switch, both_are_nucleons);

  /** Elastic collisions between two nucleons with sqrt_s below
   * low_snn_cut can not happen*/
  const bool reject_by_nucleon_elastic_cutoff =
      both_are_nucleons && t1.antiparticle_sign() == t2.antiparticle_sign() &&
      sqrt_s_ < low_snn_cut;
  if (two_to_two_switch && !reject_by_nucleon_elastic_cutoff) {
    process_list.emplace_back(elastic(elastic_parameter));
  }
  if (is_pythia) {
    /* string excitation */
    append_list(process_list, string_excitation(string_process));
  } else {
    if (two_to_one_switch) {
      /* resonance formation (2->1) */
      append_list(process_list, two_to_one());
    }
    if (two_to_two_switch) {
      /* 2->2 (inelastic) */
      append_list(process_list, two_to_two());
    }
  }
  /** NNbar annihilation thru NNbar → ρh₁(1170); combined with the decays
   *  ρ → ππ and h₁(1170) → πρ, this gives a final state of 5 pions.
   *  Only use in cases when detailed balance MUST happen, i.e. in a box! */
  if (nnbar_treatment == NNbarTreatment::Resonances) {
    if (t1.is_nucleon() && t2.pdgcode() == t1.get_antiparticle()->pdgcode()) {
      /* Has to be called after the other processes are already determined,
       *  so that the sum of the cross sections includes all other processes. */
      process_list.emplace_back(NNbar_annihilation(sum_xs_of(process_list)));
    }
    if ((t1.pdgcode() == pdg::rho_z && t2.pdgcode() == pdg::h1) ||
        (t1.pdgcode() == pdg::h1 && t2.pdgcode() == pdg::rho_z)) {
      append_list(process_list, NNbar_creation());
    }
  }
  return process_list;
}

bool cross_sections::decide_string(bool strings_switch,
                                   const bool both_are_nucleons) const {
  // Determine the energy region of the mixed scattering type for two types of
  // scattering.
  const ParticleType& t1 = incoming_particles_[0].type();
  const ParticleType& t2 = incoming_particles_[1].type();
  bool include_pythia = false;
  double mix_scatter_type_energy;
  double mix_scatter_type_window_width;
  if (both_are_nucleons) {
    // The energy region of the mixed scattering type for nucleon-nucleon
    // collision is 4.0 - 5.0 GeV.
    mix_scatter_type_energy = 4.5;
    mix_scatter_type_window_width = 0.5;
    // nucleon-nucleon collisions are included in pythia.
    include_pythia = true;
  } else if ((t1.pdgcode().is_pion() && t2.is_nucleon()) ||
             (t1.is_nucleon() && t2.pdgcode().is_pion())) {
    // The energy region of the mixed scattering type for pion-nucleon collision
    // is 2.3 - 3.1 GeV.
    mix_scatter_type_energy = 2.7;
    mix_scatter_type_window_width = 0.4;
    // pion-nucleon collisions are included in pythia.
    include_pythia = true;
  }
  // string fragmentation is enabled when strings_switch is on and the process
  // is included in pythia.
  const bool enable_pythia = strings_switch && include_pythia;
  // Whether the scattering is through string fragmentaion
  bool is_pythia = false;
  if (enable_pythia) {
    if (sqrt_s_ > mix_scatter_type_energy + mix_scatter_type_window_width) {
      // scatterings at high energies are through string fragmentation
      is_pythia = true;
    } else if (sqrt_s_ >
               mix_scatter_type_energy - mix_scatter_type_window_width) {
      const double probability_pythia =
          (sqrt_s_ - mix_scatter_type_energy + mix_scatter_type_window_width) /
          mix_scatter_type_window_width / 2.0;
      if (probability_pythia > Random::uniform(0., 1.)) {
        // scatterings at the middle energies are through string
        // fragmentation by chance.
        is_pythia = true;
      }
    }
  }
  return is_pythia;
}

CollisionBranchPtr cross_sections::elastic(double elast_par) {
  double elastic_xs = 0.;
  if (elast_par >= 0.) {
    // use constant elastic cross section from config file
    elastic_xs = elast_par;
  } else {
    // use parametrization
    elastic_xs = elastic_parametrization();
  }
  return make_unique<CollisionBranch>(incoming_particles_[0].type(),
                                      incoming_particles_[1].type(), elastic_xs,
                                      ProcessType::Elastic);
}

double cross_sections::elastic_parametrization() {
  const auto& pdg_a = incoming_particles_[0].pdgcode();
  const auto& pdg_b = incoming_particles_[1].pdgcode();
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
  }
  return elastic_xs;
}

double cross_sections::nn_el() {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();

  const double s = sqrt_s_ * sqrt_s_;

  /* Use parametrized cross sections. */
  double sig_el;
  if (pdg_a == pdg_b) { /* pp */
    sig_el = pp_elastic(s);
  } else if (pdg_a.is_antiparticle_of(pdg_b)) { /* ppbar */
    sig_el = ppbar_elastic(s);
  } else { /* np */
    sig_el = np_elastic(s);
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

double cross_sections::npi_el() {
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

double cross_sections::nk_el() {
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

CollisionBranchList cross_sections::two_to_one() {
  const auto& log = logger<LogArea::CrossSections>();
  CollisionBranchList resonance_process_list;
  /* There is no resonance formation out of two baryons: Return empty list. */
  if (!(incoming_particles_[0].is_baryon() &&
        incoming_particles_[1].is_baryon())) {
    const ParticleType& type_particle_a = incoming_particles_[0].type();
    const ParticleType& type_particle_b = incoming_particles_[1].type();

    const double m1 = incoming_particles_[0].effective_mass();
    const double m2 = incoming_particles_[1].effective_mass();
    const double p_cm_sqr = pCM_sqr(sqrt_s_, m1, m2);

    /* Find all the possible resonances */
    for (const ParticleType& type_resonance : ParticleType::list_all()) {
      /* Not a resonance, go to next type of particle */
      if (type_resonance.is_stable()) {
        continue;
      }

      /* Same resonance as in the beginning, ignore */
      if ((!type_particle_a.is_stable() &&
           type_resonance.pdgcode() == type_particle_a.pdgcode()) ||
          (!type_particle_b.is_stable() &&
           type_resonance.pdgcode() == type_particle_b.pdgcode())) {
        continue;
      }

      double resonance_xsection = formation(type_resonance, p_cm_sqr);

      /* If cross section is non-negligible, add resonance to the list */
      if (resonance_xsection > really_small) {
        resonance_process_list.push_back(make_unique<CollisionBranch>(
            type_resonance, resonance_xsection, ProcessType::TwoToOne));
        log.debug("Found resonance: ", type_resonance);
        log.debug("2->1 with original particles: ", type_particle_a,
                  type_particle_b);
      }
    }
  }
  return resonance_process_list;
}

double cross_sections::formation(const ParticleType& type_resonance,
                                 double cm_momentum_sqr) {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();
  /* Check for charge conservation. */
  if (type_resonance.charge() !=
      type_particle_a.charge() + type_particle_b.charge()) {
    return 0.;
  }

  /* Check for baryon-number conservation. */
  if (type_resonance.baryon_number() !=
      type_particle_a.baryon_number() + type_particle_b.baryon_number()) {
    return 0.;
  }

  /* Calculate partial in-width. */
  const double partial_width = type_resonance.get_partial_in_width(
      sqrt_s_, incoming_particles_[0], incoming_particles_[1]);
  if (partial_width <= 0.) {
    return 0.;
  }

  /* Calculate spin factor */
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

CollisionBranchList cross_sections::two_to_two() {
  CollisionBranchList process_list;
  const ParticleData& data_a = incoming_particles_[0];
  const ParticleData& data_b = incoming_particles_[1];
  const auto& pdg_a = data_a.pdgcode();
  const auto& pdg_b = data_b.pdgcode();
  if (data_a.is_baryon() && data_b.is_baryon()) {
    if (pdg_a.is_nucleon() && pdg_b.is_nucleon() &&
        pdg_a.antiparticle_sign() == pdg_b.antiparticle_sign()) {
      // Nucleon Nucleon Scattering
      process_list = nn_xx();
    } else {
      // Baryon Baryon Scattering
      process_list = bb_xx_except_nn();
    }
  } else if (data_a.is_baryon() || data_b.is_baryon()) {
    // Baryon Meson Scattering
    if ((pdg_a.is_nucleon() && pdg_b.is_kaon()) ||
        (pdg_b.is_nucleon() && pdg_a.is_kaon())) {
      // Nucleon Kaon Scattering
      process_list = nk_xx();
    } else if ((pdg_a.is_hyperon() && pdg_b.is_pion()) ||
               (pdg_b.is_hyperon() && pdg_a.is_pion())) {
      // Hyperon Pion Scattering
      process_list = ypi_xx();
    } else if ((pdg_a.is_Delta() && pdg_b.is_kaon()) ||
               (pdg_b.is_Delta() && pdg_a.is_kaon())) {
      // Delta Kaon Scattering
      process_list = deltak_xx();
    }
  }
  return process_list;
}

CollisionBranchList cross_sections::bb_xx_except_nn() {
  CollisionBranchList process_list;
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();

  if (type_a.is_nucleon() || type_a.is_Delta() || type_b.is_nucleon() ||
      type_b.is_Delta()) {
    if (type_a.antiparticle_sign() == 1 && type_b.antiparticle_sign() == 1) {
      /* N R → N N, Δ R → N N */
      process_list = bar_bar_to_nuc_nuc(false);
    } else if (type_a.antiparticle_sign() == -1 &&
               type_b.antiparticle_sign() == -1) {
      /* N̅ R → N̅ N̅, Δ̅ R → N̅ N̅ */
      process_list = bar_bar_to_nuc_nuc(true);
    }
  }

  return process_list;
}

CollisionBranchList cross_sections::nn_xx() {
  CollisionBranchList process_list, channel_list;

  const double sqrts = sqrt_s_;

  /* Find whether colliding particles are nucleons or anti-nucleons;
   * adjust lists of produced particles. */
  const ParticleTypePtrList& nuc_or_anti_nuc =
      incoming_particles_[0].type().antiparticle_sign() == -1 &&
              incoming_particles_[1].type().antiparticle_sign() == -1
          ? ParticleType::list_anti_nucleons()
          : ParticleType::list_nucleons();
  const ParticleTypePtrList& delta_or_anti_delta =
      incoming_particles_[0].type().antiparticle_sign() == -1 &&
              incoming_particles_[1].type().antiparticle_sign() == -1
          ? ParticleType::list_anti_Deltas()
          : ParticleType::list_Deltas();
  /* First: Find N N → N R channels. */
  channel_list = find_nn_xsection_from_type(
      ParticleType::list_baryon_resonances(), nuc_or_anti_nuc,
      [&sqrts](const ParticleType& type_res_1, const ParticleType&) {
        return type_res_1.iso_multiplet()->get_integral_NR(sqrts);
      });
  process_list.reserve(process_list.size() + channel_list.size());
  std::move(channel_list.begin(), channel_list.end(),
            std::inserter(process_list, process_list.end()));
  channel_list.clear();

  /* Second: Find N N → Δ R channels. */
  channel_list = find_nn_xsection_from_type(
      ParticleType::list_baryon_resonances(), delta_or_anti_delta,
      [&sqrts](const ParticleType& type_res_1, const ParticleType& type_res_2) {
        return type_res_1.iso_multiplet()->get_integral_RR(type_res_2, sqrts);
      });
  process_list.reserve(process_list.size() + channel_list.size());
  std::move(channel_list.begin(), channel_list.end(),
            std::inserter(process_list, process_list.end()));
  channel_list.clear();

  return process_list;
}

CollisionBranchList cross_sections::nk_xx() {
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

  CollisionBranchList process_list;
  switch (pdg_kaon) {
    case pdg::K_m: {
      // All inelastic K- N channels here are strangeness exchange, plus one
      // charge exchange.
      switch (pdg_nucleon) {
        case pdg::p: {
          const auto& type_n = ParticleType::find(pdg::n);
          const auto& type_pi_z = ParticleType::find(pdg::pi_z);
          const auto& type_pi_m = ParticleType::find(pdg::pi_m);
          const auto& type_pi_p = ParticleType::find(pdg::pi_p);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
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
          add_channel(process_list, [&] { return kminusp_pi0sigma0(sqrt_s_); },
                      sqrt_s_, type_pi_z, type_Sigma_z);
          add_channel(process_list, [&] { return kminusp_pi0lambda(sqrt_s_); },
                      sqrt_s_, type_pi_z, type_Lambda);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                      type_Kbar_z, type_n);
          break;
        }
        case pdg::n: {
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
          break;
        }
        case -pdg::p: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_pp_bar = ParticleType::find(-pdg::Delta_pp);
          const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_pp_bar);
                      },
                      sqrt_s_, type_Kbar_z, type_Delta_pp_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp * kplusn_ratios.get_ratio(
                                                  type_nucleon, type_kaon,
                                                  type_K_m, type_Delta_p_bar);
                      },
                      sqrt_s_, type_K_m, type_Delta_p_bar);
          break;
        }
        case -pdg::n: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
          const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
          const auto& type_p_bar = ParticleType::find(-pdg::p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_p_bar);
                      },
                      sqrt_s_, type_Kbar_z, type_Delta_p_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn * kplusn_ratios.get_ratio(
                                                  type_nucleon, type_kaon,
                                                  type_K_m, type_Delta_z_bar);
                      },
                      sqrt_s_, type_K_m, type_Delta_z_bar);
          add_channel(process_list, [&] { return kplusn_k0p(s); }, sqrt_s_,
                      type_Kbar_z, type_p_bar);
          break;
        }
      }
      break;
    }
    case pdg::K_p: {
      // All inelastic channels are K+ N -> K Delta -> K pi N, with identical
      // cross section, weighted by the isospin factor.
      switch (pdg_nucleon) {
        case pdg::p: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_Delta_pp = ParticleType::find(pdg::Delta_pp);
          const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_pp);
                      },
                      sqrt_s_, type_K_z, type_Delta_pp);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_p);
                      },
                      sqrt_s_, type_K_p, type_Delta_p);
          break;
        }
        case pdg::n: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_p = ParticleType::find(pdg::p);
          const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
          const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_p);
                      },
                      sqrt_s_, type_K_z, type_Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_z);
                      },
                      sqrt_s_, type_K_p, type_Delta_z);
          add_channel(process_list, [&] { return kplusn_k0p(s); }, sqrt_s_,
                      type_K_z, type_p);
          break;
        }
        case -pdg::p: {
          const auto& type_n_bar = ParticleType::find(-pdg::n);
          const auto& type_pi_z = ParticleType::find(pdg::pi_z);
          const auto& type_pi_m = ParticleType::find(pdg::pi_m);
          const auto& type_pi_p = ParticleType::find(pdg::pi_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
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
          add_channel(process_list, [&] { return kminusp_pi0sigma0(sqrt_s_); },
                      sqrt_s_, type_pi_z, type_Sigma_z_bar);
          add_channel(process_list, [&] { return kminusp_pi0lambda(sqrt_s_); },
                      sqrt_s_, type_pi_z, type_Lambda_bar);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                      type_K_z, type_n_bar);
          break;
        }
        case -pdg::n: {
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
          break;
        }
      }
      break;
    }
    case pdg::K_z: {
      // K+ and K0 have the same isospin projection, they are assumed to have
      // the same cross section here.

      switch (pdg_nucleon) {
        case pdg::p: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_n = ParticleType::find(pdg::n);
          const auto& type_Delta_p = ParticleType::find(pdg::Delta_p);
          const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_p);
                      },
                      sqrt_s_, type_K_z, type_Delta_p);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_z);
                      },
                      sqrt_s_, type_K_p, type_Delta_z);
          add_channel(process_list,
                      [&] {
                        return kplusn_k0p(s) *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_n);
                      },
                      sqrt_s_, type_K_p, type_n);
          break;
        }
        case pdg::n: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_K_z = ParticleType::find(pdg::K_z);
          const auto& type_Delta_z = ParticleType::find(pdg::Delta_z);
          const auto& type_Delta_m = ParticleType::find(pdg::Delta_m);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_z, type_Delta_z);
                      },
                      sqrt_s_, type_K_z, type_Delta_z);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_p, type_Delta_m);
                      },
                      sqrt_s_, type_K_p, type_Delta_m);
          break;
        }
        case -pdg::n: {
          const auto& type_K_p = ParticleType::find(pdg::K_p);
          const auto& type_p_bar = ParticleType::find(-pdg::p);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                      type_K_p, type_p_bar);
          break;
        }
      }
      break;
    }
    case pdg::Kbar_z:
      switch (pdg_nucleon) {
        case pdg::n: {
          const auto& type_p = ParticleType::find(pdg::p);
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          add_channel(process_list, [&] { return kminusp_kbar0n(s); }, sqrt_s_,
                      type_K_m, type_p);
          break;
        }
        case -pdg::p: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_p_bar = ParticleType::find(-pdg::Delta_p);
          const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
          const auto& type_n_bar = ParticleType::find(-pdg::n);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_p_bar);
                      },
                      sqrt_s_, type_Kbar_z, type_Delta_p_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusp * kplusn_ratios.get_ratio(
                                                  type_nucleon, type_kaon,
                                                  type_K_m, type_Delta_z_bar);
                      },
                      sqrt_s_, type_K_m, type_Delta_z_bar);
          add_channel(process_list,
                      [&] {
                        return kplusn_k0p(s) *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_K_m, type_n_bar);
                      },
                      sqrt_s_, type_K_m, type_n_bar);
          break;
        }
        case -pdg::n: {
          const auto& type_K_m = ParticleType::find(pdg::K_m);
          const auto& type_Kbar_z = ParticleType::find(pdg::Kbar_z);
          const auto& type_Delta_z_bar = ParticleType::find(-pdg::Delta_z);
          const auto& type_Delta_m_bar = ParticleType::find(-pdg::Delta_m);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn *
                               kplusn_ratios.get_ratio(type_nucleon, type_kaon,
                                                       type_Kbar_z,
                                                       type_Delta_z_bar);
                      },
                      sqrt_s_, type_Kbar_z, type_Delta_z_bar);
          add_channel(process_list,
                      [&] {
                        return sigma_kplusn * kplusn_ratios.get_ratio(
                                                  type_nucleon, type_kaon,
                                                  type_K_m, type_Delta_m_bar);
                      },
                      sqrt_s_, type_K_m, type_Delta_m_bar);
          break;
        }
      }
      break;
  }

  return process_list;
}

CollisionBranchList cross_sections::deltak_xx() {
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_delta = a.pdgcode().is_Delta() ? a : b;
  const ParticleType& type_kaon = a.pdgcode().is_Delta() ? b : a;

  const auto pdg_delta = type_delta.pdgcode().code();
  const auto pdg_kaon = type_kaon.pdgcode().code();

  const double s = sqrt_s_ * sqrt_s_;
  const double pcm = cm_momentum();

  CollisionBranchList process_list;
  // The cross sections are determined from the backward reactions via
  // detailed balance. The same isospin factors as for the backward reaction
  // are used.
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
                           kplusn_ratios.get_ratio(type_p, type_K_p, type_kaon,
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
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p_bar,
                                                      type_K_m) *
                           kplusn_ratios.get_ratio(type_p_bar, type_K_m,
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
      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_n,
                                                      type_K_p) *
                           kplusn_ratios.get_ratio(type_n, type_K_p, type_kaon,
                                                   type_delta) *
                           kplusn_inelastic_background(s);
                  },
                  sqrt_s_, type_n, type_K_p);

      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p,
                                                      type_K_z) *
                           kplusn_ratios.get_ratio(type_p, type_K_z, type_kaon,
                                                   type_delta) *
                           kplusp_inelastic_background(s);
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
                           kplusn_ratios.get_ratio(type_n_bar, type_K_m,
                                                   type_kaon, type_delta) *
                           kplusn_inelastic_background(s);
                  },
                  sqrt_s_, type_n_bar, type_K_m);

      add_channel(process_list,
                  [&] {
                    return detailed_balance_factor_RK(sqrt_s_, pcm, type_delta,
                                                      type_kaon, type_p_bar,
                                                      type_Kbar_z) *
                           kplusn_ratios.get_ratio(type_p_bar, type_Kbar_z,
                                                   type_kaon, type_delta) *
                           kplusp_inelastic_background(s);
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
                           kplusn_ratios.get_ratio(type_n, type_K_z, type_kaon,
                                                   type_delta) *
                           kplusn_inelastic_background(s);
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
                           kplusn_ratios.get_ratio(type_n_bar, type_Kbar_z,
                                                   type_kaon, type_delta) *
                           kplusn_inelastic_background(s);
                  },
                  sqrt_s_, type_n_bar, type_Kbar_z);
      break;
    }
    default:
      break;
  }

  return process_list;
}

CollisionBranchList cross_sections::ypi_xx() {
  const ParticleType& a = incoming_particles_[0].type();
  const ParticleType& b = incoming_particles_[1].type();
  const ParticleType& type_hyperon = a.pdgcode().is_hyperon() ? a : b;
  const ParticleType& type_pion = a.pdgcode().is_hyperon() ? b : a;

  const auto pdg_hyperon = type_hyperon.pdgcode().code();
  const auto pdg_pion = type_pion.pdgcode().code();

  const double s = sqrt_s_ * sqrt_s_;

  CollisionBranchList process_list;

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

double cross_sections::high_energy() const {
  const PdgCode& pdg_a = incoming_particles_[0].type().pdgcode();
  const PdgCode& pdg_b = incoming_particles_[1].type().pdgcode();
  const double s = sqrt_s_ * sqrt_s_;

  /* Currently all BB collisions use the nucleon-nucleon parametrizations. */
  if (pdg_a.is_baryon() && pdg_b.is_baryon()) {
    if (pdg_a == pdg_b) {
      return pp_high_energy(s);  // pp, nn
    } else if (pdg_a.is_antiparticle_of(pdg_b)) {
      return ppbar_high_energy(s);  // ppbar, nnbar
    } else if (pdg_a.antiparticle_sign() * pdg_b.antiparticle_sign() == 1) {
      return np_high_energy(s);  // np, nbarpbar
    } else {
      return npbar_high_energy(s);  // npbar, nbarp
    }
  }

  /* Pion nucleon interaction. */
  if ((pdg_a == pdg::pi_p && pdg_b == pdg::p) ||
      (pdg_b == pdg::pi_p && pdg_a == pdg::p) ||
      (pdg_a == pdg::pi_m && pdg_b == pdg::n) ||
      (pdg_b == pdg::pi_m && pdg_a == pdg::n)) {
    return piplusp_high_energy(s);  // pi+ p, pi- n
  } else if ((pdg_a == pdg::pi_m && pdg_b == pdg::p) ||
             (pdg_b == pdg::pi_m && pdg_a == pdg::p) ||
             (pdg_a == pdg::pi_p && pdg_b == pdg::n) ||
             (pdg_b == pdg::pi_p && pdg_a == pdg::n)) {
    return piminusp_high_energy(s);  // pi- p, pi+ n
  } else {
    return 0;
  }
}

CollisionBranchList cross_sections::string_excitation(
    StringProcess* string_process) {
  const auto& log = logger<LogArea::CrossSections>();
  /* Calculate string-excitation cross section:
   * Parametrized total minus all other present channels. */
  double sig_string_all =
      std::max(0., high_energy() - elastic_parametrization());

  /* get PDG id for evaluation of the parametrized cross sections
   * for diffractive processes.
   * (anti-)proton is used for (anti-)baryons and
   * pion is used for mesons.
   * This must be rescaled according to additive quark model
   * in the case of exotic hadrons. */
  std::array<int, 2> pdgid;
  for (int i = 0; i < 2; i++) {
    PdgCode pdg = incoming_particles_[i].type().pdgcode();
    pdg.deexcite();
    if (pdg.baryon_number() == 1) {
      pdgid[i] = 2212;
    } else if (pdg.baryon_number() == -1) {
      pdgid[i] = -2212;
    } else {
      pdgid[i] = 211;
    }
  }

  CollisionBranchList channel_list;
  if (sig_string_all > 0.) {
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
    if (!string_process) {
      throw std::runtime_error("string_process should be initialized.");
    }
    std::array<double, 3> xs =
        string_process->cross_sections_diffractive(pdgid[0], pdgid[1], sqrt_s_);
    double single_diffr_AX = xs[0], single_diffr_XB = xs[1],
           double_diffr = xs[2];
    double single_diffr = single_diffr_AX + single_diffr_XB;
    double diffractive = single_diffr + double_diffr;
    const double nondiffractive_all =
        std::max(0., sig_string_all - diffractive);
    diffractive = sig_string_all - nondiffractive_all;
    double_diffr = std::max(0., diffractive - single_diffr);
    const double a = (diffractive - double_diffr) / single_diffr;
    single_diffr_AX *= a;
    single_diffr_XB *= a;
    assert(std::abs(single_diffr_AX + single_diffr_XB + double_diffr +
                    nondiffractive_all - sig_string_all) < 1.e-6);

    /* Hard string process is added by hard cross section
     * in conjunction with multipartion interaction picture
     * \iref{Sjostrand:1987su}. */
    const double hard_xsec = string_hard_cross_section();
    const double nondiffractive_soft =
        nondiffractive_all * std::exp(-hard_xsec / nondiffractive_all);
    const double nondiffractive_hard = nondiffractive_all - nondiffractive_soft;
    log.debug("String cross sections [mb] are");
    log.debug("Single-diffractive AB->AX: ", single_diffr_AX);
    log.debug("Single-diffractive AB->XB: ", single_diffr_XB);
    log.debug("Double-diffractive AB->XX: ", double_diffr);
    log.debug("Soft non-diffractive: ", nondiffractive_soft);
    log.debug("Hard non-diffractive: ", nondiffractive_hard);

    /* cross section of soft string excitation */
    const double sig_string_soft = sig_string_all - nondiffractive_hard;

    /* fill cross section arrays */
    std::array<double, 5> string_sub_cross_sections;
    std::array<double, 6> string_sub_cross_sections_sum;
    string_sub_cross_sections[0] = single_diffr_AX;
    string_sub_cross_sections[1] = single_diffr_XB;
    string_sub_cross_sections[2] = double_diffr;
    string_sub_cross_sections[3] = nondiffractive_soft;
    string_sub_cross_sections[4] = nondiffractive_hard;
    string_sub_cross_sections_sum[0] = 0.;
    for (int i = 0; i < 5; i++) {
      string_sub_cross_sections_sum[i + 1] =
          string_sub_cross_sections_sum[i] + string_sub_cross_sections[i];
    }

    /* soft subprocess selection */
    int iproc = -1;
    double r_xsec = string_sub_cross_sections_sum[4] * Random::uniform(0., 1.);
    for (int i = 0; i < 4; i++) {
      if ((r_xsec >= string_sub_cross_sections_sum[i]) &&
          (r_xsec < string_sub_cross_sections_sum[i + 1])) {
        iproc = i;
        break;
      }
    }
    if (iproc == -1) {
      throw std::runtime_error("soft string subprocess is not specified.");
    }

    string_process->set_iproc(iproc);

    /* fill the list of process channels */
    if (sig_string_soft > 0.) {
      channel_list.push_back(make_unique<CollisionBranch>(
          sig_string_soft, ProcessType::StringSoft));
    }
    if (nondiffractive_hard > 0.) {
      channel_list.push_back(make_unique<CollisionBranch>(
          nondiffractive_hard, ProcessType::StringHard));
    }
  }
  return channel_list;
}

double cross_sections::string_hard_cross_section() const {
  double cross_sec = 0.;
  const ParticleData& data_a = incoming_particles_[0];
  const ParticleData& data_b = incoming_particles_[1];
  if (data_a.is_baryon() && data_b.is_baryon()) {
    /**
     * Currently nucleon-nucleon cross section is used for all baryon-baryon
     * casees. This will be changed later by applying additive quark model.
     */
    cross_sec = NN_string_hard(sqrt_s_ * sqrt_s_);
  } else if (data_a.is_baryon() || data_b.is_baryon()) {
    /**
     * Currently nucleon-pion cross section is used for all baryon-meson cases.
     * This will be changed later by applying additive quark model.
     */
    cross_sec = Npi_string_hard(sqrt_s_ * sqrt_s_);
  } else {
    /**
     * Currently pion-pion cross section is used for all meson-meson cases.
     * This will be changed later by applying additive quark model.
     */
    cross_sec = pipi_string_hard(sqrt_s_ * sqrt_s_);
  }
  return cross_sec;
}

CollisionBranchPtr cross_sections::NNbar_annihilation(const double current_xs) {
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

CollisionBranchList cross_sections::NNbar_creation() {
  const auto& log = logger<LogArea::CrossSections>();
  CollisionBranchList channel_list;
  /* Calculate NNbar reverse cross section:
   * from reverse reaction (see NNbar_annihilation_cross_section).*/
  const double s = sqrt_s_ * sqrt_s_;
  const double sqrts = sqrt_s_;
  const double pcm = cm_momentum();

  const auto& type_N = ParticleType::find(pdg::p);
  const auto& type_Nbar = ParticleType::find(-pdg::p);

  // Check available energy
  if (sqrts - 2 * type_N.mass() < 0) {
    return channel_list;
  }

  double xsection = detailed_balance_factor_RR(
                        sqrts, pcm, incoming_particles_[0].type(),
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

CollisionBranchList cross_sections::bar_bar_to_nuc_nuc(
    const bool is_anti_particles) {
  const ParticleType& type_a = incoming_particles_[0].type();
  const ParticleType& type_b = incoming_particles_[1].type();
  CollisionBranchList process_list;

  const double s = sqrt_s_ * sqrt_s_;
  /* CM momentum in final state */
  double p_cm_final = std::sqrt(s - 4. * nucleon_mass * nucleon_mass) / 2.;

  ParticleTypePtrList nuc_or_anti_nuc;
  if (is_anti_particles) {
    nuc_or_anti_nuc = ParticleType::list_anti_nucleons();
  } else {
    nuc_or_anti_nuc = ParticleType::list_nucleons();
  }

  /* Loop over all nucleon or anti-nucleon charge states. */
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
        /* If Clebsch-Gordan coefficient is zero, don't bother with the rest */
        if (std::abs(isospin_factor) < really_small) {
          continue;
        }

        /* Calculate matrix element for inverse process. */
        const double matrix_element =
            nn_to_resonance_matrix_element(type_a, type_b, twoI);
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

double cross_sections::nn_to_resonance_matrix_element(
    const ParticleType& type_a, const ParticleType& type_b,
    const int twoI) const {
  const double m_a = type_a.mass();
  const double m_b = type_b.mass();
  const double msqr = 2. * (m_a * m_a + m_b * m_b);
  /* If the c.m. energy is larger than the sum of the pole masses of the
   * outgoing particles plus three times of the sum of the widths plus 3 GeV,
   * the collision will be neglected.*/
  const double w_a = type_a.width_at_pole();
  const double w_b = type_b.width_at_pole();
  const double uplmt = m_a + m_b + 3.0 * (w_a + w_b) + 3.0;
  if (sqrt_s_ > uplmt) {
    return 0.;
  }
  /** NN → NΔ: fit sqrt(s)-dependence to OBE model [\iref{Dmitriev:1986st}] */
  if (((type_a.is_Delta() && type_b.is_nucleon()) ||
       (type_b.is_Delta() && type_a.is_nucleon())) &&
      (type_a.antiparticle_sign() == type_b.antiparticle_sign())) {
    return 68. / std::pow(sqrt_s_ - 1.104, 1.951);
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
      /* pn → pnη cross section is known to be larger than the corresponding
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
  }
  // all cases not listed: zero!
  return 0.;
}

template <class IntegrationMethod>
CollisionBranchList cross_sections::find_nn_xsection_from_type(
    const ParticleTypePtrList& list_res_1,
    const ParticleTypePtrList& list_res_2, const IntegrationMethod integrator) {
  const ParticleType& type_particle_a = incoming_particles_[0].type();
  const ParticleType& type_particle_b = incoming_particles_[1].type();

  const auto& log = logger<LogArea::CrossSections>();
  CollisionBranchList channel_list;
  const double s = sqrt_s_ * sqrt_s_;

  /* Loop over specified first resonance list */
  for (ParticleTypePtr type_res_1 : list_res_1) {
    /* Loop over specified second resonance list */
    for (ParticleTypePtr type_res_2 : list_res_2) {
      /* Check for charge conservation. */
      if (type_res_1->charge() + type_res_2->charge() !=
          type_particle_a.charge() + type_particle_b.charge()) {
        continue;
      }

      // loop over total isospin
      for (const int twoI : I_tot_range(type_particle_a, type_particle_b)) {
        const double isospin_factor = isospin_clebsch_gordan_sqr_2to2(
            type_particle_a, type_particle_b, *type_res_1, *type_res_2, twoI);
        /* If Clebsch-Gordan coefficient is zero, don't bother with the rest. */
        if (std::abs(isospin_factor) < really_small) {
          continue;
        }

        /* Integration limits. */
        const double lower_limit = type_res_1->min_mass_kinematic();
        const double upper_limit = sqrt_s_ - type_res_2->mass();
        /* Check the available energy (requiring it to be a little above the
         * threshold, because the integration will not work if it's too close).
         */
        if (upper_limit - lower_limit < 1E-3) {
          continue;
        }

        /* Calculate matrix element. */
        const double matrix_element =
            nn_to_resonance_matrix_element(*type_res_1, *type_res_2, twoI);
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

double cross_sections::cm_momentum() const {
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  return pCM(sqrt_s_, m1, m2);
}

}  // namespace smash
