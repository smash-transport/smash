/*
 *
 *    Copyright (c) 2014-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteraction.h"

#include <cmath>

#include "Pythia8/Pythia.h"

#include "smash/angles.h"
#include "smash/constants.h"
#include "smash/crosssections.h"
#include "smash/fpenvironment.h"
#include "smash/logging.h"
#include "smash/pdgcode.h"
#include "smash/pow.h"
#include "smash/random.h"

namespace smash {
static constexpr int LScatterAction = LogArea::ScatterAction::id;

ScatterAction::ScatterAction(const ParticleData &in_part_a,
                             const ParticleData &in_part_b, const double time,
                             const bool isotropic,
                             const double string_formation_time,
                             const double box_length,
                             const bool is_total_parametrized,
                             const SpinInteractionType spin_interaction_type)
    : Action({in_part_a, in_part_b}, time),
      sum_of_partial_cross_sections_(0.),
      isotropic_(isotropic),
      string_formation_time_(string_formation_time),
      is_total_parametrized_(is_total_parametrized),
      spin_interaction_type_(spin_interaction_type) {
  box_length_ = box_length;
  if (is_total_parametrized_) {
    parametrized_total_cross_section_ = smash_NaN<double>;
  }
}

void ScatterAction::add_collision(CollisionBranchPtr p) {
  add_process<CollisionBranch>(p, collision_channels_,
                               sum_of_partial_cross_sections_);
}

void ScatterAction::add_collisions(CollisionBranchList pv) {
  add_processes<CollisionBranch>(std::move(pv), collision_channels_,
                                 sum_of_partial_cross_sections_);
}

void ScatterAction::generate_final_state() {
  logg[LScatterAction].debug("Incoming particles: ", incoming_particles_);

  const CollisionBranch *proc = choose_channel<CollisionBranch>(
      collision_channels_, is_total_parametrized_
                               ? *parametrized_total_cross_section_
                               : sum_of_partial_cross_sections_);
  process_type_ = proc->get_type();
  outgoing_particles_ = proc->particle_list();
  partial_cross_section_ = proc->weight();

  logg[LScatterAction].debug("Chosen channel: ", process_type_,
                             outgoing_particles_);

  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  switch (process_type_) {
    case ProcessType::Elastic:
      /* 2->2 elastic scattering */
      elastic_scattering();
      spin_interaction();
      break;
    case ProcessType::TwoToOne:
      /* resonance formation */
      resonance_formation();
      break;
    case ProcessType::TwoToTwo:
      /* 2->2 inelastic scattering */
      /* Sample the particle momenta in CM system. */
      inelastic_scattering();
      spin_interaction();
      break;
    case ProcessType::TwoToThree:
    case ProcessType::TwoToFour:
    case ProcessType::TwoToFive:
      /* 2->m scattering */
      two_to_many_scattering();
      break;
    case ProcessType::StringSoftSingleDiffractiveAX:
    case ProcessType::StringSoftSingleDiffractiveXB:
    case ProcessType::StringSoftDoubleDiffractive:
    case ProcessType::StringSoftAnnihilation:
    case ProcessType::StringSoftNonDiffractive:
    case ProcessType::StringHard:
      string_excitation();
      break;
    default:
      throw InvalidScatterAction(
          "ScatterAction::generate_final_state: Invalid process type " +
          std::to_string(static_cast<int>(process_type_)) + " was requested. " +
          "(PDGcode1=" + incoming_particles_[0].pdgcode().string() +
          ", PDGcode2=" + incoming_particles_[1].pdgcode().string() + ")");
  }

  const bool core_in_incoming =
      std::any_of(incoming_particles_.begin(), incoming_particles_.end(),
                  [](const ParticleData &p) { return p.is_core(); });
  for (ParticleData &new_particle : outgoing_particles_) {
    // Boost to the computational frame
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
    /* Set positions of the outgoing particles */
    if (proc->get_type() != ProcessType::Elastic) {
      new_particle.set_4position(middle_point);
      if (core_in_incoming) {
        new_particle.fluidize();
      }
      if (new_particle.type().pdgcode().is_heavy_flavor()) {
        // New particle HF weight is the product of incoming weights
        const double HF_weight = std::accumulate(
            incoming_particles_.begin(), incoming_particles_.end(), 1.0,
            [](double w, const ParticleData &p) {
              return w * p.perturbative_HF_weight();
            });
        new_particle.set_perturbative_HF_weight(HF_weight);
      }
    }
  }
}

void ScatterAction::add_all_scatterings(
    const ScatterActionsFinderParameters &finder_parameters) {
  if (were_processes_added_) {
    logg[LScatterAction].fatal() << "Trying to add processes again.";
    throw std::logic_error(
        "add_all_scatterings should be called only once per ScatterAction "
        "instance");
  } else {
    were_processes_added_ = true;
  }
  CrossSections xs(incoming_particles_, sqrt_s(),
                   get_potential_at_interaction_point());
  CollisionBranchList processes =
      xs.generate_collision_list(finder_parameters, string_process_);

  /* Add various subprocesses.*/
  add_collisions(std::move(processes));

  /* If the string processes are not triggered by a probability, then they
   * always happen as long as the parametrized total cross section is larger
   * than the sum of the cross sections of the non-string processes, and the
   * square root s exceeds the threshold by at least 0.9 GeV. The cross section
   * of the string processes are counted by taking the difference between the
   * parametrized total and the sum of the non-strings. */
  if (!finder_parameters.strings_with_probability &&
      xs.string_probability(finder_parameters) > 0) {
    const double xs_diff =
        xs.high_energy(finder_parameters) - sum_of_partial_cross_sections_;
    if (xs_diff > 0.) {
      add_collisions(
          xs.string_excitation(xs_diff, string_process_, finder_parameters));
    }
  }

  ParticleTypePtr pseudoresonance =
      try_find_pseudoresonance(finder_parameters.pseudoresonance_method,
                               finder_parameters.transition_high_energy);
  if (pseudoresonance && finder_parameters.two_to_one) {
    const double xs_total = is_total_parametrized_
                                ? *parametrized_total_cross_section_
                                : xs.high_energy(finder_parameters);
    const double xs_gap = xs_total - sum_of_partial_cross_sections_;
    // The pseudo-resonance is only created if there is a (positive) cross
    // section gap
    if (xs_gap > really_small) {
      auto pseudoresonance_branch = std::make_unique<CollisionBranch>(
          *pseudoresonance, xs_gap, ProcessType::TwoToOne);
      add_collision(std::move(pseudoresonance_branch));
      logg[LScatterAction].debug()
          << "Pseudoresonance between " << incoming_particles_[0].type().name()
          << " and " << incoming_particles_[1].type().name() << "is"
          << pseudoresonance->name() << " with cross section " << xs_gap
          << "mb.";
    }
  }
  were_processes_added_ = true;
  // Rescale the branches so that their sum matches the parametrization
  if (is_total_parametrized_) {
    rescale_outgoing_branches();
  }
}

void ScatterAction::rescale_outgoing_branches() {
  if (!were_processes_added_) {
    logg[LScatterAction].fatal()
        << "Trying to rescale branches before adding processes.";
    throw std::logic_error(
        "This function can only be called after having added processes.");
  }
  if (sum_of_partial_cross_sections_ < really_small) {
    const ParticleTypePtr type_a = &incoming_particles_[0].type();
    const ParticleTypePtr type_b = &incoming_particles_[1].type();
    // This is a std::set instead of std::pair because the order of particles
    // does not matter here
    const std::set<ParticleTypePtr> pair{type_a, type_b};
    if (!warned_no_rescaling_available.count(pair)) {
      logg[LScatterAction].warn()
          << "Total cross section between " << type_a->name() << " and "
          << type_b->name() << "is roughly zero at sqrt(s) = " << sqrt_s()
          << " GeV, and no rescaling to match the parametrized value will be "
             "done.\nAn elastic process will be added, instead, to match the "
             "total cross section.\nFor this pair of particles, this warning "
             "will be subsequently suppressed.";
      warned_no_rescaling_available.insert(pair);
    }
    auto elastic_branch = std::make_unique<CollisionBranch>(
        *type_a, *type_b, *parametrized_total_cross_section_,
        ProcessType::Elastic);
    add_collision(std::move(elastic_branch));
  } else {
    const double reweight =
        *parametrized_total_cross_section_ / sum_of_partial_cross_sections_;
    logg[LScatterAction].debug("Reweighting ", sum_of_partial_cross_sections_,
                               " to ", *parametrized_total_cross_section_);
    for (auto &proc : collision_channels_) {
      proc->set_weight(proc->weight() * reweight);
    }
  }
}

ParticleTypePtr ScatterAction::try_find_pseudoresonance(
    const PseudoResonance method,
    const StringTransitionParameters &transition) const {
  const ParticleTypePtr type_a = &incoming_particles_[0].type();
  const ParticleTypePtr type_b = &incoming_particles_[1].type();
  const double desired_mass = sqrt_s();

  double string_offset = 0.5 * transition.sqrts_add_lower;
  const bool nucleon_and_pion = (type_a->is_nucleon() && type_b->is_pion()) ||
                                (type_a->is_pion() && type_b->is_nucleon());
  const bool two_nucleons = type_a->is_nucleon() && type_b->is_nucleon();
  const bool two_pions = type_a->is_pion() && type_b->is_pion();
  const bool nucleon_and_kaon = (type_a->is_nucleon() && type_b->is_kaon()) ||
                                (type_a->is_kaon() && type_b->is_nucleon());
  if (nucleon_and_pion) {
    string_offset = transition.sqrts_range_Npi.first - pion_mass - nucleon_mass;
  } else if (two_nucleons) {
    string_offset = transition.sqrts_range_NN.first - 2 * nucleon_mass;
  } else if (two_pions) {
    string_offset = transition.pipi_offset;
  } else if (nucleon_and_kaon) {
    string_offset = transition.KN_offset;
  }
  /*
   * Artificial cutoff to create a pseudo-resonance only close to the string
   * transition, where data or first-principle models are unhelpful
   */
  if (desired_mass < incoming_particles_[0].effective_mass() +
                         incoming_particles_[1].effective_mass() +
                         string_offset) {
    return {};
  }

  if (method == PseudoResonance::None) {
    return {};
  } else if (method == PseudoResonance::LargestFromUnstable ||
             method == PseudoResonance::ClosestFromUnstable) {
    if (type_a->is_stable() && type_b->is_stable()) {
      return {};
    }
  }

  // If this list is empty, there are no possible pseudo-resonances.
  ParticleTypePtrList list = list_possible_resonances(type_a, type_b);
  if (std::empty(list)) {
    return {};
  }

  if (method == PseudoResonance::Largest ||
      method == PseudoResonance::LargestFromUnstable) {
    auto largest = *std::max_element(list.begin(), list.end(),
                                     [](ParticleTypePtr a, ParticleTypePtr b) {
                                       return a->mass() < b->mass();
                                     });
    return largest;
  } else if (method == PseudoResonance::Closest ||
             method == PseudoResonance::ClosestFromUnstable) {
    auto comparison = [&desired_mass](ParticleTypePtr a, ParticleTypePtr b) {
      return std::abs(a->mass() - desired_mass) <
             std::abs(b->mass() - desired_mass);
    };
    auto closest = *std::min_element(list.begin(), list.end(), comparison);
    return closest;
  } else {
    throw std::logic_error("Unknown method for selecting pseudoresonance.");
  }
}

void ScatterAction::set_parametrized_total_cross_section(
    const ScatterActionsFinderParameters &finder_parameters) {
  CrossSections xs(incoming_particles_, sqrt_s(),
                   get_potential_at_interaction_point());

  if (is_total_parametrized_) {
    parametrized_total_cross_section_ =
        xs.parametrized_total(finder_parameters);
  } else {
    logg[LScatterAction].fatal()
        << "Trying to parametrize total cross section when it shouldn't be.";
    throw std::logic_error(
        "This function can only be called on ScatterAction objects with "
        "parametrized cross section.");
  }
}

double ScatterAction::get_total_weight() const {
  return sum_of_partial_cross_sections_ *
         incoming_particles_[0].xsec_scaling_factor() *
         incoming_particles_[1].xsec_scaling_factor();
}

double ScatterAction::get_partial_weight() const {
  return partial_cross_section_ * incoming_particles_[0].xsec_scaling_factor() *
         incoming_particles_[1].xsec_scaling_factor();
}

ThreeVector ScatterAction::beta_cm() const {
  return total_momentum().velocity();
}

double ScatterAction::gamma_cm() const {
  return (1. / std::sqrt(1.0 - beta_cm().sqr()));
}

double ScatterAction::mandelstam_s() const { return total_momentum().sqr(); }

double ScatterAction::cm_momentum() const {
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  return pCM(sqrt_s(), m1, m2);
}

double ScatterAction::cm_momentum_squared() const {
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  return pCM_sqr(sqrt_s(), m1, m2);
}

double ScatterAction::relative_velocity() const {
  const double m1 = incoming_particles()[0].effective_mass();
  const double m2 = incoming_particles()[1].effective_mass();
  const double m_s = mandelstam_s();
  const double lamb = lambda_tilde(m_s, m1 * m1, m2 * m2);
  return std::sqrt(lamb) / (2. * incoming_particles()[0].momentum().x0() *
                            incoming_particles()[1].momentum().x0());
}

double ScatterAction::transverse_distance_sqr() const {
  // local copy of particles (since we need to boost them)
  ParticleData p_a = incoming_particles_[0];
  ParticleData p_b = incoming_particles_[1];
  /* Boost particles to center-of-momentum frame. */
  const ThreeVector velocity = beta_cm();
  p_a.boost(velocity);
  p_b.boost(velocity);
  const ThreeVector pos_diff =
      p_a.position().threevec() - p_b.position().threevec();
  const ThreeVector mom_diff =
      p_a.momentum().threevec() - p_b.momentum().threevec();

  logg[LScatterAction].debug("Particle ", incoming_particles_,
                             " position difference [fm]: ", pos_diff,
                             ", momentum difference [GeV]: ", mom_diff);

  const double dp2 = mom_diff.sqr();
  const double dr2 = pos_diff.sqr();
  /* Zero momentum leads to infite distance. */
  if (dp2 < really_small) {
    return dr2;
  }
  const double dpdr = pos_diff * mom_diff;

  /* UrQMD squared distance criterion, in the center of momentum frame:
   * position of particle a: x_a
   * position of particle b: x_b
   * momentum of particle a: p_a
   * momentum of particle b: p_b
   * d^2_{coll} = (x_a - x_b)^2 - ((x_a - x_b) . (p_a - p_b))^2 / (p_a - p_b)^2
   */
  const double result = dr2 - dpdr * dpdr / dp2;
  return result > 0.0 ? result : 0.0;
}

double ScatterAction::cov_transverse_distance_sqr() const {
  // local copy of particles (since we need to boost them)
  ParticleData p_a = incoming_particles_[0];
  ParticleData p_b = incoming_particles_[1];

  const FourVector delta_x = p_a.position() - p_b.position();
  const double mom_diff_sqr =
      (p_a.momentum().threevec() - p_b.momentum().threevec()).sqr();
  const double x_sqr = delta_x.sqr();

  if (mom_diff_sqr < really_small) {
    return -x_sqr;
  }

  const double p_a_sqr = p_a.momentum().sqr();
  const double p_b_sqr = p_b.momentum().sqr();
  const double p_a_dot_x = p_a.momentum().Dot(delta_x);
  const double p_b_dot_x = p_b.momentum().Dot(delta_x);
  const double p_a_dot_p_b = p_a.momentum().Dot(p_b.momentum());

  const double b_sqr =
      -x_sqr -
      (p_a_sqr * std::pow(p_b_dot_x, 2) + p_b_sqr * std::pow(p_a_dot_x, 2) -
       2 * p_a_dot_p_b * p_a_dot_x * p_b_dot_x) /
          (std::pow(p_a_dot_p_b, 2) - p_a_sqr * p_b_sqr);
  return b_sqr > 0.0 ? b_sqr : 0.0;
}

/**
 * Computes the B coefficients from the STAR fit, see fig. (6) in
 * \iref{STAR:2020phn}.
 *
 * \param[in] plab Lab momentum in GeV.
 *
 * \return B coefficients of high-energy elastic proton-proton scatterings.
 */
static double high_energy_bpp(double plab) {
  double mandelstam_s = s_from_plab(plab, nucleon_mass, nucleon_mass);
  return 7.6 + 0.66 * std::log(mandelstam_s);
}

/**
 * Computes the B coefficients from the Cugnon parametrization of the angular
 * distribution in elastic pp scattering.
 *
 * See equation (8) in \iref{Cugnon:1996kh}.
 * Note: The original Cugnon parametrization is only applicable for
 * plab < 6 GeV and keeps rising above that.
 *
 * \param[in] plab Lab momentum in GeV.
 *
 * \return Cugnon B coefficient for elastic proton-proton scatterings.
 */
static double Cugnon_bpp(double plab) {
  if (plab < 2.) {
    double p8 = pow_int(plab, 8);
    return 5.5 * p8 / (7.7 + p8);
  } else {
    return std::min(high_energy_bpp(plab), 5.334 + 0.67 * (plab - 2.));
  }
}

/**
 * Computes the B coefficients from the Cugnon parametrization of the angular
 * distribution in elastic np scattering.
 *
 * See equation (10) in \iref{Cugnon:1996kh}.
 *
 * \param[in] plab Lab momentum in GeV.
 *
 * \return Cugnon B coefficient for elastic proton-neutron scatterings.
 */
static double Cugnon_bnp(double plab) {
  if (plab < 0.225) {
    return 0.;
  } else if (plab < 0.6) {
    return 16.53 * (plab - 0.225);
  } else if (plab < 1.6) {
    return -1.63 * plab + 7.16;
  } else {
    return Cugnon_bpp(plab);
  }
}

void ScatterAction::sample_angles(std::pair<double, double> masses,
                                  double kinetic_energy_cm) {
  if (is_string_soft_process(process_type_) ||
      (process_type_ == ProcessType::StringHard)) {
    // We potentially have more than two particles, so the following angular
    // distributions don't work. Instead we just keep the angular
    // distributions generated by string fragmentation.
    return;
  }
  assert(outgoing_particles_.size() == 2);

  // NN scattering is anisotropic currently
  const bool nn_scattering = incoming_particles_[0].type().is_nucleon() &&
                             incoming_particles_[1].type().is_nucleon();
  /* Elastic process is anisotropic and
   * the angular distribution is based on the NN elastic scattering. */
  const bool el_scattering = process_type_ == ProcessType::Elastic;

  const double mass_in_a = incoming_particles_[0].effective_mass();
  const double mass_in_b = incoming_particles_[1].effective_mass();

  ParticleData *p_a = &outgoing_particles_[0];
  ParticleData *p_b = &outgoing_particles_[1];

  const double mass_a = masses.first;
  const double mass_b = masses.second;

  const std::array<double, 2> t_range = get_t_range<double>(
      kinetic_energy_cm, mass_in_a, mass_in_b, mass_a, mass_b);
  Angles phitheta;
  if (el_scattering && !isotropic_) {
    /** NN → NN: Choose angular distribution according to Cugnon
     * parametrization,
     * see \iref{Cugnon:1996kh}. */
    double mandelstam_s_new = 0.;
    if (nn_scattering) {
      mandelstam_s_new = mandelstam_s();
    } else {
      /* In the case of elastic collisions other than NN collisions,
       * there is an ambiguity on how to get the lab-frame momentum (plab),
       * since the incoming particles can have different masses.
       * Right now, we first obtain the center-of-mass momentum
       * of the collision (pcom_now).
       * Then, the lab-frame momentum is evaluated from the mandelstam s,
       * which yields the original center-of-mass momentum
       * when nucleon mass is assumed. */
      const double pcm_now = pCM_from_s(mandelstam_s(), mass_in_a, mass_in_b);
      mandelstam_s_new =
          4. * std::sqrt(pcm_now * pcm_now + nucleon_mass * nucleon_mass);
    }
    double bb, a, plab = plab_from_s(mandelstam_s_new);
    if (nn_scattering &&
        p_a->pdgcode().antiparticle_sign() ==
            p_b->pdgcode().antiparticle_sign() &&
        std::abs(p_a->type().charge() + p_b->type().charge()) == 1) {
      // proton-neutron and antiproton-antineutron
      bb = std::max(Cugnon_bnp(plab), really_small);
      a = (plab < 0.8) ? 1. : 0.64 / (plab * plab);
    } else {
      /* all others including pp, nn and AQM elastic processes
       * This is applied for all particle pairs, which are allowed to
       * interact elastically. */
      bb = std::max(Cugnon_bpp(plab), really_small);
      a = 1.;
    }
    double t = random::expo(bb, t_range[0], t_range[1]);
    if (random::canonical() > 1. / (1. + a)) {
      t = t_range[0] + t_range[1] - t;
    }
    // determine scattering angles in center-of-mass frame
    phitheta = Angles(2. * M_PI * random::canonical(),
                      1. - 2. * (t - t_range[0]) / (t_range[1] - t_range[0]));
  } else if (nn_scattering && p_a->pdgcode().is_Delta() &&
             p_b->pdgcode().is_nucleon() &&
             p_a->pdgcode().antiparticle_sign() ==
                 p_b->pdgcode().antiparticle_sign() &&
             !isotropic_) {
    /** NN → NΔ: Sample scattering angles in center-of-mass frame from an
     * anisotropic angular distribution, using the same distribution as for
     * elastic pp scattering, as suggested in \iref{Cugnon:1996kh}. */
    const double plab = plab_from_s(mandelstam_s());
    const double bb = std::max(Cugnon_bpp(plab), really_small);
    double t = random::expo(bb, t_range[0], t_range[1]);
    if (random::canonical() > 0.5) {
      t = t_range[0] + t_range[1] - t;  // symmetrize
    }
    phitheta = Angles(2. * M_PI * random::canonical(),
                      1. - 2. * (t - t_range[0]) / (t_range[1] - t_range[0]));
  } else if (nn_scattering && p_b->pdgcode().is_nucleon() && !isotropic_ &&
             (p_a->type().is_Nstar() || p_a->type().is_Deltastar())) {
    /** NN → NR: Fit to HADES data, see \iref{Agakishiev:2014wqa}. */
    const std::array<double, 4> p{1.46434, 5.80311, -6.89358, 1.94302};
    const double a = p[0] + mass_a * (p[1] + mass_a * (p[2] + mass_a * p[3]));
    /*  If the resonance is so heavy that the index "a" exceeds 30,
     *  the power function turns out to be too sharp. Take t directly to be
     *  t_0 in such a case. */
    double t = t_range[0];
    if (a < 30) {
      t = random::power(-a, t_range[0], t_range[1]);
    }
    if (random::canonical() > 0.5) {
      t = t_range[0] + t_range[1] - t;  // symmetrize
    }
    phitheta = Angles(2. * M_PI * random::canonical(),
                      1. - 2. * (t - t_range[0]) / (t_range[1] - t_range[0]));
  } else {
    /* isotropic angular distribution */
    phitheta.distribute_isotropically();
  }

  ThreeVector pscatt = phitheta.threevec();
  // 3-momentum of first incoming particle in center-of-mass frame
  ThreeVector pcm =
      incoming_particles_[0].momentum().lorentz_boost(beta_cm()).threevec();
  pscatt.rotate_z_axis_to(pcm);

  // final-state CM momentum
  const double p_f = pCM(kinetic_energy_cm, mass_a, mass_b);
  if (!(p_f > 0.0)) {
    logg[LScatterAction].warn("Particle: ", p_a->pdgcode(),
                              " radial momentum: ", p_f);
    logg[LScatterAction].warn("Etot: ", kinetic_energy_cm, " m_a: ", mass_a,
                              " m_b: ", mass_b);
  }
  p_a->set_4momentum(mass_a, pscatt * p_f);
  p_b->set_4momentum(mass_b, -pscatt * p_f);

  /* Debug message is printed before boost, so that p_a and p_b are
   * the momenta in the center of mass frame and thus opposite to
   * each other.*/
  logg[LScatterAction].debug("p_a: ", *p_a, "\np_b: ", *p_b);
}

void ScatterAction::elastic_scattering() {
  // copy initial particles into final state
  outgoing_particles_[0] = incoming_particles_[0];
  outgoing_particles_[1] = incoming_particles_[1];
  // resample momenta
  sample_angles({outgoing_particles_[0].effective_mass(),
                 outgoing_particles_[1].effective_mass()},
                sqrt_s());
}

void ScatterAction::inelastic_scattering() {
  // create new particles
  sample_2body_phasespace();
  assign_formation_time_to_outgoing_particles();
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    assign_unpolarized_spin_vector_to_outgoing_particles();
  }
}

void ScatterAction::two_to_many_scattering() {
  sample_manybody_phasespace();
  assign_formation_time_to_outgoing_particles();
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    assign_unpolarized_spin_vector_to_outgoing_particles();
  }
  logg[LScatterAction].debug("2->", outgoing_particles_.size(),
                             " scattering:", incoming_particles_, " -> ",
                             outgoing_particles_);
}

void ScatterAction::resonance_formation() {
  if (outgoing_particles_.size() != 1) {
    std::string s =
        "resonance_formation: "
        "Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + " (";
    s += incoming_particles_[0].pdgcode().string() + " + ";
    s += incoming_particles_[1].pdgcode().string() + ")";
    throw InvalidResonanceFormation(s);
  }
  // Set the momentum of the formed resonance in its rest frame.
  outgoing_particles_[0].set_4momentum(
      total_momentum_of_outgoing_particles().abs(), 0., 0., 0.);
  assign_formation_time_to_outgoing_particles();
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    assign_unpolarized_spin_vector_to_outgoing_particles();
  }
  /* this momentum is evaluated in the computational frame. */
  logg[LScatterAction].debug("Momentum of the new particle: ",
                             outgoing_particles_[0].momentum());
}

/* This function generates the outgoing state when
 * ScatterAction::string_excitation() is used */

void ScatterAction::create_string_final_state() {
  outgoing_particles_ = string_process_->get_final_state();
  assign_formation_time_to_outgoing_particles();
  /* Check momentum difference for debugging */
  FourVector out_mom;
  for (ParticleData data : outgoing_particles_) {
    out_mom += data.momentum();
  }
  logg[LPythia].debug("Incoming momenta string:", total_momentum());
  logg[LPythia].debug("Outgoing momenta string:", out_mom);
}

/* This function will generate outgoing particles in computational frame
 * from a hard process.
 * The way to excite soft strings is based on the UrQMD model */

void ScatterAction::string_excitation() {
  assert(incoming_particles_.size() == 2);
  // Disable floating point exception trap for Pythia
  {
    DisableFloatTraps guard;
    /* initialize the string_process_ object for this particular collision */
    string_process_->init(incoming_particles_, time_of_execution_);
    /* implement collision */
    bool success = false;
    int ntry = 0;
    const int ntry_max = 10000;
    while (!success && ntry < ntry_max) {
      ntry++;
      switch (process_type_) {
        case ProcessType::StringSoftSingleDiffractiveAX:
          /* single diffractive to A+X */
          success = string_process_->next_SDiff(true);
          break;
        case ProcessType::StringSoftSingleDiffractiveXB:
          /* single diffractive to X+B */
          success = string_process_->next_SDiff(false);
          break;
        case ProcessType::StringSoftDoubleDiffractive:
          /* double diffractive */
          success = string_process_->next_DDiff();
          break;
        case ProcessType::StringSoftNonDiffractive:
          /* soft non-diffractive */
          success = string_process_->next_NDiffSoft();
          break;
        case ProcessType::StringSoftAnnihilation:
          /* soft BBbar 2 mesonic annihilation */
          success = string_process_->next_BBbarAnn();
          break;
        case ProcessType::StringHard:
          success = string_process_->next_NDiffHard();
          break;
        default:
          logg[LPythia].error("Unknown string process required.");
          success = false;
      }
    }
    if (ntry == ntry_max) {
      /* If pythia fails to form a string, it is usually because the energy
       * is not large enough. In this case, annihilation is then enforced. If
       * this process still does not not produce any results, it defaults to
       * an elastic collision. */
      bool success_newtry = false;

      /* Check if the initial state is a baryon-antibaryon state.*/
      PdgCode part1 = incoming_particles_[0].pdgcode(),
              part2 = incoming_particles_[1].pdgcode();
      bool is_BBbar_Pair = (part1.baryon_number() != 0) &&
                           (part1.baryon_number() == -part2.baryon_number());

      /* Decide on the new process .*/
      if (is_BBbar_Pair) {
        process_type_ = ProcessType::StringSoftAnnihilation;
      } else {
        process_type_ = ProcessType::StringSoftDoubleDiffractive;
      }
      /* Perform the new process*/
      int ntry_new = 0;
      while (!success_newtry && ntry_new < ntry_max) {
        ntry_new++;
        if (is_BBbar_Pair) {
          success_newtry = string_process_->next_BBbarAnn();
        } else {
          success_newtry = string_process_->next_DDiff();
        }
      }

      if (success_newtry) {
        create_string_final_state();
      }

      if (!success_newtry) {
        /* If annihilation fails:
         * Particles are normally added after process selection for
         * strings, outgoing_particles is still uninitialized, and memory
         * needs to be allocated. We also shift the process_type_ to elastic
         * so that sample_angles does a proper treatment. */
        outgoing_particles_.reserve(2);
        outgoing_particles_.push_back(ParticleData{incoming_particles_[0]});
        outgoing_particles_.push_back(ParticleData{incoming_particles_[1]});
        process_type_ = ProcessType::FailedString;
        elastic_scattering();
      }
    } else {
      create_string_final_state();
    }
  }
}

static void boost_spin_vectors_after_elastic_scattering(
    ParticleData &outgoing_particle_a, ParticleData &outgoing_particle_b) {
  // Boost spin vectors
  outgoing_particle_a.set_spin_vector(
      outgoing_particle_a.spin_vector().lorentz_boost(
          outgoing_particle_a.velocity()));
  outgoing_particle_b.set_spin_vector(
      outgoing_particle_b.spin_vector().lorentz_boost(
          outgoing_particle_b.velocity()));
}

void ScatterAction::spin_interaction() {
  if (spin_interaction_type_ != SpinInteractionType::Off) {
    /* 2->2 elastic scattering */
    if (process_type_ == ProcessType::Elastic) {
      // Final boost to the outgoing particle momenta
      boost_spin_vectors_after_elastic_scattering(outgoing_particles_[0],
                                                  outgoing_particles_[1]);
    }

    /* 2->1 resonance formation */
    if (process_type_ == ProcessType::TwoToOne) {
      /*
       * @brief Λ+π → Σ* resonance formation with Λ–spin bookkeeping.
       *
       * We do not simulate a direct inelastic Λ+π scattering; instead we form a
       * Σ* resonance and let it decay later. To preserve Λ polarization per
       * arXiv:2404.15890v2, we treat the Σ* spin vector as a proxy for the
       * would-be outgoing Λ spin: at formation, we set the Σ* spin to the
       * incoming Λ spin and apply a possible spin flip according to the
       * Λ–flip/non-flip fractions extracted from the paper. On Σ* → Λ+π decay,
       * the Σ* spin vector is copied to the Λ, thus transporting Λ polarization
       * through the resonance stage.
       */
      // Identify if the outgoing resonance is a Σ*
      if (outgoing_particles_[0].is_sigmastar()) {
        // Check that one of the incoming particles is a Λ and the other a π
        const bool has_lambda = incoming_particles_[0].pdgcode().is_Lambda() ||
                                incoming_particles_[1].pdgcode().is_Lambda();
        const bool has_pion = incoming_particles_[0].is_pion() ||
                              incoming_particles_[1].is_pion();
        if (has_lambda && has_pion) {
          auto &lambda = (incoming_particles_[0].pdgcode().is_Lambda())
                             ? incoming_particles_[0]
                             : incoming_particles_[1];
          auto &sigma_star = outgoing_particles_[0];

          // Perform spin flip with probability of 2/9
          int random_int = random::uniform_int(1, 9);
          FourVector final_spin_vector = lambda.spin_vector();

          if (random_int <= 7) {
            // No spin flip
            final_spin_vector = final_spin_vector.lorentz_boost(
                outgoing_particles_[0].velocity());
            outgoing_particles_[0].set_spin_vector(final_spin_vector);
          } else {
            // Spin flip in Lambda rest frame
            ThreeVector lambda_velocity = lambda.velocity();
            final_spin_vector =
                final_spin_vector.lorentz_boost(lambda_velocity);

            // Flip the spatial spin vector components
            final_spin_vector[1] = -final_spin_vector[1];
            final_spin_vector[2] = -final_spin_vector[2];
            final_spin_vector[3] = -final_spin_vector[3];

            // Boost back to computational frame and to Sigma* frame
            final_spin_vector =
                final_spin_vector.lorentz_boost(-lambda_velocity);
            final_spin_vector =
                final_spin_vector.lorentz_boost(sigma_star.velocity());
            sigma_star.set_spin_vector(final_spin_vector);
          }
        }
      }
    }
  }
}

void ScatterAction::format_debug_output(std::ostream &out) const {
  out << "Scatter of " << incoming_particles_;
  if (outgoing_particles_.empty()) {
    out << " (not performed)";
  } else {
    out << " to " << outgoing_particles_;
  }
}

}  // namespace smash
