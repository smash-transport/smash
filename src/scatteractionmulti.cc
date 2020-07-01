/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionmulti.h"

#include "smash/integrate.h"
#include "smash/logging.h"

namespace smash {
static constexpr int LScatterActionMulti = LogArea::ScatterActionMulti::id;

ScatterActionMulti::ScatterActionMulti(const ParticleList& in_plist,
                                       double time)
    : Action(in_plist, time), total_probability_(0.) {}

void ScatterActionMulti::add_reaction(CollisionBranchPtr p) {
  add_process<CollisionBranch>(p, reaction_channels_, total_probability_);
}

void ScatterActionMulti::add_reactions(CollisionBranchList pv) {
  add_processes<CollisionBranch>(std::move(pv), reaction_channels_,
                                 total_probability_);
}

double ScatterActionMulti::get_total_weight() const {
  double xsec_scaling = 1.0;
  for (const ParticleData& in_part : incoming_particles_) {
    xsec_scaling *= in_part.xsec_scaling_factor();
  }
  return total_probability_ * xsec_scaling;
}

double ScatterActionMulti::get_partial_weight() const {
  double xsec_scaling = 1.0;
  for (const ParticleData& in_part : incoming_particles_) {
    xsec_scaling *= in_part.xsec_scaling_factor();
  }
  return partial_probability_ * xsec_scaling;
}

void ScatterActionMulti::add_possible_reactions(double dt,
                                                const double gcell_vol,
                                                const bool three_to_one) {
  if (three_to_one && incoming_particles().size() == 3) {
    if (three_different_pions(incoming_particles()[0], incoming_particles()[1],
                              incoming_particles()[2])) {
      // 3pi -> omega
      const ParticleTypePtr type_omega = ParticleType::try_find(0x223);
      if (type_omega) {
        add_reaction(make_unique<CollisionBranch>(
            *type_omega,
            probability_three_meson_to_one(*type_omega, dt, gcell_vol),
            ProcessType::MultiParticleThreeMesonsToOne));
      }
      // 3pi -> phi
      const ParticleTypePtr type_phi = ParticleType::try_find(0x333);
      if (type_phi) {
        add_reaction(make_unique<CollisionBranch>(
            *type_phi, probability_three_meson_to_one(*type_phi, dt, gcell_vol),
            ProcessType::MultiParticleThreeMesonsToOne));
      }
    } else if (two_pions_eta(incoming_particles()[0], incoming_particles()[1],
                             incoming_particles()[2])) {
      // eta2pi -> eta-prime
      const ParticleTypePtr type_eta_prime = ParticleType::try_find(0x331);
      if (type_eta_prime) {
        // TODO(stdnmr) Do we need a symmetry factor if we have two pi0?
        add_reaction(make_unique<CollisionBranch>(
            *type_eta_prime,
            probability_three_meson_to_one(*type_eta_prime, dt, gcell_vol),
            ProcessType::MultiParticleThreeMesonsToOne));
      }
    }
  }
}

void ScatterActionMulti::generate_final_state() {
  logg[LScatterActionMulti].debug("Incoming particles: ", incoming_particles_);

  /* Decide for a particular final state. */
  const CollisionBranch* proc =
      choose_channel<CollisionBranch>(reaction_channels_, total_probability_);
  process_type_ = proc->get_type();
  outgoing_particles_ = proc->particle_list();
  partial_probability_ = proc->weight();

  logg[LScatterActionMulti].debug("Chosen channel: ", process_type_,
                                  outgoing_particles_);

  switch (process_type_) {
    case ProcessType::MultiParticleThreeMesonsToOne:
      /* n->1 annihilation */
      annihilation();
      break;
    default:
      throw InvalidScatterActionMulti(
          "ScatterActionMulti::generate_final_state: Invalid process type " +
          std::to_string(static_cast<int>(process_type_)) + " was requested.");
  }

  /* The production point of the new particles.  */
  FourVector middle_point = get_interaction_point();

  for (ParticleData& new_particle : outgoing_particles_) {
    // Boost to the computational frame
    new_particle.boost_momentum(
        -total_momentum_of_outgoing_particles().velocity());
    /* Set positions of the outgoing particles */
    new_particle.set_4position(middle_point);
  }
}

double ScatterActionMulti::calculate_I3(const double sqrts) const {
  static Integrator integrate;
  const double m1 = incoming_particles()[0].effective_mass();
  const double m2 = incoming_particles()[1].effective_mass();
  const double m3 = incoming_particles()[2].effective_mass();
  const double lower_bound = (m1 + m2) * (m1 + m2);
  const double upper_bound = (sqrts - m3) * (sqrts - m3);
  const auto result = integrate(lower_bound, upper_bound, [&](double m12_sqr) {
    const double m12 = std::sqrt(m12_sqr);
    const double e2_star = (m12_sqr - m1 * m1 + m2 * m2) / (2 * m12);
    const double e3_star = (sqrts * sqrts - m12_sqr - m3 * m3) / (2 * m12);
    const double m23_sqr_min =
        (e2_star + e3_star) * (e2_star + e3_star) -
        std::pow(std::sqrt(e2_star * e2_star - m2 * m2) +
                     std::sqrt(e3_star * e3_star - m3 * m3),
                 2.0);
    const double m23_sqr_max =
        (e2_star + e3_star) * (e2_star + e3_star) -
        std::pow(std::sqrt(e2_star * e2_star - m2 * m2) -
                     std::sqrt(e3_star * e3_star - m3 * m3),
                 2.0);
    return m23_sqr_max - m23_sqr_min;
  });

  return result;
}

double ScatterActionMulti::probability_three_meson_to_one(
    const ParticleType& type_out, double dt, const double gcell_vol) const {
  const double e1 = incoming_particles()[0].momentum().x0();
  const double e2 = incoming_particles()[1].momentum().x0();
  const double e3 = incoming_particles()[2].momentum().x0();
  const double sqrts = sqrt_s();

  const double gamma_decay = type_out.get_partial_width(
      sqrts, {&incoming_particles()[0].type(), &incoming_particles()[1].type(),
              &incoming_particles()[2].type()});

  // Spin degneracy of outgoing particles (incoming p. assumed to have no spin)
  const int spin_deg_out = type_out.spin_degeneracy();
  const double I_3 = calculate_I3(sqrts);
  const double ph_sp_3 =
      1. / (8 * M_PI * M_PI * M_PI) * 1. / (16 * sqrts * sqrts) * I_3;

  const double spec_f_val = type_out.spectral_function(sqrts);

  // Symmetry factor for incoming particles
  int sym_factor_in = 1;
  if (incoming_particles()[0].type() == incoming_particles()[1].type() &&
      incoming_particles()[1].type() == incoming_particles()[2].type()) {
    sym_factor_in = 6;  // 3!
  } else if (incoming_particles()[0].type() == incoming_particles()[1].type() ||
             incoming_particles()[1].type() == incoming_particles()[2].type() ||
             incoming_particles()[2].type() == incoming_particles()[0].type()) {
    sym_factor_in = 2;  // 2!
  }

  return dt / (gcell_vol * gcell_vol) * M_PI / (4. * e1 * e2 * e3) *
         gamma_decay / ph_sp_3 * spec_f_val * std::pow(hbarc, 5.0) *
         spin_deg_out * sym_factor_in;
}

void ScatterActionMulti::annihilation() {
  if (outgoing_particles_.size() != 1) {
    std::string s =
        "Annihilation: "
        "Incorrect number of particles in final state: ";
    s += std::to_string(outgoing_particles_.size()) + ".";
    throw InvalidScatterActionMulti(s);
  }
  // Set the momentum of the formed particle in its rest frame.
  outgoing_particles_[0].set_4momentum(
      total_momentum_of_outgoing_particles().abs(), 0., 0., 0.);
  // Make sure to assign formation times before boost to the computational frame
  assign_formation_time_to_outgoing_particles();

  logg[LScatterActionMulti].debug("Momentum of the new particle: ",
                                  outgoing_particles_[0].momentum());
}

bool ScatterActionMulti::three_different_pions(
    const ParticleData& data_a, const ParticleData& data_b,
    const ParticleData& data_c) const {
  // We want a combination of pi+, pi- and pi0
  const PdgCode pdg_a = data_a.pdgcode();
  const PdgCode pdg_b = data_b.pdgcode();
  const PdgCode pdg_c = data_c.pdgcode();

  return (pdg_a.is_pion() && pdg_b.is_pion() && pdg_c.is_pion()) &&
         (pdg_a != pdg_b && pdg_b != pdg_c && pdg_c != pdg_a);
}

bool ScatterActionMulti::two_pions_eta(const ParticleData& data_a,
                                       const ParticleData& data_b,
                                       const ParticleData& data_c) const {
  // We want a combination of pi0, pi0 and eta or pi+, pi- and eta
  const PdgCode pdg_a = data_a.pdgcode();
  const PdgCode pdg_b = data_b.pdgcode();
  const PdgCode pdg_c = data_c.pdgcode();

  return (pdg_a == pdg::pi_z && pdg_b == pdg::pi_z && pdg_c == pdg::eta) ||
         (pdg_a == pdg::pi_z && pdg_b == pdg::eta && pdg_c == pdg::pi_z) ||
         (pdg_a == pdg::eta && pdg_b == pdg::pi_z && pdg_c == pdg::pi_z) ||
         (pdg_a == pdg::eta && pdg_b == pdg::pi_m && pdg_c == pdg::pi_p) ||
         (pdg_a == pdg::eta && pdg_b == pdg::pi_p && pdg_c == pdg::pi_m) ||
         (pdg_a == pdg::pi_m && pdg_b == pdg::pi_p && pdg_c == pdg::eta) ||
         (pdg_a == pdg::pi_m && pdg_b == pdg::eta && pdg_c == pdg::pi_p) ||
         (pdg_a == pdg::pi_p && pdg_b == pdg::pi_m && pdg_c == pdg::eta) ||
         (pdg_a == pdg::pi_p && pdg_b == pdg::eta && pdg_c == pdg::pi_m);
}

void ScatterActionMulti::format_debug_output(std::ostream& out) const {
  out << "MultiParticleScatter of " << incoming_particles_;
  if (outgoing_particles_.empty()) {
    out << " (not performed)";
  } else {
    out << " to " << outgoing_particles_;
  }
}

}  // namespace smash
