/*
 *
 *    Copyright (c) 2014-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionmulti.h"

#include "smash/crosssections.h"
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
                                                const bool three_to_one,
                                                const bool two_to_three) {
  // 3 -> m
  if (incoming_particles_.size() == 3) {
    // 3 -> 1
    if (three_to_one) {
      if (three_different_pions(incoming_particles_[0], incoming_particles_[1],
                                incoming_particles_[2])) {
        // 3pi -> omega
        const ParticleTypePtr type_omega = ParticleType::try_find(0x223);
        if (type_omega) {
          add_reaction(make_unique<CollisionBranch>(
              *type_omega,
              probability_three_meson_to_one(*type_omega, dt, gcell_vol,
                                             type_omega->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
        // 3pi -> phi
        const ParticleTypePtr type_phi = ParticleType::try_find(0x333);
        if (type_phi) {
          add_reaction(make_unique<CollisionBranch>(
              *type_phi,
              probability_three_meson_to_one(*type_phi, dt, gcell_vol,
                                             type_phi->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
      } else if (two_pions_eta(incoming_particles_[0], incoming_particles_[1],
                               incoming_particles_[2])) {
        // eta2pi -> eta-prime
        const ParticleTypePtr type_eta_prime = ParticleType::try_find(0x331);

        int sym_factor_in = 1;
        if (incoming_particles_[0].type() == incoming_particles_[1].type() ||
            incoming_particles_[1].type() == incoming_particles_[2].type() ||
            incoming_particles_[2].type() == incoming_particles_[0].type()) {
          sym_factor_in = 2;  // 2 factorial
        }

        if (type_eta_prime) {
          add_reaction(make_unique<CollisionBranch>(
              *type_eta_prime,
              probability_three_meson_to_one(
                  *type_eta_prime, dt, gcell_vol,
                  sym_factor_in * type_eta_prime->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
      }
    }
    if (two_to_three) {
      // 3 -> 2
      if (possible_three_to_two_reaction(incoming_particles_[0], incoming_particles_[1],
                                         incoming_particles_[2])) {

        // nppi -> dpi
        const ParticleTypePtr type_out1 = ParticleType::try_find(PdgCode::from_decimal(pdg::decimal_d));

        PdgCode pdg_of_inc_pion;
        for (auto &part : incoming_particles_) {
          if (part.is_pion()) {
            pdg_of_inc_pion = part.pdgcode();
            break;
          }
        }
        const ParticleTypePtr type_out2 = ParticleType::try_find(pdg_of_inc_pion);

        const int degen = 1;  // TODO(stdnmr) think about later

        add_reaction(make_unique<CollisionBranch>(
            *type_out1, *type_out2,
            probability_three_to_two(*type_out1, *type_out2, dt, gcell_vol, degen),
            ProcessType::MultiParticleThreeToTwo));
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
    case ProcessType::MultiParticleThreeToTwo:
      /* 3->2 scattering */
      three_to_two();
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
  const double m1 = incoming_particles_[0].effective_mass();
  const double m2 = incoming_particles_[1].effective_mass();
  const double m3 = incoming_particles_[2].effective_mass();
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
    const ParticleType& type_out, double dt, const double gcell_vol,
    const int degen_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double sqrts = sqrt_s();

  const double gamma_decay = type_out.get_partial_width(
      sqrts, {&incoming_particles_[0].type(), &incoming_particles_[1].type(),
              &incoming_particles_[2].type()});

  const double I_3 = calculate_I3(sqrts);
  const double ph_sp_3 =
      1. / (8 * M_PI * M_PI * M_PI) * 1. / (16 * sqrts * sqrts) * I_3;

  const double spec_f_val = type_out.spectral_function(sqrts);

  return dt / (gcell_vol * gcell_vol) * M_PI / (4. * e1 * e2 * e3) *
         gamma_decay / ph_sp_3 * spec_f_val * std::pow(hbarc, 5.0) *
         degen_factor;
}

double ScatterActionMulti::probability_three_to_two(
    const ParticleType& type_out1, const ParticleType& type_out2, double dt,
    const double gcell_vol, const int degen_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double m4 = type_out1.mass();
  const double m5 = type_out2.mass();

  const double sqrts = sqrt_s();
  const double xs =
      CrossSections::two_to_three_xs(type_out1, type_out2, sqrts) / gev2_mb;
  const double lamb = lambda_tilde(sqrts * sqrts, m4 * m4, m5 * m5);

  const double I_3 = calculate_I3(sqrts);
  const double ph_sp_3 =
      1. / (8 * M_PI * M_PI * M_PI) * 1. / (16 * sqrts * sqrts) * I_3;

  return dt / (gcell_vol * gcell_vol) * 1. / (4. * e1 * e2 * e3) * lamb /
         (ph_sp_3 * 8 * M_PI * sqrts * sqrts) * xs * std::pow(hbarc, 5.0) *
         degen_factor;
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

void ScatterActionMulti::three_to_two() {
  sample_2body_phasespace();
  // Make sure to assign formation times before boost to the computational frame
  assign_formation_time_to_outgoing_particles();
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

bool ScatterActionMulti::possible_three_to_two_reaction(const ParticleData& data_a, const ParticleData& data_b,
                                                        const ParticleData& data_c) const {
  // TODO(stdnmr) Rename function for deuterons

  const PdgCode pdg_a = data_a.pdgcode();
  const PdgCode pdg_b = data_b.pdgcode();
  const PdgCode pdg_c = data_c.pdgcode();

  // We want nppi
  return ((pdg_a.is_pion() && ((pdg_b == pdg::p && pdg_c == pdg::n) || (pdg_b == pdg::n && pdg_c == pdg::p))) ||
          (pdg_b.is_pion() && ((pdg_a == pdg::p && pdg_c == pdg::n) || (pdg_a == pdg::n && pdg_c == pdg::p))) ||
          (pdg_c.is_pion() && ((pdg_b == pdg::p && pdg_a == pdg::n) || (pdg_b == pdg::n && pdg_a == pdg::p))));
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
