/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/scatteractionmulti.h"

#include "smash/crosssections.h"
#include "smash/integrate.h"
#include "smash/logging.h"
#include "smash/parametrizations.h"

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

void ScatterActionMulti::add_possible_reactions(
    double dt, const double gcell_vol,
    const MultiParticleReactionsBitSet incl_multi) {
  // 3 -> m
  if (incoming_particles_.size() == 3) {
    // 3 -> 1
    if (incl_multi[IncludedMultiParticleReactions::Meson_3to1] == 1) {
      if (three_different_pions(incoming_particles_[0], incoming_particles_[1],
                                incoming_particles_[2])) {
        // 3pi -> omega
        const ParticleTypePtr type_omega = ParticleType::try_find(0x223);
        if (type_omega) {
          add_reaction(make_unique<CollisionBranch>(
              *type_omega,
              probability_three_to_one(*type_omega, dt, gcell_vol,
                                       type_omega->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
        // 3pi -> phi
        const ParticleTypePtr type_phi = ParticleType::try_find(0x333);
        if (type_phi) {
          add_reaction(make_unique<CollisionBranch>(
              *type_phi,
              probability_three_to_one(*type_phi, dt, gcell_vol,
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
              probability_three_to_one(
                  *type_eta_prime, dt, gcell_vol,
                  sym_factor_in * type_eta_prime->spin_degeneracy()),
              ProcessType::MultiParticleThreeMesonsToOne));
        }
      }
    }
    // 3 -> 2
    if (incl_multi[IncludedMultiParticleReactions::Deuteron_3to2] == 1) {
      const PdgCode pdg_a = incoming_particles_[0].pdgcode();
      const PdgCode pdg_b = incoming_particles_[1].pdgcode();
      const PdgCode pdg_c = incoming_particles_[2].pdgcode();
      const ParticleTypePtr type_deuteron =
          ParticleType::try_find(PdgCode::from_decimal(pdg::decimal_d));
      const ParticleTypePtr type_anti_deuteron =
          ParticleType::try_find(PdgCode::from_decimal(pdg::decimal_antid));

      const int spin_factor_inc = pdg_a.spin_degeneracy() *
                                  pdg_b.spin_degeneracy() *
                                  pdg_c.spin_degeneracy();

      if (type_deuteron && type_anti_deuteron) {
        // πpn → πd
        if ((pdg_a.is_pion() && pdg_b == pdg::p && pdg_c == pdg::n) ||
            (pdg_a.is_pion() && pdg_b == pdg::n && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b.is_pion() && pdg_c == pdg::n) ||
            (pdg_a == pdg::n && pdg_b.is_pion() && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b == pdg::n && pdg_c.is_pion()) ||
            (pdg_a == pdg::n && pdg_b == pdg::p && pdg_c.is_pion())) {
          // Get type of incoming π
          ParticleList::iterator it = std::find_if(
              incoming_particles_.begin(), incoming_particles_.end(),
              [](ParticleData x) { return x.is_pion(); });
          const ParticleType& type_pi = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_pi.spin_degeneracy(),
                                 type_deuteron->spin_degeneracy());

          add_reaction(make_unique<CollisionBranch>(
              type_pi, *type_deuteron,
              probability_three_to_two(type_pi, *type_deuteron, dt, gcell_vol,
                                       spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }

        // πp̅n̅ → πd̅
        if ((pdg_a.is_pion() && pdg_b == -pdg::p && pdg_c == -pdg::n) ||
            (pdg_a.is_pion() && pdg_b == -pdg::n && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b.is_pion() && pdg_c == -pdg::n) ||
            (pdg_a == -pdg::n && pdg_b.is_pion() && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b == -pdg::n && pdg_c.is_pion()) ||
            (pdg_a == -pdg::n && pdg_b == -pdg::p && pdg_c.is_pion())) {
          // Get type of incoming π
          ParticleList::iterator it = std::find_if(
              incoming_particles_.begin(), incoming_particles_.end(),
              [](ParticleData x) { return x.is_pion(); });
          const ParticleType& type_pi = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_pi.spin_degeneracy(),
                                 type_anti_deuteron->spin_degeneracy());

          add_reaction(make_unique<CollisionBranch>(
              type_pi, *type_anti_deuteron,
              probability_three_to_two(type_pi, *type_anti_deuteron, dt,
                                       gcell_vol, spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }

        // Nnp → Nd, N̅np → N̅d
        if ((pdg_a.is_nucleon() && pdg_b == pdg::p && pdg_c == pdg::n) ||
            (pdg_a.is_nucleon() && pdg_b == pdg::n && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b.is_nucleon() && pdg_c == pdg::n) ||
            (pdg_a == pdg::n && pdg_b.is_nucleon() && pdg_c == pdg::p) ||
            (pdg_a == pdg::p && pdg_b == pdg::n && pdg_c.is_nucleon()) ||
            (pdg_a == pdg::n && pdg_b == pdg::p && pdg_c.is_nucleon())) {
          int symmetry_factor = 1;  // already true for N̅np → N̅d case

          ParticleList::iterator it =
              std::find_if(incoming_particles_.begin(),
                           incoming_particles_.end(), [](ParticleData x) {
                             return x.pdgcode().antiparticle_sign() == -1;
                           });
          if (it == incoming_particles_.end()) {
            /* Meaning no anti-N found by find_if,
             * therefore not N̅np → N̅d, but Nnp → Nd. */
            symmetry_factor = 2;  // for Nnp → Nd (2 factorial)
            // It is already clear here that we have a double of two N
            if (pdg_a == pdg_b) {
              it = incoming_particles_.begin();
            } else {
              // If a and b are not the double, then c has to be part of it
              it = incoming_particles_.begin() + 2;
            }
          }
          const ParticleType& type_N = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_N.spin_degeneracy(),
                                 type_deuteron->spin_degeneracy());

          add_reaction(make_unique<CollisionBranch>(
              type_N, *type_deuteron,
              probability_three_to_two(type_N, *type_deuteron, dt, gcell_vol,
                                       symmetry_factor * spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }

        // Np̅n̅ → Nd̅, N̅p̅n̅ → N̅d̅
        if ((pdg_a.is_nucleon() && pdg_b == -pdg::p && pdg_c == -pdg::n) ||
            (pdg_a.is_nucleon() && pdg_b == -pdg::n && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b.is_nucleon() && pdg_c == -pdg::n) ||
            (pdg_a == -pdg::n && pdg_b.is_nucleon() && pdg_c == -pdg::p) ||
            (pdg_a == -pdg::p && pdg_b == -pdg::n && pdg_c.is_nucleon()) ||
            (pdg_a == -pdg::n && pdg_b == -pdg::p && pdg_c.is_nucleon())) {
          int symmetry_factor = 1;  // already true for Np̅n̅ → Nd̅ case

          ParticleList::iterator it =
              std::find_if(incoming_particles_.begin(),
                           incoming_particles_.end(), [](ParticleData x) {
                             return x.pdgcode().antiparticle_sign() == 1;
                           });
          if (it == incoming_particles_.end()) {
            /* Meaning no N found by find_if,
             * therefore not Np̅n̅ → Nd̅, but N̅p̅n̅ → N̅d̅. */
            symmetry_factor = 2;  // for N̅p̅n̅ → N̅d̅ (2 factorial)
            // It is already clear here that we have a double of two N̅
            if (pdg_a == pdg_b) {
              it = incoming_particles_.begin();
            } else {
              // If a and b are not the double, then c has to be part of it
              it = incoming_particles_.begin() + 2;
            }
          }
          const ParticleType& type_N = it->type();

          const double spin_degn =
              react_degen_factor(spin_factor_inc, type_N.spin_degeneracy(),
                                 type_anti_deuteron->spin_degeneracy());

          add_reaction(make_unique<CollisionBranch>(
              type_N, *type_anti_deuteron,
              probability_three_to_two(type_N, *type_anti_deuteron, dt,
                                       gcell_vol, symmetry_factor * spin_degn),
              ProcessType::MultiParticleThreeToTwo));
        }
      }
    }
  }
  // 5 -> 2
  if (incoming_particles_.size() == 5) {
    // TODO(stdnmr) Introduce config flag for 5-to-2 here
    if (true) {
      // TODO(stdnmr) finalize the inc. pion if statement
      if (all_incoming_particles_are_pions_and_have_charge_zero_together(incoming_particles_[0], incoming_particles_[1],
                                incoming_particles_[2], incoming_particles_[3], incoming_particles_[4])) {
        // TODO(stdnmr) calculate correct symmetry and spin factors
        const double spin_degn = 1.0;
        const double symmetry_factor = 1.0;

        const ParticleTypePtr type_N = ParticleType::try_find(pdg::p);
        const ParticleTypePtr type_anti_N = ParticleType::try_find(-pdg::p);
        if (type_N) {
          // TODO(stdnmr) probably also need to add the neutron reactions here
          add_reaction(make_unique<CollisionBranch>(
              *type_N, *type_anti_N,
              probability_five_to_two(*type_N, dt, gcell_vol,
                                      symmetry_factor * spin_degn),
              ProcessType::MultiParticleFiveToTwo));
        }
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
    case ProcessType::MultiParticleFiveToTwo:
      /* 5->2 scattering */
      five_to_two();
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

double ScatterActionMulti::probability_three_to_one(
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
    const double gcell_vol, const double degen_factor) const {
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

double ScatterActionMulti::probability_five_to_two(
    const ParticleType& type_out, double dt, const double gcell_vol,
    const double degen_factor) const {
  const double e1 = incoming_particles_[0].momentum().x0();
  const double e2 = incoming_particles_[1].momentum().x0();
  const double e3 = incoming_particles_[2].momentum().x0();
  const double e4 = incoming_particles_[3].momentum().x0();
  const double e5 = incoming_particles_[4].momentum().x0();
  const double mout = type_out.mass();

  const double man_s = sqrt_s() * sqrt_s();
  const double lamb = lambda_tilde(man_s, mout * mout, mout * mout);

  // Oscars parametrization for phi5
  const double s_zero = 25 * pion_mass * pion_mass;
  const double fit_a = 5.02560248e-11;
  const double fit_alpha = 1.982;
  const double ph_sp_5 = fit_a * std::pow(man_s - s_zero, 5.0) * std::pow(1 + man_s / s_zero, -fit_alpha);

  // TODO(stdnmr) Clarify if want to account for other baryons than p and if
  // this is the same cross section as for the inverse process
  const double xs = xs_ppbar_annihilation(man_s) / gev2_mb;

  return dt / std::pow(gcell_vol, 4.0) * 1. / (32. * e1 * e2 * e3 * e4 * e5) *
         xs / (4. * M_PI * man_s) * lamb / ph_sp_5 * std::pow(hbarc, 11.0) * degen_factor;
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
  logg[LScatterActionMulti].debug("3->2 scattering:", incoming_particles_,
                                  " -> ", outgoing_particles_);
}

void ScatterActionMulti::five_to_two() {
  sample_2body_phasespace();
  // Make sure to assign formation times before boost to the computational frame
  assign_formation_time_to_outgoing_particles();
  logg[LScatterActionMulti].debug("5->2 scattering:", incoming_particles_,
                                  " -> ", outgoing_particles_);
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

bool ScatterActionMulti::all_incoming_particles_are_pions_and_have_charge_zero_together(const ParticleData& data_a,
                                       const ParticleData& data_b,
                                       const ParticleData& data_c,
                                       const ParticleData& data_d,
                                       const ParticleData& data_e) const {
     const PdgCode pdg_a = data_a.pdgcode();
     const PdgCode pdg_b = data_b.pdgcode();
     const PdgCode pdg_c = data_c.pdgcode();
     const PdgCode pdg_d = data_d.pdgcode();
     const PdgCode pdg_e = data_e.pdgcode();
  return (pdg_a.is_pion() && pdg_b.is_pion() && pdg_c.is_pion() && pdg_d.is_pion() && pdg_e.is_pion()) &&
         (pdg_a.charge() + pdg_b.charge() + pdg_c.charge() + pdg_d.charge() + pdg_e.charge() == 0);

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
