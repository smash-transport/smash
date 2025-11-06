/*
 *
 *    Copyright (c) 2013-2021
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decayaction.h"

#include "smash/decaymodes.h"
#include "smash/logging.h"
#include "smash/pdgcode.h"
#include "smash/potential_globals.h"

namespace smash {

static constexpr int LDecayModes = LogArea::DecayModes::id;

DecayAction::DecayAction(const ParticleData &p, double time,
                         SpinInteractionType spin_interaction_type)
    : Action({p}, time),
      total_width_(0.),
      spin_interaction_type_(spin_interaction_type) {}

void DecayAction::add_decays(DecayBranchList pv) {
  add_processes<DecayBranch>(std::move(pv), decay_channels_, total_width_);
}

void DecayAction::add_decay(DecayBranchPtr p) {
  add_process<DecayBranch>(p, decay_channels_, total_width_);
}

bool DecayAction::sample_outgoing_particles() {
  logg[LDecayModes].debug("Process: Resonance decay. ");
  /* Execute a decay process for the selected particle.
   *
   * randomly select one of the decay modes of the particle
   * according to their relative weights. Then decay the particle
   * by calling function sample_2body_phasespace or sample_manybody_phasespace.
   */
  const DecayBranch *proc =
      choose_channel<DecayBranch>(decay_channels_, total_width_);
  outgoing_particles_ = proc->particle_list();
  // set positions of the outgoing particles
  for (auto &p : outgoing_particles_) {
    p.set_4position(incoming_particles_[0].position());
  }
  process_type_ = proc->get_type();
  L_ = proc->angular_momentum();
  partial_width_ = proc->weight();

  switch (outgoing_particles_.size()) {
    case 2:
      sample_2body_phasespace();
      if (pot_pointer) {
        return was_2body_phase_space_sampled_with_potentials_as_valid_.value();
      } else {
        return true;
      }
    case 3:
      sample_manybody_phasespace();
      return true;
    default:
      throw InvalidDecay(
          "DecayAction::perform: Only 1->2 or 1->3 processes are supported. "
          "Decay from 1->" +
          std::to_string(outgoing_particles_.size()) +
          " was requested. (PDGcode=" +
          incoming_particles_[0].pdgcode().string() + ", mass=" +
          std::to_string(incoming_particles_[0].effective_mass()) + ")");
  }
}

void DecayAction::generate_final_state() {
  int n_try = 1000;
  while (n_try--) {
    if (sample_outgoing_particles())
      break;
  }

  const bool core_in_incoming =
      std::any_of(incoming_particles_.begin(), incoming_particles_.end(),
                  [](const ParticleData &p) { return p.is_core(); });
  // Set formation time.
  for (auto &p : outgoing_particles_) {
    logg[LDecayModes].debug("particle momenta in lrf ", p);
    // assuming decaying particles are always fully formed
    p.set_formation_time(time_of_execution_);
    // Boost to the computational frame
    p.boost_momentum(-total_momentum_of_outgoing_particles().velocity());
    logg[LDecayModes].debug("particle momenta in comp ", p);
    if (core_in_incoming) {
      p.fluidize();
    }
  }

  /*
   * @brief Σ* → Λ + π decay: propagate Λ polarization from the intermediate
   * resonance.
   *
   * During Λ+π → Σ* formation we stored the incoming Λ polarization by
   * writing it into the Σ* spin 4-vector (optionally applying a Λ spin-flip
   * probability, cf. arXiv:2404.15890v2). At decay, we must hand this
   * polarization back to the outgoing Λ to transport Λ polarization through the
   * resonance stage.
   */
  if (spin_interaction_type_ != SpinInteractionType::Off &&
      outgoing_particles_.size() == 2) {
    // Check for Σ* → Λ + π decay channel
    int lambda_idx = -1;
    int pion_idx = -1;
    const bool is_sigmastar_decay = is_sigmastar_to_lambda_pion_decay(
        incoming_particles_[0], outgoing_particles_[0], outgoing_particles_[1],
        lambda_idx, pion_idx);

    if (is_sigmastar_decay) {
      auto &sigma_star = incoming_particles_[0];
      auto &lambda = outgoing_particles_[lambda_idx];
      auto &pion = outgoing_particles_[pion_idx];
      // Copy spin vector from Σ* to Λ and boost it to the Λ frame
      FourVector final_spin_vector =
          sigma_star.spin_vector().lorentz_boost(lambda.velocity());
      lambda.set_spin_vector(final_spin_vector);
      pion.set_spin_vector(FourVector{0., 0., 0., 0.});
    } else {
      // Set unpolarized spin vectors
      assign_unpolarized_spin_vector_to_outgoing_particles();
    }
  } else if (spin_interaction_type_ != SpinInteractionType::Off &&
             outgoing_particles_.size() != 2) {
    // Set unpolarized spin vectors
    assign_unpolarized_spin_vector_to_outgoing_particles();
  }
}

void DecayAction::sample_2body_phasespace() {
  assert(outgoing_particles_.size() == 2);
  const FourVector p_tot = total_momentum_of_outgoing_particles();
  const double cm_kin_energy = p_tot.abs();
  const std::pair<double, double> masses = sample_masses(cm_kin_energy);

  const bool is_valid = !std::isnan(masses.first) && !std::isnan(masses.second);

  if (pot_pointer) {
    was_2body_phase_space_sampled_with_potentials_as_valid_ = is_valid;
  }

  sample_angles(masses, cm_kin_energy);
}

/* This is overridden from the Action class in order to
 * take care of the angular momentum L_. */
std::pair<double, double> DecayAction::sample_masses(
    double kinetic_energy_cm) const {
  const ParticleType &t_a = outgoing_particles_[0].type();
  const ParticleType &t_b = outgoing_particles_[1].type();

  // start with pole masses
  std::pair<double, double> masses = {t_a.mass(), t_b.mass()};

  const bool below_threshold_energy =
      kinetic_energy_cm < t_a.min_mass_kinematic() + t_b.min_mass_kinematic();

  const bool return_nan_on_failure = pot_pointer != nullptr;

  if (below_threshold_energy) {
    if (return_nan_on_failure) {
      return {smash_NaN<double>, smash_NaN<double>};
    } else {
      const std::string reaction =
          incoming_particles_[0].type().name() + "→" + t_a.name() + t_b.name();
      throw InvalidResonanceFormation(
          reaction + ": not enough energy, " +
          std::to_string(kinetic_energy_cm) + " < " +
          std::to_string(t_a.min_mass_kinematic()) + " + " +
          std::to_string(t_b.min_mass_kinematic()));
    }
  }
  // If one of the particles is a resonance, sample its mass.
  if (!t_a.is_stable() && t_b.is_stable()) {
    masses.first = t_a.sample_resonance_mass(t_b.mass(), kinetic_energy_cm, L_);
  } else if (!t_b.is_stable() && t_a.is_stable()) {
    masses.second =
        t_b.sample_resonance_mass(t_a.mass(), kinetic_energy_cm, L_);
  } else if (!t_a.is_stable() && !t_b.is_stable()) {
    // two resonances in final state
    masses = t_a.sample_resonance_masses(t_b, kinetic_energy_cm, L_);
  }

  return masses;
}

void DecayAction::format_debug_output(std::ostream &out) const {
  out << "Decay of " << incoming_particles_ << " to " << outgoing_particles_
      << ", sqrt(s)=" << format(sqrt_s(), "GeV", 11, 9);
}

}  // namespace smash
