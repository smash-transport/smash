/*
 *
 *    Copyright (c) 2014-2020
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/decayaction.h"

#include "smash/decaymodes.h"
#include "smash/logging.h"
#include "smash/pdgcode.h"

namespace smash {
static constexpr int LDecayModes = LogArea::DecayModes::id;

DecayAction::DecayAction(const ParticleData &p, double time)
    : Action({p}, time), total_width_(0.) {}

void DecayAction::add_decays(DecayBranchList pv) {
  add_processes<DecayBranch>(std::move(pv), decay_channels_, total_width_);
}

void DecayAction::add_decay(DecayBranchPtr p) {
  add_process<DecayBranch>(p, decay_channels_, total_width_);
}

void DecayAction::generate_final_state() {
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
      break;
    case 3:
      sample_manybody_phasespace();
      break;
    default:
      throw InvalidDecay(
          "DecayAction::perform: Only 1->2 or 1->3 processes are supported. "
          "Decay from 1->" +
          std::to_string(outgoing_particles_.size()) +
          " was requested. (PDGcode=" +
          incoming_particles_[0].pdgcode().string() + ", mass=" +
          std::to_string(incoming_particles_[0].effective_mass()) + ")");
  }

  // Set formation time.
  for (auto &p : outgoing_particles_) {
    logg[LDecayModes].debug("particle momenta in lrf ", p);
    // assuming decaying particles are always fully formed
    p.set_formation_time(time_of_execution_);
    // Boost to the computational frame
    p.boost_momentum(-total_momentum_of_outgoing_particles().velocity());
    logg[LDecayModes].debug("particle momenta in comp ", p);
  }
}

/* This is overridden from the Action class in order to
 * take care of the angular momentum L_. */
std::pair<double, double> DecayAction::sample_masses(
    double kinetic_energy_cm) const {
  const ParticleType &t_a = outgoing_particles_[0].type();
  const ParticleType &t_b = outgoing_particles_[1].type();

  // start with pole masses
  std::pair<double, double> masses = {t_a.mass(), t_b.mass()};

  if (kinetic_energy_cm < t_a.min_mass_kinematic() + t_b.min_mass_kinematic()) {
    const std::string reaction =
        incoming_particles_[0].type().name() + "→" + t_a.name() + t_b.name();
    throw InvalidResonanceFormation(
        reaction + ": not enough energy, " + std::to_string(kinetic_energy_cm) +
        " < " + std::to_string(t_a.min_mass_kinematic()) + " + " +
        std::to_string(t_b.min_mass_kinematic()));
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
