/*
 *
 *    Copyright (c) 2022-2025
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/fluidizationaction.h"

#include "smash/logging.h"
#include "smash/quantumnumbers.h"

namespace smash {
static constexpr int LFluidization = LogArea::HyperSurfaceCrossing::id;

void FluidizationAction::generate_final_state() {
  logg[LFluidization].debug("Process: Fluidization. ");

  // check that there is only 1 incoming particle
  assert(incoming_particles_.size() == 1);

  if (remove_particle_) {
    // Return empty list because we want to remove the incoming particle
    outgoing_particles_ = {};
  } else {
    outgoing_particles_ = incoming_particles_;
    outgoing_particles_[0].fluidize();
  }
}

double FluidizationAction::check_conservation(
    [[maybe_unused]] const uint32_t id_process) const {
  if (remove_particle_) {
    if (unlikely(outgoing_particles_.size() != 0)) {
      throw std::runtime_error(
          "Particle was not removed successfully in fluidization action.");
    }
  } else {
    QuantumNumbers before(incoming_particles_);
    QuantumNumbers after(outgoing_particles_);
    if (unlikely(before != after)) {
      throw std::runtime_error(
          "Conservation laws not obeyed during fluidization, but they should "
          "since supposedly no removal was done.");
    }
    if (unlikely(outgoing_particles_.size() == 0)) {
      throw std::runtime_error(
          "Particle was removed in a FluidizationNoRemoval process.");
    }
  }

  return 0.;
}

}  // namespace smash
