/*
 *
 *    Copyright (c) 2023
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

  // Return empty list because we want to remove the incoming particle
  outgoing_particles_ = {};
}

double FluidizationAction::check_conservation(const uint32_t id_process) const {
  QuantumNumbers before(incoming_particles_);
  QuantumNumbers after(outgoing_particles_);
  if (before == after) {
    throw std::runtime_error(
        "Conservation laws obeyed during fluidization, which should not happen "
        "as particles are removed. Particle was not properly removed in "
        "process: " +
        std::to_string(id_process));
  }

  if (outgoing_particles_.size() != 0) {
    throw std::runtime_error(
        "Particle was not removed successfully in fluidization action.");
  }

  return 0.;
}

}  // namespace smash
