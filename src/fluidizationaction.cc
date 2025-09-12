/*
 *
 *    Copyright (c) 2022-2024
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
  if (unlikely(outgoing_particles_.size() != 0)) {
    throw std::runtime_error(
        "Particle was not removed successfully in fluidization action.");
  }

  return 0.;
}

}  // namespace smash
