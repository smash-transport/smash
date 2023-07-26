/*
 *
 *    Copyright (c) 2023
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_FLUIDIZATIONACTION_H_
#define SRC_INCLUDE_SMASH_FLUIDIZATIONACTION_H_

#include "action.h"

namespace smash {

/**
 * \ingroup action
 * FluidizationAction is a special action indicating that a particle will be
 * removed from the hadronic evolution, and considered a fluid to be evolved
 * with an external hydrodynamics model. This can be done with a fixed iso-tau
 * prescription applicable to high beam energy collisions, or dynamically
 * according to the local energy density and using the fluidized particles as
 * sources, which is suitable for low beam energies. */
class FluidizationAction : public Action {
 public:
  /**
   * Construct action of fluidization.
   * \param[in] in_part Incoming particle which surpass the energy density.
   * \param[in] out_part Same particle as above, but propagated to the point of
   * crossing the hypersurface, in case of an iso-tau condition.
   * \param[in] time_until How long until the action is performed \unit{in fm}.
   */
  FluidizationAction(const ParticleData &in_part,
                             const ParticleData &out_part)
      : Action(in_part, out_part, time_until,
               ProcessType::HyperSurfaceCrossing) {}
  double get_total_weight() const override { return 0.0; };
  double get_partial_weight() const override { return 0.0; };
  void format_debug_output(std::ostream &out) const override {
    out << "Fluidization of " << incoming_particles_;
  }

  /**
   * Generate the final state of particles to be fluidized and removes them from
   * the evolution.
   */
  void generate_final_state() override;

  /**
   * Conservation laws should not be obeyed, since particles are being removed.
   *
   * \param[in] id_process process id only used for debugging output.
   */
  double check_conservation(uint32_t id_process) const override;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_
