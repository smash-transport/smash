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
 * FluidizationAction is a special action indicating that a particle will be removed from the hadronic evolution, and considered a fluid to be evolved with an external hydrodynamics model. This can be done with a fixed iso-tau prescription applicable to high beam energy collisions, or dynamically according to the local energy density and using the fluidized particles as sources, which is suitable for low beam energies. */
class FluidizationAction : public Action {
 public:
  /**
   * Construct action of fluidization.
   * \param[in] in_part Data of incoming particle.
   * \param[in] out_part Data of particles which surpass the energy density threshold
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
   * Generate the final state of particles to be fluidized and removes them from the
   * evolution.
   */
  void generate_final_state() override;

  /** 
   * This function overrides Action::check_conservation that returns
   * the amount of energy density violation due to Pythia processes,
   * which is 0. here.
   *
   * \param[in] id_process process id only used for debugging output.
   * \return 0.
   */
  double check_conservation(uint32_t id_process) const override;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_
