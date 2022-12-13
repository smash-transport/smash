/*
 *
 *    Copyright (c) 2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_FREEFORALLACTION_H_
#define SRC_INCLUDE_SMASH_FREEFORALLACTION_H_

#include "action.h"

namespace smash {

/**
 * \ingroup action
 * \brief Action class to create any incoming/outgoing particle combination freely.
 * This class is in particular designed to add and remove particles from the 
 * evolution. This introduces violations of the conservation laws, but it is
 * e.g. needed for a concurrent running of transport and hydrodynamics.
 */
class FreeforallAction : public Action {
 public:
  /**
   * The FreeforallAction is able to add particles by providing an empty
   * particle list for in_part and a list with particles which are supposed to
   * be added to the evolution as out_part.
   * If particles should be removed a filled particle list with the
   * corresponding particles should be provided as in_part and an empty list as
   * out_part.
   * \param[in] in_part List of incoming particles
   * \param[in] out_part List of outgoing particles
   * \param[in] absolute_labframe_time Absolute time at which the action is
   *                                             supposed to take place
   */
  FreeforallAction(const ParticleList &in_part, const ParticleList &out_part,
                   double absolute_labframe_time)
      : Action(in_part, out_part, absolute_labframe_time,
               ProcessType::Freeforall) {}

  /// Outgoing particles are set in prinicple in constructor
  void generate_final_state() {
    // Set time to time for arbitrary outgoing particles to time of action
    // TODO(#977) Should the position also scrolled back here?
    for (auto &r : outgoing_particles_) {
      r.set_4position({time_of_execution(), r.position().threevec()});
    }
  }

  double get_total_weight() const { return 0.0; }

  double get_partial_weight() const { return 0.0; }

  /**
   * Function for debug output of incoming and outgoing particles from
   * freeforall action
   * \param[in] out Location of the output stream
   */
  void format_debug_output(std::ostream &out) const {
    out << "Freeforall action of " << incoming_particles_.size() << " to "
        << outgoing_particles_.size() << " particles.";
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_FREEFORALLACTION_H_
