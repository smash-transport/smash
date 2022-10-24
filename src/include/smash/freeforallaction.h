/*
 *
 *    Copyright (c) 2017-
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
 * Action class to create any incoming/outgoing particle combination freely.
 */
class FreeforallAction : public Action {
 public:
  /**
   */
  FreeforallAction(const ParticleList &in_part, const ParticleList &out_part,
                   double absolute_labframe_time)
      : Action(in_part, out_part, absolute_labframe_time, ProcessType::Freeforall) {}
  /// Outgoing particles are set in prinicple in constructor
  void generate_final_state() {
    // Set time to time for arbitrary outgoing particles to time of action
    // TODO(#977) Should the position also scrolled back here?
    for (auto &p : outgoing_particles_) {
      p.set_4position({time_of_execution(), p.position().threevec()});
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
