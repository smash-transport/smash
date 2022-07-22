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
  /// No need to do anything, because outgoing particles are set in constructor
  void generate_final_state() {}
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
