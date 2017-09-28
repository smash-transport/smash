/*
 *
 *    Copyright (c) 2017-
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_THERMALIZATIONACTION_H_
#define SRC_INCLUDE_THERMALIZATIONACTION_H_

#include "action.h"
#include "grandcan_thermalizer.h"

namespace Smash {

/**
 * \ingroup action
 * ThermalizationAction implements forced thermalization as an Action class.
 * Particles before thermalization are treated as incoming, after
 * thermalization - as outgoing.
 */
class ThermalizationAction : public Action {
 public:
  ThermalizationAction(const GrandCanThermalizer& gct,
                       double absolute_labframe_time);
  // No need to do anything, because outgoing particles are set in constructor
  void generate_final_state() {}
  double raw_weight_value() const { return 0.0; }
  double partial_weight() const { return 0.0; }
  bool any_particles_thermalized() const {
    return (incoming_particles_.size() > 0);
  }
  void format_debug_output(std::ostream& out) const {
    out << " Thermalization action of " << incoming_particles_.size() << " to "
        << outgoing_particles_.size() << " particles.";
  }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_THERMALIZATIONACTION_H_
