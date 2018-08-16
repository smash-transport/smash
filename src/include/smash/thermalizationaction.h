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

namespace smash {

/**
 * \ingroup action
 * ThermalizationAction implements forced thermalization as an Action class.
 * Particles before thermalization are treated as incoming, after
 * thermalization as outgoing. This is an N->M action.
 */
class ThermalizationAction : public Action {
 public:
  /**
   * The inherited class
   * \param[in] gct The thermalization object taking care of removing and
   * sampling new particles \see GrandCanThermalizer.
   * \param[in] absolute_labframe_time Current time in the computational frame.
   */
  ThermalizationAction(const GrandCanThermalizer& gct,
                       double absolute_labframe_time);
  /// No need to do anything, because outgoing particles are set in constructor
  void generate_final_state() {}
  double get_total_weight() const { return 0.0; }
  double get_partial_weight() const { return 0.0; }
  /// This method checks, if there are particles in the region to be thermalized
  bool any_particles_thermalized() const {
    return (incoming_particles_.size() > 0);
  }
  /**
   * Function for debug output of incoming and outgoing particles from
   * thermalization action
   * \param[in] out Location of the output stream
   */
  void format_debug_output(std::ostream& out) const {
    out << " Thermalization action of " << incoming_particles_.size() << " to "
        << outgoing_particles_.size() << " particles.";
  }
};

}  // namespace smash

#endif  // SRC_INCLUDE_THERMALIZATIONACTION_H_
