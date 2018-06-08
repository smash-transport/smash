/*
 *    Copyright (c) 2015-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */

#ifndef SRC_INCLUDE_ADAPTIVEPARAMETERS_H_
#define SRC_INCLUDE_ADAPTIVEPARAMETERS_H_

#include <utility>

#include "action.h"
#include "actions.h"
#include "configuration.h"

namespace smash {

/**
 * This class implements the updating of adaptive timestep based on the
 * estimate scattering rate. If the rate is too high for a given timestep then
 * the timestep is decreased and vice vesa.
 */
class AdaptiveParameters {
 public:
  /**
   * Construct parameters for adaptive timesteps.
   * \param[in] conf Configuration build from config.yaml.
   * \param[in] dt Time step [fm]
   */
  explicit AdaptiveParameters(Configuration conf, double dt)
      : smoothing_factor_(conf.take({"Smoothing_Factor"}, 0.1)),
        target_missed_actions_(conf.take({"Target_Missed_Actions"}, 0.01)),
        deviation_factor_(conf.take({"Allowed_Deviation"}, 2.5)) {
    initialize(dt);
  }

  /**
   * Calculate the rate that would lead to the given time step size.
   * \param[in] dt time step size [fm]
   */
  void initialize(double dt) { rate_ = target_missed_actions_ / dt; }

  /**
   * Updates timestep if necessary.
   *
   * \param[in] actions The actions from the current time step.
   * \param[in] N_particles The number of particles in the system.
   * \param[out] dt timestep, which can be updated [dt]
   * \return if timestep was changed or not.
   */
  bool update_timestep(const Actions &actions, size_t N_particles, double *dt) {
    const auto &log = logger<LogArea::AdaptiveTS>();
    bool changed_timestep = false;

    uint32_t N_inc = 0u;
    for (const auto &a : actions) {
      N_inc += a->incoming_particles().size();
    }
    // factor for the incoming particles
    const double f_inc =
        (N_inc == 0u) ? 0. : static_cast<double>(N_inc) / actions.size();
    const double fraction_missed = 0.5 * N_inc * f_inc / N_particles;
    const double allowed_deviation =
        deviation_factor_ * rate_ *
        std::sqrt(0.5 * f_inc / target_missed_actions_ / N_particles);
    const double current_rate = fraction_missed / (*dt);
    const double rate_deviation = current_rate - rate_;
    log.debug("Factor for the incoming particles: ", f_inc,
              ", fraction of missed actions: ", fraction_missed,
              ", rate estimate: ", rate_, ", rate deviation = ", rate_deviation,
              ", allowed = ", allowed_deviation);
    if (rate_deviation > allowed_deviation) {
      changed_timestep = true;
      *dt = target_missed_actions_ / rate_;
    }
    // update the estimate of the rate
    rate_ += smoothing_factor_ * rate_deviation;
    return changed_timestep;
  }

  /// Get estimate of current scattering rate
  double rate() const { return rate_; }

  /**
   * The smoothing factor \f$ \alpha \f$ is responsible for smoothing the
   * estimate of the rate \f$ r \f$. It is used in the following equation
   * \f[ r = r_{old} + \alpha (r_{current} - r_{old}) \f]. An
   * \f$ \alpha \f$ of 1 corresponds to not taking an average at all and 0
   * corresponds to not changing the estimate of the rate at all with the new
   * data.
   */
  const double smoothing_factor_;

  /**
   * The fraction of missed actions that is targeted by the algorithm. A smaller
   * value will lead to smaller time steps but also less missed actions.
   */
  const double target_missed_actions_;

  /**
   * The deviation factor \f$ f_D \f$ sets the limit by how much the criterion
   * can be exceeded before the time step is simply aborted. A value of 0
   * corresponds to no room for overshooting at all which will lead to many
   * aborted time steps and consequently longer runtime. The condition is
   * \f[ r - \bar{r} < f_D \sqrt{\bar{r}/(dt \cdot N_P)} \f].
   */
  const double deviation_factor_;

 private:
  /// Estimate of current scattering rate
  double rate_;
};

/**
 * Structured output for the adaptive parameters.
 * \param[in] o Outputstream
 * \param[in] a The set of adaptive parameters.
 */
inline std::ostream &operator<<(std::ostream &o, const AdaptiveParameters &a) {
  return o << "Adaptive time step:\n"
           << "  Smoothing factor: " << a.smoothing_factor_ << "\n"
           << "  Target missed actions: " << 100 * a.target_missed_actions_
           << "%\n"
           << "  Allowed deviation: " << a.deviation_factor_ << "\n";
}

}  // namespace smash

#endif  // SRC_INCLUDE_ADAPTIVEPARAMETERS_H_
