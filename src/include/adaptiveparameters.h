/*
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 */
#ifndef SRC_INCLUDE_ADAPTIVEPARAMETERS_H_
#define SRC_INCLUDE_ADAPTIVEPARAMETERS_H_

#include <utility>

#include "action.h"
#include "actions.h"

namespace Smash {

/**
 * This struct contains all the additional parameters that are needed when
 * using adaptive time steps.
 */
struct AdaptiveParameters {
  /**
   * Calculate the recommended time step size.
   */
  float new_dt(float rate) const {
    return target_missed_actions / rate;
  }

  /**
   * Calculate the rate that would lead to the given time step size.
   */
  float rate_from_dt(float dt) const {
    return target_missed_actions / dt;
  }

  /**
   * Calculate the fraction of missed actions and the allowed deviation.
   *
   * \param actions The actions from the current time step.
   * \param previous_rate The smoothed rate from the previous time step.
   * \param N_particles The number of particles in the system.
   * \return Fraction of missed actions and allowed deviation.
   */
  std::pair<float, float> calc_missed_actions_allowed_deviation(
      const Actions &actions, float previous_rate, size_t N_particles) const {
    const auto &log = logger<LogArea::AdaptiveTS>();
    uint32_t N_inc = 0u;
    for (const auto &a : actions) {
      N_inc += a->incoming_particles().size();
    }
    // factor for the incoming particles
    const float f_inc =
        N_inc == 0u ? 0.f : static_cast<float>(N_inc) / actions.size();
    log.debug("Factor for the incoming particles: ", f_inc);
    const float fr_mi = 0.5f * N_inc * f_inc / N_particles;
    log.debug("Fraction of missed actions: ", fr_mi);
    const float al_dev =
        deviation_factor * previous_rate *
        std::sqrt(0.5f * f_inc / target_missed_actions / N_particles);
    log.debug("Allowed deviation: ", al_dev);
    return std::make_pair(fr_mi, al_dev);
  }

  /**
   * The smoothing factor \f$ \alpha \f$ is responsible for smoothing the
   * estimate of the rate \f$ r \f$. It is used in the following equation
   * \f[ r = r_{old} + \alpha (r_{current} - r_{old}) \f]. An
   * \f$ \alpha \f$ of 1 corresponds to not taking an average at all and 0
   * corresponds to not changing the estimate of the rate at all with the new
   * data.
   *
   * Default value: 0.1
   */
  float smoothing_factor = 0.1f;

  /**
   * The fraction of missed actions that is targeted by the algorithm. A smaller
   * value will lead to smaller time steps but also less missed actions.
   *
   * Default value: 1%
   */
  float target_missed_actions = 0.01f;

  /**
   * The deviation factor \f$ f_D \f$ sets the limit by how much the criterion
   * can be exceeded before the time step is simply aborted. A value of 0
   * corresponds to no room for overshooting at all which will lead to many
   * aborted time steps and consequently longer runtime. The condition is
   * \f[ r - \bar{r} < f_D \sqrt{\bar{r}/(dt \cdot N_P)} \f].
   *
   * Default value: 2.5
   */
  float deviation_factor = 2.5f;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_ADAPTIVEPARAMETERS_H_
