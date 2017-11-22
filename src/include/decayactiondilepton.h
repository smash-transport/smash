/*
 *
 *    Copyright (c) 2015-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONDILEPTON_H_
#define SRC_INCLUDE_DECAYACTIONDILEPTON_H_

#include "decayaction.h"

namespace smash {

class DecayActionDilepton : public DecayAction {
 public:
  DecayActionDilepton(const ParticleData &p, double time_of_execution,
                      double shining_weight);

  double raw_weight_value() const override { return shining_weight_; }
  double partial_weight() const override {
    return shining_weight_ * branching_;
  }

  void one_to_three() override;

 private:
  /**
   * The shining weight is a weight you apply to every Dilepton Decay. Because
   * we radiate dileptons at every timestep to increase statistics, we
   * afterwards weight them to correct the dilepton decay yields.
   */
  const double shining_weight_;
  /**
   * An additional branching factor that is multiplied with the shining weight.
   * For Dalitz decays, the primary shining weight is based on the integrated
   * width for the channel, and the branching factor corrects for the
   * differential width (evaluated at a particular dilepton mass), relative
   * to the integrated width. It is determined after the dilepton mass is fixed.
   * For direct (2-body) decays, the branching factor equals one.
   */
  double branching_ = 1.;
};

}  // namespace smash

#endif  // SRC_INCLUDE_DECAYACTIONDILEPTON_H_
