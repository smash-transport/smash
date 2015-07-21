/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_DECAYACTIONDILEPTON_H_
#define SRC_INCLUDE_DECAYACTIONDILEPTON_H_

#include "decayaction.h"

namespace Smash {

class DecayActionDilepton : public DecayAction {
 public:
   DecayActionDilepton(const ParticleData &p, float time_of_execution,
                        float shining_weight, float dilepton_mass);

  float raw_weight_value() const override {
    return shining_weight_;
  }

  void one_to_three() override;

 private:
   /** The shining weight is a weight you apply to every Dilepton Decay. Because
   * we radiate dileptons at every timestep to increase statistics, we
   * afterwards weight them to correct the dilepton decay yields.
   */
   float shining_weight_;
   // #CleanUp (doc)
   float dilepton_mass_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONDILEPTON_H_