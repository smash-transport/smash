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
                        float shining_weight_);

   // do not perform any dilepton actions = leave it empty
   void perform(Particles *, size_t &) override {};

   /* generate_final_state stayed the same for a first version
      but now for the dalitz decays one have to change it for
      new three body decay behaviours
    */
    void generate_final_state() override;

    float raw_weight_value() const override {
      return shining_weight_;
    }

 private:
   float shining_weight_;

   void dalitz_cinematics();

   void one_to_two() {
     /* Sample the masses and momenta. */
     sample_cms_momenta();
   }
};

}  // namespace Smash

#endif  // SRC_INCLUDE_DECAYACTIONDILEPTON_H_
