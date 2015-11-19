/*
 *
 *    Copyright (c) 2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_WALLCROSSINGACTION_H_
#define SRC_INCLUDE_WALLCROSSINGACTION_H_

#include "action.h"

namespace Smash {


/**
 * \ingroup action
 * WallcrossingAction is a special action which indicates that a particle
 * has crossed a box wall.
 */
class WallcrossingAction : public Action {
 public:
  WallcrossingAction(const ParticleData &in_part, const ParticleData &out_part)
                      : Action(in_part, out_part, 0., ProcessType::Wall ) {}
  float raw_weight_value() const override { return 1; };
  void generate_final_state() override {};
  double sqrt_s() const override {
    return incoming_particles_[0].momentum().abs();
  }
  void format_debug_output(std::ostream &out) const override {
    out << "Wall crossing of " << incoming_particles_;
  }
};


}  // namespace Smash

#endif  // SRC_INCLUDE_WALLCROSSINGACTION_H_
