/*
 *
 *    Copyright (c) 2016-2017
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_WALLCROSSINGACTION_H_
#define SRC_INCLUDE_WALLCROSSINGACTION_H_

#include "action.h"
#include "actionfinderfactory.h"

namespace Smash {


/**
 * \ingroup action
 * WallcrossingAction is a special action which indicates that a particle
 * has crossed a box wall.
 */
class WallcrossingAction : public Action {
 public:
  WallcrossingAction(const ParticleData &in_part, const ParticleData &out_part,
                     const double time_until = 0.0)
                 : Action(in_part, out_part, time_until, ProcessType::Wall) {}
  double raw_weight_value() const override { return 1; };
  void generate_final_state() override {};
  double sqrt_s() const override {
    return incoming_particles_[0].momentum().abs();
  }
  void format_debug_output(std::ostream &out) const override {
    out << "Wall crossing of " << incoming_particles_;
  }
};

class WallCrossActionsFinder : public ActionFinderInterface {
 public:
  explicit WallCrossActionsFinder(float l) : l_{l, l, l} {};

  /// Find the next wall crossings for every particle before time t_max
  ActionList find_actions_in_cell(const ParticleList &plist,
                                  double t_max) const override;

  /// Ignore the neighbor searches for wall crossing
  ActionList find_actions_with_neighbors(const ParticleList &,
                                         const ParticleList &,
                                         double) const override {
    return {};
  }
  /// Ignore the surrounding searches for wall crossing
  ActionList find_actions_with_surrounding_particles(const ParticleList &,
                                                     const Particles &,
                                                     double) const override {
    return {};
  }

  /// No final actions for wall crossing
  ActionList find_final_actions(const Particles &, bool) const override {
    return {};
  }

 private:
  /// Periods in x,y,z directions in fm.
  const std::array<float, 3> l_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_WALLCROSSINGACTION_H_
