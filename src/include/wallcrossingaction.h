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
  float raw_weight_value() const override { return 1; };
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
  WallCrossActionsFinder(float l) : l_{l, l, l} {};

  ActionList find_actions_in_cell(const ParticleList &plist,
                                  float t_max) const override {
    std::vector<ActionPtr> actions;
    for (const ParticleData &p : plist) {
      const ThreeVector r = p.position().threevec();
      const ThreeVector v = p.velocity();
      double time_until_crossing = t_max;
      int i_cross = -1;
      for (int i = 0; i < 3; i++) {
        double t = t_max + 1.0f;
        if (v[i] > really_small) {
          t = (l_[i] - r[i])/v[i];
        } else if (v[i] < -really_small) {
          t = -r[i]/v[i];
        }
        if (t < time_until_crossing) {
          time_until_crossing = t;
          i_cross = i;
        }
      }
      // No crossing
      if (i_cross == -1) {
        continue;
      }
      FourVector crossing_point(p.position().x0() + time_until_crossing,
                                r + v*time_until_crossing);
      crossing_point[i_cross + 1] = ((v[i_cross] > 0.0) ? 0.0 : l_[i_cross]);

      ParticleData outgoing_particle(p);
      outgoing_particle.set_4position(crossing_point);
      ActionPtr action = make_unique<WallcrossingAction>(p, outgoing_particle,
                                                         time_until_crossing);
      actions.push_back(std::move(action));
    }
    return actions;
  }

  /// Ignore the neighbor searches for wall crossing
  ActionList find_actions_with_neighbors(const ParticleList &,
                                         const ParticleList &,
                                         float) const override {
    return {};
  }
  /// Ignore the surrounding searches for wall crossing
  ActionList find_actions_with_surrounding_particles(const ParticleList &,
                                                     const Particles &,
                                                     float) const override {
    return {};
  }

  /// No final actions for wall crossing
  ActionList find_final_actions(const Particles &, bool) const override {
    return {};
  }

 private:
  const std::array<float, 3> l_;
};

}  // namespace Smash

#endif  // SRC_INCLUDE_WALLCROSSINGACTION_H_
