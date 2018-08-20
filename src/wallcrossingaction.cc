/*
 *
 *    Copyright (c) 2017-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#include "smash/wallcrossingaction.h"

namespace smash {

ActionList WallCrossActionsFinder::find_actions_in_cell(
    const ParticleList& plist, double t_max) const {
  std::vector<ActionPtr> actions;
  for (const ParticleData& p : plist) {
    const ThreeVector& r = p.position().threevec();
    const ThreeVector& v = p.velocity();
    double time_until_crossing = t_max;
    int i_cross = -1;
    for (int i = 0; i < 3; i++) {
      double t = t_max + 1.;
      if (v[i] > really_small) {
        t = (l_[i] - r[i]) / v[i];
      } else if (v[i] < -really_small) {
        t = -r[i] / v[i];
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
                              r + v * time_until_crossing);
    crossing_point[i_cross + 1] = ((v[i_cross] > 0.0) ? 0.0 : l_[i_cross]);

    ParticleData outgoing_particle(p);
    outgoing_particle.set_4position(crossing_point);
    ActionPtr action = make_unique<WallcrossingAction>(p, outgoing_particle,
                                                       time_until_crossing);
    actions.emplace_back(std::move(action));
  }
  return actions;
}

}  // namespace smash
