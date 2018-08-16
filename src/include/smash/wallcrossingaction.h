/*
 *
 *    Copyright (c) 2016-2018
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_WALLCROSSINGACTION_H_
#define SRC_INCLUDE_WALLCROSSINGACTION_H_

#include "action.h"
#include "actionfinderfactory.h"

namespace smash {

/**
 * \ingroup action
 * WallcrossingAction is a special action which indicates that a particle
 * has crossed a box wall.
 */
class WallcrossingAction : public Action {
 public:
  /**
   * Construct wallcrossing action.
   * \param[in] in_part Data of incoming particle.
   * \param[in] out_part Data of outgoing particle. Same as in_part, but with
   * updated position.
   * \param[in] time_until Time when the crossing takes place (relative to
   * current time). Usually no delay therefore optional. [fm]
   */
  WallcrossingAction(const ParticleData &in_part, const ParticleData &out_part,
                     const double time_until = 0.0)
      : Action(in_part, out_part, time_until, ProcessType::Wall) {}
  double get_total_weight() const override { return 0.0; };
  double get_partial_weight() const override { return 0.0; };
  void generate_final_state() override{};
  void format_debug_output(std::ostream &out) const override {
    out << "Wall crossing of " << incoming_particles_;
  }
};

/**
 * \ingroup action
 * Finder for wall crossing actions, when using peridic boundary conditons.
 * Loops through all particles and checks if they cross the box wall during the
 * next timestep.
 */
class WallCrossActionsFinder : public ActionFinderInterface {
 public:
  /**
   * Construct wallcrossing actionfinder.
   * \param[in] l Box edge length. Box is assumbed to be a cube. [fm]
   */
  explicit WallCrossActionsFinder(double l) : l_{l, l, l} {};

  /**
   * Find the next wall crossings for every particle before time t_max.
   * \param[in] plist List of all particles.
   * \param[in] t_max Time until crossing can appear. [fm]
   * \return List of all found wall crossings.
   */
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
  const std::array<double, 3> l_;
};

}  // namespace smash

#endif  // SRC_INCLUDE_WALLCROSSINGACTION_H_
