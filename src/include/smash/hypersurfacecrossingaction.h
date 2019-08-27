/*
 *
 *    Copyright (c) 2016-2019
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_HYPERSURFACECROSSINGACTION_H_
#define SRC_INCLUDE_HYPERSURFACECROSSINGACTION_H_

#include "action.h"
#include "actionfinderfactory.h"

namespace smash {

/**
 * \ingroup action
 * Hypersurfacecrossingaction is a special action which indicates that a
 * particle has crossed a hypersurface of given proper time. This is necessary
 * to generate initial conditions for hybrid models.
 */
class HypersurfacecrossingAction : public Action {
 public:
  /**
   * Construct hypersurfacecrossing action.
   * \param[in] in_part Data of incoming particle.
   * \param[in] out_part Data of particles leaving hypersurface
   * \param[in] time_until Time when the crossing takes place[fm]
   */
  HypersurfacecrossingAction(const ParticleData &in_part,
                             const ParticleData &out_part,
                             const double time_until)
      : Action(in_part, out_part, time_until,
               ProcessType::HyperSurfaceCrossing) {}
  double get_total_weight() const override { return 0.0; };
  double get_partial_weight() const override { return 0.0; };
  void format_debug_output(std::ostream &out) const override {
    out << "Hypersurface crossing of " << incoming_particles_;
  }

  /**
   * Generate the final state of the hypersurface crossing particles.
   * Removes all particles from the evolution and redirects them to the output.
   *
   */
  void generate_final_state() override;
};

/**
 * \ingroup action
 * Finder for hypersurface crossing actions.
 * Loops through all particles and checks if they cross the hypersurface
 * during the next timestep.
 */
class HyperSurfaceCrossActionsFinder : public ActionFinderInterface {
 public:
  /**
   * Construct hypersurfacecrossing action finder.
   * \param[in] tau Proper time of the hypersurface. [fm]
   */
  explicit HyperSurfaceCrossActionsFinder(double tau) : prop_time_{tau} {};

  /**
   * Find the next hypersurface crossings for every particle before time t_max.
   * \param[in] plist List of all particles.
   * \param[in] t_max Time until crossing can appear. [fm]
   * \return List of all found wall crossings.
   */
  ActionList find_actions_in_cell(const ParticleList &plist,
                                  double t_max) const override;

  /// Ignore the neighbor searches for hypersurface crossing
  ActionList find_actions_with_neighbors(const ParticleList &,
                                         const ParticleList &,
                                         double) const override {
    return {};
  }

  /// Ignore the surrounding searches for hypersurface crossing
  ActionList find_actions_with_surrounding_particles(const ParticleList &,
                                                     const Particles &,
                                                     double) const override {
    return {};
  }

  /// No final actions for hypersurface crossing
  ActionList find_final_actions(const Particles &, bool) const override {
    return {};
  }

 private:
  /// Proper time of the hypersurface in fm.
  const double prop_time_;

  bool crosses_hypersurface(ParticleData &pdata_before_propagation,
                            ParticleData &pdata_after_propagation,
                            const double tau) const;
  FourVector coordinates_on_hypersurface(ParticleData &pdata_before_propagation,
                                         ParticleData &pdata_after_propagation,
                                         const double tau) const;
};

}  // namespace smash

#endif  // SRC_INCLUDE_HYPERSURFACECROSSINGACTION_H_
