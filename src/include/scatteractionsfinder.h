/*
 *
 *    Copyright (c) 2014-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SCATTERACTIONSFINDER_H_
#define SRC_INCLUDE_SCATTERACTIONSFINDER_H_

#include "actionfinderfactory.h"
#include "action.h"
#include "configuration.h"

namespace Smash {

/**
 * \ingroup action
 * A simple scatter finder:
 * Just loops through all particles and checks each pair for a collision.  */
class ScatterActionsFinder : public ActionFinderInterface {
 public:
  /** Initialize the finder with the given parameters. */
  ScatterActionsFinder(Configuration config,
                       const ExperimentParameters &parameters);
  ScatterActionsFinder(float elastic_parameter, int testparticles);

  /** Determine the collision time of the two particles. */
  static double collision_time(const ParticleData &p_a,
                               const ParticleData &p_b);
  /** Check the whole particle list for collisions
   * and return a list with the corrsponding Action objects. */
  ActionList find_possible_actions(
      const ParticleList &search_list,
      const std::vector<const ParticleList *> &neighbors_list,
      float dt) const override;
  /** Find some final collisions at the end of the simulation.
   * Currently does nothing. */
  ActionList find_final_actions(const ParticleList &) const override
    { return ActionList(); }

 private:
  /** Check for a single pair of particles (id_a, id_b) if a collision will happen
   * in the next timestep and create a corresponding Action object in that case. */
  ActionPtr check_collision(const ParticleData &data_a,
                            const ParticleData &data_b, float dt) const;
  /** Elastic cross section parameter (in mb). */
  float elastic_parameter_ = 0.0;
  /** Number of test particles. */
  int testparticles_ = 1;
};

#if 0
/* An advanced scatter finder:
 * Sets up a grid and sorts the particles into grid cells. */
class GridScatterFinder : public ScatterActionsFinder {
 public:
  GridScatterFinder(float length);
  void find_possible_actions(std::vector<ActionPtr> &actions,
                             Particles *particles,
                             const ExperimentParameters &parameters,
                             CrossSections *cross_sections = nullptr)
                             const override;
 private:
  /* Cube edge length. */
  const float length_;
};
#endif

}  // namespace Smash

#endif  // SRC_INCLUDE_SCATTERACTIONSFINDER_H_
