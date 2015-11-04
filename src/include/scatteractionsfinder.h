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

#include <vector>

#include "action.h"
#include "actionfinderfactory.h"
#include "configuration.h"
#include "constants.h"

namespace Smash {

/**
 * \ingroup action
 * A simple scatter finder:
 * Just loops through all particles and checks each pair for a collision.  */
class ScatterActionsFinder : public ActionFinderInterface {
 public:
  /** Initialize the finder with the given parameters. */
  ScatterActionsFinder(Configuration config,
                       const ExperimentParameters &parameters,
                       bool two_to_one, bool two_to_two);
  /** Constructor for testing purposes. */
  ScatterActionsFinder(float elastic_parameter, int testparticles);

  /** Determine the collision time of the two particles [fm/c].
   *  Time of the closest approach is taken as collision time.
   *
   * \fpPrecision Why \c double?
   */
  static inline double collision_time(const ParticleData &p1,
                                      const ParticleData &p2) {
    /** UrQMD collision time in computational frame,
    * see \iref{Bass:1998ca} (3.28):
    * position of particle 1: r_1 [fm]
    * position of particle 2: r_2 [fm]
    * velocity of particle 1: v_1
    * velocity of particle 1: v_2
    * t_{coll} = - (r_1 - r_2) . (v_1 - v_2) / (v_1 - v_2)^2 [fm/c]
    */
    const ThreeVector dv_times_e1e2 =
            p1.momentum().threevec() * p2.momentum().x0() -
            p2.momentum().threevec() * p1.momentum().x0();
    const double dv_times_e1e2_sqr = dv_times_e1e2.sqr();
    /* Zero relative velocity . particles are not approaching. */
    if (dv_times_e1e2_sqr < really_small) {
      return -1.0;
    }
    const ThreeVector dr = p1.position().threevec() - p2.position().threevec();
    return -(dr*dv_times_e1e2) *
             (p1.momentum().x0() * p2.momentum().x0() / dv_times_e1e2_sqr);
  }
  /** Check the whole particle list for collisions
   * and return a list with the corrsponding Action objects. */
  ActionList find_actions_in_cell(const ParticleList &search_list,
                                  float dt) const override;
  ActionList find_actions_with_neighbors(const ParticleList &search_list,
                                         const ParticleList &neighbors_list,
                                         float dt) const override;
  ActionList find_actions_with_surrounding_particles(
      const ParticleList &search_list, const Particles &surrounding_list,
      float dt) const override;
  /** Find some final collisions at the end of the simulation.
   * Currently does nothing. */
  ActionList find_final_actions(const Particles &) const override {
    return ActionList();
  }

  /**
   * If there is only one particle sort, no decays
   * (only elastic scatterings are possible),
   * scatterings are isotropic and cros-section fixed to elastic_parameter_
   * independently on momenta, then maximal cross-section is elastic_parameter_.
   * This knowledge can be used for improving performance.
   */
  inline bool is_constant_elastic_isotropic() const {
    return ParticleType::list_all().size() == 1 &&
           !two_to_one_ &&
           isotropic_ &&
           elastic_parameter_ > 0.0f;
  }

  /**
   * Returns the maximal transverse distance squared.
   *
   * Particle pairs whose transverse distance is larger then this, are not
   * checked for collisions.
   */
  float max_transverse_distance_sqr(int testparticles) const {
    return (is_constant_elastic_isotropic() ? elastic_parameter_ :
            maximum_cross_section) / testparticles * fm2_mb * M_1_PI;
  }

 private:
  /* Construct a ScatterAction object,
   * based on the types of the incoming particles. */
  ScatterActionPtr construct_scatter_action(const ParticleData &data_a,
                                            const ParticleData &data_b,
                                            float time_until_collision) const;
  /** Check for a single pair of particles (id_a, id_b) if a collision will happen
   * in the next timestep and create a corresponding Action object in that case. */
  ActionPtr check_collision(const ParticleData &data_a,
                            const ParticleData &data_b, float dt) const;
  /** Elastic cross section parameter (in mb). */
  const float elastic_parameter_;
  /** Number of test particles. */
  const int testparticles_;
  /** Do all collisions isotropically. */
  const bool isotropic_;
  /** Enable 2->1 processes. */
  const bool two_to_one_;
  /** Enable 2->2 processes. */
  const bool two_to_two_;
};

#if 0
/* An advanced scatter finder:
 * Sets up a grid and sorts the particles into grid cells. */
class GridScatterFinder : public ScatterActionsFinder {
 public:
  explicit GridScatterFinder(float length);
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
