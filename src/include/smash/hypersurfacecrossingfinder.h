/*
 *
 *    Copyright (c) 2019-2020,2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_
#define SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_

#include <vector>

#include "action.h"
#include "actionfinderfactory.h"

namespace smash {

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
   * \param[in] y Value for rapidity cut: absolute value of momentum space
   *            rapidity up to which particles are considered for initial
   *            conditions
   * \param[in] pT Value for transverse momentum cut: maximum transverse
   *            momentum up to which particles are considered for initial
   *            conditions
   */
  explicit HyperSurfaceCrossActionsFinder(double tau, double y, double pT)
      : prop_time_{tau}, rap_cut_{y}, pT_cut_{pT} {};

  /**
   * Find the next hypersurface crossings for each particle that occur within
   * the timestepless propagation.
   * \param[in] plist List of all particles.
   * \param[in] dt Time until crossing can appear (until end of timestep). [fm]
   * \param[in] beam_momentum [GeV] List of beam momenta for each particle;
   * only necessary for frozen Fermi motion
   * necessary if frozen Fermi Motion is activated \return List of all found
   * wall crossings.
   */
  ActionList find_actions_in_cell(
      const ParticleList &plist, double dt, const double,
      const std::vector<FourVector> &beam_momentum) const override;

  /// Ignore the neighbor searches for hypersurface crossing
  ActionList find_actions_with_neighbors(
      const ParticleList &, const ParticleList &, double,
      const std::vector<FourVector> &) const override {
    return {};
  }

  /// Ignore the surrounding searches for hypersurface crossing
  ActionList find_actions_with_surrounding_particles(
      const ParticleList &, const Particles &, double,
      const std::vector<FourVector> &) const override {
    return {};
  }

  /// No final actions for hypersurface crossing
  ActionList find_final_actions(const Particles &, bool) const override {
    return {};
  }

 private:
  /// Proper time of the hypersurface in fm.
  const double prop_time_;

  /**
   * Rapidity (momentum space) cut for the particles contributing to the initial
   * conditions for hydrodynamics.
   * If applied, only particles characterized by a rapidity between
   * [-y_cut, y_cut] are printed to the hypersurface.
   */
  const double rap_cut_;

  /**
   * Transverse momentum cut for the particles contributing to the initial
   * conditions for hydrodynamics.
   * If applied, only particles characterized by a transverse momentum between
   * [0, pT_cut] are
   * printed to the hypersurface.
   */
  const double pT_cut_;

  /**
   * Determine whether particle crosses hypersurface within next timestep
   * during propagation
   * \param[in] pdata_before_propagation Particle data at the beginning of time
   *            step in question
   * \param[in] pdata_after_propagation Particle data at the end of time step
   *            in question
   * \param[in] tau Proper time of the hypersurface that is tested
   * \return Does particle cross the hypersurface?
   */
  bool crosses_hypersurface(ParticleData &pdata_before_propagation,
                            ParticleData &pdata_after_propagation,
                            const double tau) const;

  /**
   * Find the coordinates where particle crosses hypersurface
   * \param[in] pdata_before_propagation Particle data at the beginning of time
   *            in question
   * \param[in] pdata_after_propagation Particle data at the end of time step
   *            in question
   * \param[in] tau Proper time of the hypersurface that is crossed
   * \return Fourvector of the crossing position
   */
  FourVector coordinates_on_hypersurface(ParticleData &pdata_before_propagation,
                                         ParticleData &pdata_after_propagation,
                                         const double tau) const;
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_HYPERSURFACECROSSINGACTION_H_
