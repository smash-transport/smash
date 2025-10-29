/*
 *
 *    Copyright (c) 2014-2022
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_SCATTERACTIONSFINDER_H_
#define SRC_INCLUDE_SMASH_SCATTERACTIONSFINDER_H_

#include <memory>
#include <set>
#include <vector>

#include "action.h"
#include "actionfinderfactory.h"
#include "configuration.h"
#include "scatteraction.h"
#include "scatteractionsfinderparameters.h"

namespace smash {

/**
 * \ingroup action
 * A simple scatter finder:
 * Just loops through all particles and checks each pair for a collision.
 * It supports two collision criteria: a geometric and stochastic criterion.
 */
class ScatterActionsFinder : public ActionFinderInterface {
 public:
  /**
   * Constructor of the finder with the given parameters.
   *
   * \param[inout] config Configuration of smash from which we take:
   *               1) A global elastic cross section [mb]. It will be used
   *                  regardless of the species of the colliding particles.
   *                  It won't be used if the value is negative.
   *               2) An option determining whether all the scatterings are
   *                  isotropic
   *               3) Parameters of the string process
   * \param[in] parameters Struct of parameters determining whether to
   *            exclude some certain types of scatterings and switching
   *            among the methods to treat with the NNbar collisions.
   */
  ScatterActionsFinder(Configuration &config,
                       const ExperimentParameters &parameters);

  /**
   * Determine the collision time of the two particles.
   * Time of the closest approach is taken as collision time, if the geometric
   * collision criterion is used. For stochastic criterion the time is
   * distributed uniformly within the timestep.
   *
   * \param[in] p1 First incoming particle
   * \param[in] p2 Second incoming particle
   * \param[in] dt The maximum time interval at the current time step [fm]
   * \param[in] beam_momentum [GeV] List of beam momenta for each particle;
   * only necessary for frozen Fermi motion
   * \return How long does it take for the two incoming particles
   *         to propagate before scattering [fm]. It's set equal
   *         to -1 if the two particles are not moving relative to each
   *         other.
   */
  inline double collision_time(
      const ParticleData &p1, const ParticleData &p2, double dt,
      const std::vector<FourVector> &beam_momentum) const {
    if (finder_parameters_.coll_crit == CollisionCriterion::Stochastic) {
      return dt * random::uniform(0., 1.);
    } else {
      /*
       * For frozen Fermi motion:
       * If particles have not yet interacted and are the initial nucleons,
       * perform action finding with beam momentum instead of Fermi motion
       * corrected momentum. That is because the particles are propagated with
       * the beam momentum until they interact.
       */
      if (p1.id() < 0 || p2.id() < 0) {
        throw std::runtime_error("Invalid particle ID for Fermi motion");
      }
      const bool p1_has_no_prior_interactions =
          (static_cast<uint64_t>(p1.id()) <                 // particle from
           static_cast<uint64_t>(beam_momentum.size())) &&  // initial nucleus
          (p1.get_history().collisions_per_particle == 0);

      const bool p2_has_no_prior_interactions =
          (static_cast<uint64_t>(p2.id()) <                 // particle from
           static_cast<uint64_t>(beam_momentum.size())) &&  // initial nucleus
          (p2.get_history().collisions_per_particle == 0);

      const FourVector p1_mom = (p1_has_no_prior_interactions)
                                    ? beam_momentum[p1.id()]
                                    : p1.momentum();
      const FourVector p2_mom = (p2_has_no_prior_interactions)
                                    ? beam_momentum[p2.id()]
                                    : p2.momentum();
      if (finder_parameters_.coll_crit == CollisionCriterion::Covariant) {
        /**
         * JAM collision times from the closest approach
         * in the two-particle center-of-mass-framem,
         * see \iref{Hirano:2012yy} (5.13) and (5.14).
         * The scatteraction is performed at the mean of these two times.
         */
        const FourVector delta_x = p1.position() - p2.position();
        const double p1_sqr = p1_mom.sqr();
        const double p2_sqr = p2_mom.sqr();
        const double p1_dot_x = p1_mom.Dot(delta_x);
        const double p2_dot_x = p2_mom.Dot(delta_x);
        const double p1_dot_p2 = p1_mom.Dot(p2_mom);
        const double denominator = std::pow(p1_dot_p2, 2) - p1_sqr * p2_sqr;
        if (unlikely(std::abs(denominator) < really_small * really_small)) {
          return -1.0;
        }

        const double time_1 = (p2_sqr * p1_dot_x - p1_dot_p2 * p2_dot_x) *
                              p1_mom.x0() / denominator;
        const double time_2 = -(p1_sqr * p2_dot_x - p1_dot_p2 * p1_dot_x) *
                              p2_mom.x0() / denominator;
        return (time_1 + time_2) / 2;
      } else {
        /**
         * UrQMD collision time in computational frame,
         * see \iref{Bass:1998ca} (3.28):
         * position of particle 1: \f$r_1\f$ [fm]
         * position of particle 2: \f$r_2\f$ [fm]
         * velocity of particle 1: \f$v_1\f$
         * velocity of particle 1: \f$v_2\f$
         * \f[t_{coll} = - (r_1 - r_2) . (v_1 - v_2) / (v_1 - v_2)^2\f] [fm]
         */
        const ThreeVector dv_times_e1e2 =
            p1_mom.threevec() * p2_mom.x0() - p2_mom.threevec() * p1_mom.x0();
        const double dv_times_e1e2_sqr = dv_times_e1e2.sqr();
        if (dv_times_e1e2_sqr < really_small) {
          return -1.0;
        }
        const ThreeVector dr =
            p1.position().threevec() - p2.position().threevec();
        return -(dr * dv_times_e1e2) *
               (p1_mom.x0() * p2_mom.x0() / dv_times_e1e2_sqr);
      }
    }
  }

  /**
   * Search for all the possible collisions within one cell. This function is
   * only used for counting the primary collisions at the beginning of each
   * time step. (Although it's also called afterwards for searching the
   * secondary collisions among the outgoing particles, no new actions will be
   * found since the scattered pairs cannot scatter again.)
   *
   * \param[in] search_list A list of particles within one cell
   * \param[in] dt The maximum time interval at the current time step [fm]
   * \param[in] gcell_vol Volume of searched grid cell [fm^3]
   * \param[in] beam_momentum [GeV] List of beam momenta for each particle;
   * only necessary for frozen Fermi motion
   * \return A list of possible scatter actions
   */
  ActionList find_actions_in_cell(
      const ParticleList &search_list, double dt, const double gcell_vol,
      const std::vector<FourVector> &beam_momentum) const override;

  /**
   * Search for all the possible collisions among the neighboring cells. This
   * function is only used for counting the primary collisions at the beginning
   * of each time step.
   *
   * \param[in] search_list A list of particles within the current cell
   * \param[in] neighbors_list A list of particles within the neighboring cell
   * \param[in] dt The maximum time interval at the current time step [fm]
   * \param[in] beam_momentum [GeV] List of beam momenta for each particle;
   * only necessary for frozen Fermi motion
   * \return A list of possible scatter actions
   */
  ActionList find_actions_with_neighbors(
      const ParticleList &search_list, const ParticleList &neighbors_list,
      double dt, const std::vector<FourVector> &beam_momentum) const override;

  /**
   * Search for all the possible secondary collisions between the outgoing
   * particles and the rest.
   *
   * \param[in] search_list A list of particles within the current cell
   * \param[in] surrounding_list The whole particle list
   * \param[in] dt The maximum time interval at the current time step [fm]
   * \param[in] beam_momentum [GeV] List of beam momenta for each particle;
   * only necessary for frozen Fermi motion
   * \return A list of possible scatter actions
   */
  ActionList find_actions_with_surrounding_particles(
      const ParticleList &search_list, const Particles &surrounding_list,
      double dt, const std::vector<FourVector> &beam_momentum) const override;

  /// No scatterings should be found when the event is over.
  ActionList find_final_actions(const Particles &) const override { return {}; }

  /**
   * If there is only one particle sort, no decays
   * (only elastic scatterings are possible),
   * scatterings are isotropic and cross-section fixed to elastic_parameter_
   * independently on momenta, then maximal cross-section is elastic_parameter_.
   * This knowledge can be used for improving performance.
   *
   * \return A boolean indicating whether all the scatterings are elastic
   *         and isotropic
   */
  inline bool is_constant_elastic_isotropic() const {
    return ParticleType::list_all().size() == 1 &&
           !finder_parameters_.two_to_one && isotropic_ &&
           finder_parameters_.elastic_parameter > 0.;
  }

  /**
   * The maximal distance over which particles can interact in case of the
   * geometric criterion, related to the number of test particles and the
   * maximal cross section.
   *
   * \param[in] testparticles Number of test particles.
   *
   * \return Maximal transverse distance squared. [fm\f$^{2}\f$]
   *         Particle pairs whose transverse distance is larger than this
   *         are not checked for collisions in case of the default geometric
   *         collision criterion.
   */
  double max_transverse_distance_sqr(int testparticles) const {
    return (is_constant_elastic_isotropic()
                ? finder_parameters_.elastic_parameter
                : finder_parameters_.maximum_cross_section) /
           testparticles * fm2_mb * M_1_PI;
  }

  /**
   * Prints out all the 2-> n (n > 1) reactions with non-zero cross-sections
   * between all possible pairs of particle types.
   */
  void dump_reactions() const;

  /**
   * Print out partial cross-sections of all processes that can occur in
   * the collision of a(mass = m_a) and b(mass = m_b).
   *
   * \param[in] a The specie of the first incoming particle.
   * \param[in] b The specie of the second incoming particle.
   * \param[in] m_a Mass of species a [GeV].
   * \param[in] m_b Mass of species b [GeV].
   * \param[in] final_state Whether the final state cross sections should be
   *                        printed.
   * \param[in] plab Optional momenta in lab frame to be evaluated [GeV].
   *                 Ignored if empty.
   */
  void dump_cross_sections(const ParticleType &a, const ParticleType &b,
                           double m_a, double m_b, bool final_state,
                           std::vector<double> &plab) const;

  /**
   * \return Pointer to the string process class object.
   *         If string is turned off, the null pointer is returned.
   */
  StringProcess *get_process_string_ptr() {
    if (finder_parameters_.strings_switch) {
      return string_process_interface_.get();
    } else {
      return NULL;
    }
  }

 private:
  /**
   * Check for a single pair of particles (id_a, id_b) if a collision will
   * happen in the next timestep and create a corresponding Action object
   * in that case.
   *
   * Two criteria for the collision decision are supported: 1. The default
   * geometric criterion from UrQMD \iref{Bass:1998ca} (3.27). 2. A stochastic
   * collision criterion as introduced in \iref{Staudenmaier:2021lrg}.
   *
   * \param[in] data_a First incoming particle
   * \param[in] data_b Second incoming particle
   * \param[in] dt Maximum time interval within which a collision can happen
   * \param[in] beam_momentum [GeV] List of beam momenta for each particle;
   * only necessary for frozen Fermi motion
   * \param[in] gcell_vol (optional) volume of grid cell in which the collision
   *                                is checked
   * \return A null pointer if no collision happens or an action which contains
   *         the information of the outgoing particles.
   *
   * Note: gcell_vol is optional, since only find_actions_in_cell has (and
   * needs) this information for the stochastic collision criterion.
   */
  ActionPtr check_collision_two_part(
      const ParticleData &data_a, const ParticleData &data_b, double dt,
      const std::vector<FourVector> &beam_momentum = {},
      const double gcell_vol = 0.0) const;

  /**
   * Check for multiple i.e. more than 2 particles if a collision will happen in
   * the next timestep and create a corresponding Action object in that case.
   *
   * This is only possible for the stochastic collision criterion, which is
   * introduced in \iref{Staudenmaier:2021lrg}. Following the same general idea
   * as for the 2-particle scatterings, probabilities for multi-particle
   * scatterings can be derived.
   *
   * \param[in] plist List of incoming particles
   * \param[in] dt Maximum time interval within which a collision can happen
   * \param[in] gcell_vol volume of grid cell in which the collision is checked
   * \return A null pointer if no collision happens or an action which contains
   *         the information of the outgoing particles.
   */
  ActionPtr check_collision_multi_part(const ParticleList &plist, double dt,
                                       const double gcell_vol) const;

  /// Struct collecting several parameters.
  ScatterActionsFinderParameters finder_parameters_;
  /// Class that deals with strings, interfacing Pythia.
  std::unique_ptr<StringProcess> string_process_interface_;
  /// Do all collisions isotropically.
  const bool isotropic_;
  /**
   * Box length: needed to determine coordinates of collision
   * correctly in case of collision through the wall.
   * Ignored if negative.
   */
  const double box_length_;
  /// Parameter for formation time
  const double string_formation_time_;
};

/**
 * Gather all relevant parameters for a \c ScatterActionsFinder either getting
 * them from an \c ExperimentParameters  instance or extracting them from
 * a \c Configuration .
 *
 * \param[in] parameters The parameters of the considered experiment
 * \param[inout] config SMASH input configuration
 * \return A \c ScatterActionsFinderParameters appropriately filled.
 */
ScatterActionsFinderParameters create_finder_parameters(
    Configuration &config, const ExperimentParameters &parameters);

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_SCATTERACTIONSFINDER_H_
