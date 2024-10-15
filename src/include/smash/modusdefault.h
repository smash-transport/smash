/*
 *
 *    Copyright (c) 2013-2022,2024
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_SMASH_MODUSDEFAULT_H_
#define SRC_INCLUDE_SMASH_MODUSDEFAULT_H_

#include <memory>

#include "configuration.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "grandcan_thermalizer.h"
#include "grid.h"
#include "outputinterface.h"
#include "potentials.h"

namespace smash {
/**
 * \ingroup modus
 * Base class for Modus classes that provides default function implementations.
 *
 * This is only a base class for actual Modus classes. Meaning there will never
 * be objects, references, or pointers to ModusDefault. Therefore, it does not
 * have - and will never need any virtual functions.
 *
 * The rules for adding functions to this class are as follows:
 * - This class is empty per default.
 * - You can add a function if you have a function that is different in at least
 *   one subclass.
 * - Code that is common to all goes into ExperimentImplementation.
 *
 * \todo JB: many of these functions could/should be virtual (contradicts
 * description given above, Vinzent says (in a nice way): does not make any
 * sense whatsoever anyway)
 */
class ModusDefault {
 public:
  // Never needs a virtual destructor.

  // Missing functions for concrete Modus implementations:
  // \todo(JB): Why is the function below commented out?
  // double initial_conditions(Particles *particles);

  /** Enforces sensible positions for the particles.
   *
   * Currently, this is only needed for BoxModus; the other Modi do
   * nothing.
   *
   * \see BoxModus::impose_boundary_conditions
   */
  int impose_boundary_conditions(Particles* /*p*/,
                                 const OutputsList& /*out_list*/ = {}) {
    return 0;
  }

  /// \return Checks if modus is collider; overwritten in ColliderModus
  bool is_collider() const { return false; }
  /// \return Checks if modus is a box; overwritten in BoxModus
  bool is_box() const { return false; }
  /// \return Checks if modus is list modus; overwritten in ListModus
  bool is_list() const { return false; }
  /// \return Checks if modus is sphere modus; overwritten in SphereModus
  bool is_sphere() const { return false; }
  /// \return Center of mass energy per nucleon pair in ColliderModus
  double sqrt_s_NN() const { return 0.; }
  /// \return The impact parameter; overwritten in ColliderModus
  double impact_parameter() const { return -1.; }
  /// sample impact parameter for collider modus
  void sample_impact() const {}
  /** \return The beam velocity of the projectile required in the Collider
   * modus. In the other modus, return zero. */
  double velocity_projectile() const { return 0.0; }
  /** \return The beam velocity of the target required in the Collider modus.
   * In the other modus, return zero. */
  double velocity_target() const { return 0.0; }
  /** \return The type of Fermi motion required in the Collider modus. In the
   * other modus, just return FermiMotion::Off. */
  FermiMotion fermi_motion() const { return FermiMotion::Off; }
  /// \return Maximal timestep accepted by this modus. Negative means infinity.
  double max_timestep(double) const { return -1.; }
  /// \return equilibration time of the box; overwritten in BoxModus
  double equilibration_time() const { return -1.; }
  /// \return length of the box; overwritten in BoxModus
  double length() const { return -1.; }
  /// \return radius of the sphere; overwritten in SphereModus
  double radius() const { return -1.; }
  /** \return Whether the calculation frame is fixed target;
   *  overwritten in ColliderModus */
  bool calculation_frame_is_fixed_target() const { return false; }
  /**
   * Get the passing time of the two nuclei in a collision. This time
   * corresponds to the moment when the nuclei have just passed entirely
   * through each other and all primary collisions have occured.
   * Formula taken from: Eq. (1) in \iref{Karpenko:2015xea}
   *
   * Only used in ColliderModus for IC output.
   * \return passing_time
   */
  double nuclei_passing_time() const { return 0.0; }
  /// \return Proper time of the hypersurface for IC in ColliderModus
  std::optional<double> proper_time() const { return std::nullopt; }
  /// \return Lower bound on proper time of the hypersurface for IC in
  /// ColliderModus
  std::optional<double> lower_bound() const { return std::nullopt; }
  /// \return Maximum rapidity for IC in ColliderModus
  std::optional<double> rapidity_cut() const { return std::nullopt; }
  /// \return Maximum transverse momentum for IC in ColliderModus
  std::optional<double> pT_cut() const { return std::nullopt; }
  /**
   * Creates the Grid with normal boundary conditions.
   *
   * \param[in] particles The Particles object containing all particles of the
   *                  currently running Experiment.
   * \param[in] min_cell_length The minimal length of the grid cells.
   * \param[in] timestep_duration Duration of the timestep. It is necessary for
   * formation times treatment: if particle is fully or partially formed before
   * the end of the timestep, it has to be on the grid.
   * \param[in] crit Collision criterion (decides if cell number can be limited)
   * \param[in] include_unformed_particles include unformed particles from
                                           the grid (worsens runtime, necessary
                                           for IC output)
   * \param[in] strategy The strategy to determine the cell size \return the
   * Grid object
   *
   * \see Grid::Grid
   */
  Grid<GridOptions::Normal> create_grid(
      const Particles& particles, double min_cell_length,
      double timestep_duration, CollisionCriterion crit,
      const bool include_unformed_particles,
      CellSizeStrategy strategy = CellSizeStrategy::Optimal) const {
    CellNumberLimitation limit = CellNumberLimitation::ParticleNumber;
    if (crit == CollisionCriterion::Stochastic) {
      limit = CellNumberLimitation::None;
    }
    return {particles,
            min_cell_length,
            timestep_duration,
            limit,
            include_unformed_particles,
            strategy};
  }

  /**
   * Creates GrandCanThermalizer
   *
   * \param[in] conf configuration object
   * \return unique pointer to created thermalizer class
   */
  std::unique_ptr<GrandCanThermalizer> create_grandcan_thermalizer(
      Configuration& conf) const {
    /* Lattice is placed such that the center is 0,0,0.
       If one wants to have a central cell with center at 0,0,0 then
       number of cells should be odd (2k+1) in every direction.
     */
    const std::array<double, 3> l =
        conf.take(InputKeys::forcedThermalization_latticeSizes);
    const std::array<double, 3> origin = {-0.5 * l[0], -0.5 * l[1],
                                          -0.5 * l[2]};
    const bool periodicity = false;
    return std::make_unique<GrandCanThermalizer>(conf, l, origin, periodicity);
  }

  /**
   * \ingroup exception
   *  BadInput is an error to throw if the configuration options are invalid.
   */
  struct BadInput : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
  /**
   * \ingroup exception
   * Thrown when the requested energy is smaller than the masses
   * of two particles.
   */
  struct InvalidEnergy : public BadInput {
    using BadInput::BadInput;
  };
};

}  // namespace smash

#endif  // SRC_INCLUDE_SMASH_MODUSDEFAULT_H_
