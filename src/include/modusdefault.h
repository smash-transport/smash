/*
 *
 *    Copyright (c) 2013-2015
 *      SMASH Team
 *
 *    GNU General Public License (GPLv3 or later)
 *
 */

#ifndef SRC_INCLUDE_MODUSDEFAULT_H_
#define SRC_INCLUDE_MODUSDEFAULT_H_

#include "configuration.h"
#include "cxx14compat.h"
#include "forwarddeclarations.h"
#include "fourvector.h"
#include "grandcan_thermalizer.h"
#include "grid.h"
#include "outputinterface.h"
#include "potentials.h"

namespace Smash {

/**
 * \ingroup modus
 * Baseclass for Modus classes that provides default function implementations.
 *
 * This is only a base class for actual Modus classes. Meaning there will never
 * be objects, references, or pointers to ModusDefault. Therefore, it does not
 * have - and will never need any virtual functions.
 *
 * The rules for adding functions to this class are as follows:
 * - This class is empty per default.
 * - You can add a function if you have a function that is different in at least
 *   two subclasses and equal in at least two subclasses.
 * - Code that is common to all goes into ExperimentImplementation.
 */
class ModusDefault {
 public:
  // never needs a virtual destructor

  // Missing functions for concrete Modus implementations:
  // float initial_conditions(Particles *particles);

  /** Enforces sensible positions for the particles.
   *
   * Currently, this is only needed for BoxModus; the other Modi do
   * nothing.
   *
   * \see BoxModus::impose_boundary_conditions
   */
  int impose_boundary_conditions(Particles * /*p*/,
                         const OutputsList & /*out_list*/ = {})
  { return 0;}

  /** The following four parameters tell the number of nucleons in the
   *  colliding nuclei, the number of nucleons in the projectile, whether
   *  the first collisions within the same nucleus are forbiden, and whether
   *  the modus is Collider. They are all needed in the Collider modus. In
   *  the other modus, they take the following default value. */
  int total_N_number() const { return 0; }
  int proj_N_number() const {return 0; }
  bool cll_in_nucleus() const { return true; }
  bool is_collider() const { return false; }
  /** The beam velocity of the projectile required in the Collider modus.
   *   In the other modus, return zero. */
  double velocity_projectile() const { return 0.0; }
  /** The beam velocity of the target required in the Collider modus.
   *   In the other modus, return zero. */
  double velocity_target() const { return 0.0; }
  /** The type of Fermi motion required in the Collider modus. In the other
   *  modus, just return FermiMotion::Off. */
  FermiMotion fermi_motion() const { return FermiMotion::Off; }
  /// Maximal timestep accepted by this modus. Negative means infinity.
  float max_timestep(float ) const { return -1.f; }

  float length() const { return -1.f; }

  /**
   * Creates the Grid with normal boundary conditions.
   *
   * \param particles The Particles object containing all particles of the
   *                  currently running Experiment.
   * \param min_cell_length The minimal length of the grid cells.
   * \param strategy The strategy to determine the cell size
   *
   * \see Grid::Grid
   */
  Grid<GridOptions::Normal> create_grid(
      const Particles &particles, float min_cell_length,
      CellSizeStrategy strategy = CellSizeStrategy::Optimal) const {
    return {particles, min_cell_length, strategy};
  }

  /**
   * Creates GrandCanThermalizer
   *
   * \param[in] conf configuration object
   * \param[out] th unique pointer to created thermalizer class
   */
  std::unique_ptr<GrandCanThermalizer> create_grandcan_thermalizer(
                                               Configuration& conf) const {
    /* Lattice is placed such that the center is 0,0,0.
       If one wants to have a central cell with center at 0,0,0 then
       number of cells should be odd (2k+1) in each direction.
     */
    const std::array<float, 3> l = conf.take({"Lattice_Sizes"});
    const std::array<float, 3> origin = {-0.5f*l[0], -0.5f*l[1], -0.5f*l[2]};
    const bool periodicity = false;
    return make_unique<GrandCanThermalizer>(conf, l, origin, periodicity);
  }

  /** \ingroup exception
   *  BadInput is an error to throw if the configuration options are invalid.
   *
   **/
  struct BadInput : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };
  /// \ingroup exception
  /// Thrown when the requested energy is smaller than the masses
  /// of two particles.
  struct InvalidEnergy : public BadInput {
    using BadInput::BadInput;
  };
};

}  // namespace Smash

#endif  // SRC_INCLUDE_MODUSDEFAULT_H_
