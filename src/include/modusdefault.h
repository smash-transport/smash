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

#include <stdexcept>

#include "forwarddeclarations.h"
#include "grid.h"

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

  /**
   * Creates the Grid with normal boundary conditions.
   *
   * \param particles The Particles object containing all particles of the
   *                  currently running Experiment.
   * \param testparticles The number of testparticles.
   * \param strategy The strategy to determine the cell size
   *
   * \see Grid::Grid
   */
  Grid<GridOptions::Normal> create_grid(
      const Particles &particles, int testparticles,
      CellSizeStrategy strategy = CellSizeStrategy::Optimal) const {
    return {particles, testparticles, strategy};
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
